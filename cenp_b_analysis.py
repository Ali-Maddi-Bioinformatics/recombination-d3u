"""Independent, streaming CENP-B motif analysis for T2T/hs1.

Run: python cenp_b_analysis.py
Help: python cenp_b_analysis.py --help

Only this module's two CSV outputs are replaced after a successful scan.
No imports from main.py and no network requests are made.

Output: Genome_Details/CENP_B_sites.csv and CENP_B_chromosome_summary.csv.
Both tables use 1-based inclusive coordinates: End_Locus - Start_Locus + 1
is always 17. This differs from the legacy colonies.py end convention.
Sequence is the reference (+) sequence, including for '-' strand sites;
Motif_Oriented_Sequence is reverse-complemented for those sites.
Distances count intervening bases between consecutive loci on either strand,
within each chromosome; negative values denote overlaps. First/last gaps
are blank. Boundary-overlapping sites count as outside (not fully inside)
and are also counted separately. Missing acen annotations (chrM) produce
blank comparison counts, not a claim of absence or containment.

UCSC cytobands use 0-based half-open coordinates and are converted here:
https://genome-blog.gi.ucsc.edu/blog/2016/12/12/the-ucsc-genome-browser-coordinate-counting-systems/
The requested broad recognition motif is described in:
https://www.nature.com/articles/srep27833
Sequence matching alone does not demonstrate CENP-B binding. Cytoband acen
annotations are the comparison intervals, not experimentally mapped binding.

Each input file must contain one chromosome FASTA record. Soft masking is
ignored; X matches A/C/G/T only. IUPAC ambiguity symbols consume coordinates
but cannot form a hit. At most --chunk-mb MiB plus a 16 bp overlap are read
per worker at a time. Temporary per-chromosome tables are merged in natural
chromosome order. Existing analysis outputs are replaced only after every
selected chromosome has passed its length check. A subset run also replaces
these two outputs with the selected subset. Defaults are independent of cwd.
"""

import argparse
import csv
import gzip
import heapq
import multiprocessing
import os
import queue
import re
import shutil
import sys
import tempfile
import time
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path
from urllib.parse import urlencode


ROOT = Path(__file__).resolve().parent
DEFAULT_GENOME = ROOT / "input" / "Dataset_31"
DEFAULT_RESULTS = DEFAULT_GENOME / "Results" / "Results_Set_2026_11_01_23_42_29"
MOTIF = "XTTCGXXXXAXXCGGGX"
COMPLEMENT = str.maketrans("ACGTX", "TGCAX")
COORDINATES = "1-based inclusive"
SITE_FIELDS = [
    "Chromosome_Name", "Sequence", "Start_Locus", "End_Locus", "Strand",
    "Motif_Oriented_Sequence", "Link_to_Location", "Distance_from_previous",
    "Distance_to_next", "Centromere_Relation", "Assembly", "Coordinate_System",
]
SUMMARY_FIELDS = [
    "Chromosome_Name", "Chromosome_Length", "Total_Sites", "Forward_Sites",
    "Reverse_Sites", "First_Site_Start_Locus", "First_Site_End_Locus",
    "Last_Site_Start_Locus", "Last_Site_End_Locus", "Centromere_Start_Locus",
    "Centromere_End_Locus", "Sites_Inside_Centromere", "Sites_Outside_Centromere",
    "Sites_Overlapping_Centromere_Boundary", "All_Sites_Within_Centromere",
    "First_Site_Link", "Last_Site_Link", "Assembly", "Coordinate_System",
]
_EVENTS = None


def chromosome_name(name):
    """Use the project's chr01 spelling; browser links use UCSC chr1."""
    name = name.strip()
    if name.lower().startswith("chr"):
        name = name[3:]
    if name.isdigit() and 1 <= int(name) <= 22:
        return "chr%02d" % int(name)
    if name.upper() in {"X", "Y", "M", "MT"}:
        return "chr" + ("M" if name.upper() == "MT" else name.upper())
    raise ValueError("Unrecognized chromosome name: %r" % name)


def chromosome_key(name):
    suffix = name[3:]
    return int(suffix) if suffix.isdigit() else {"X": 23, "Y": 24, "M": 25}[suffix]


def location_link(chromosome, start, end):
    suffix = chromosome[3:]
    browser_name = "chr" + (str(int(suffix)) if suffix.isdigit() else suffix)
    return "https://genome.ucsc.edu/cgi-bin/hgTracks?" + urlencode(
        {"db": "hs1", "position": "%s:%d-%d" % (browser_name, start, end)}
    )


def read_sizes(path):
    sizes = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"Chromosome", "Size"}.issubset(reader.fieldnames or []):
            raise ValueError("Sizes CSV must have Chromosome and Size columns: %s" % path)
        for row in reader:
            name, size = chromosome_name(row["Chromosome"]), int(row["Size"])
            if size <= 0 or name in sizes:
                raise ValueError("Invalid or duplicate chromosome size: %s" % name)
            sizes[name] = size
    if not sizes:
        raise ValueError("Empty chromosome sizes file: %s" % path)
    return sizes


def read_centromeres(path, sizes):
    """Convert UCSC [start, end) acen rows to 1-based inclusive intervals."""
    bands, ends = {}, {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or not row[0].strip() or row[0].startswith("#"):
                continue
            if row[0].lower() in {"chrom", "chromosome", "chromosome_name"}:
                continue
            if len(row) != 5:
                raise ValueError("Expected five cytoband columns at line %d" % line_number)
            name = chromosome_name(row[0])
            start, end = int(row[1]), int(row[2])
            if name not in sizes or not 0 <= start < end <= sizes[name]:
                raise ValueError("Cytoband/size mismatch at line %d: %s" % (line_number, name))
            ends[name] = max(ends.get(name, 0), end)
            if row[4].strip().lower() == "acen":
                bands.setdefault(name, []).append((start, end))
    result = {}
    for name, size in sizes.items():
        if name == "chrM":
            continue
        if ends.get(name) != size:
            raise ValueError("Cytobands do not span the supplied chromosome length: %s" % name)
        intervals = sorted(bands.get(name, []))
        if not intervals:
            raise ValueError("No acen annotation for nuclear chromosome %s" % name)
        start, end = intervals[0]
        for next_start, next_end in intervals[1:]:
            if next_start > end:
                raise ValueError("Disjoint acen intervals cannot be treated as one: %s" % name)
            end = max(end, next_end)
        result[name] = (start + 1, end)
    return result


def open_fasta(path):
    return gzip.open(path, "rb") if path.suffix.lower() == ".gz" else path.open("rb")


def fasta_header(handle, path):
    header = handle.readline().lstrip(b"\xef\xbb\xbf").strip()
    if not header.startswith(b">") or not header[1:].strip():
        raise ValueError("Expected a FASTA header at the start of %s" % path)
    return chromosome_name(header[1:].split()[0].decode("ascii"))


def discover_chromosomes(genome_dir, sizes):
    found = {}
    for path in sorted(genome_dir.iterdir()):
        suffixes = [part.lower() for part in path.suffixes]
        if suffixes and suffixes[-1] == ".gz":
            suffixes.pop()
        if not path.is_file() or not suffixes or suffixes[-1] not in {".fa", ".fasta", ".fna", ".txt"}:
            continue
        with open_fasta(path) as handle:
            name = fasta_header(handle, path)
        if name not in sizes:
            raise ValueError("No expected length for FASTA chromosome %s" % name)
        if name in found:
            raise ValueError("Duplicate FASTA chromosome %s: %s and %s" % (name, found[name], path))
        found[name] = path
    missing = set(sizes) - set(found)
    if missing:
        raise ValueError("Missing chromosome FASTA files: " + ", ".join(sorted(missing)))
    return found


def compile_patterns(strand="both"):
    motifs = []
    if strand in {"both", "forward"}:
        motifs.append((MOTIF, "+"))
    if strand in {"both", "reverse"}:
        motifs.append((MOTIF.translate(COMPLEMENT)[::-1], "-"))
    return [(re.compile(("(?=(" + motif.replace("X", "[ACGT]") + "))").encode("ascii")), sign)
            for motif, sign in motifs]


def find_sites(sequence, offset, patterns):
    # Lookahead preserves overlapping matches. Merge the two strand streams by locus.
    def one_strand(pattern, sign):
        for match in pattern.finditer(sequence):
            yield offset + match.start() + 1, sign, match.group(1).decode("ascii")
    streams = [one_strand(pattern, sign) for pattern, sign in patterns]
    yield from heapq.merge(*streams)


def centromere_relation(start, end, centromere):
    if centromere is None:
        return "not_annotated"
    c_start, c_end = centromere
    if c_start <= start and end <= c_end:
        return "inside"
    if start <= c_end and end >= c_start:
        return "overlaps_boundary"
    return "outside"


def worker_init(events):
    global _EVENTS
    _EVENTS = events


def scan_chromosome(name, path, expected_length, centromere, chunk_bytes, strand, output):
    patterns = compile_patterns(strand)
    processed, count, forward, inside, boundary = 0, 0, 0, 0, 0
    tail, previous = b"", None
    first_start = first_end = last_start = last_end = ""
    with open_fasta(path) as source, output.open("w", encoding="utf-8", newline="") as destination:
        if fasta_header(source, path) != name:
            raise ValueError("FASTA header changed: %s" % path)
        writer = csv.writer(destination)
        while True:
            raw = source.read(chunk_bytes)
            if not raw:
                break
            chunk = raw.translate(None, b" \t\r\n").upper()
            if re.search(b"[^ACGTRYSWKMBDHVN]", chunk):
                raise ValueError("Invalid sequence or multiple FASTA records in %s" % path)
            sequence = tail + chunk
            for start, sign, reference in find_sites(sequence, processed - len(tail), patterns):
                end = start + len(MOTIF) - 1
                relation = centromere_relation(start, end, centromere)
                row = [name, reference, start, end, sign,
                       reference if sign == "+" else reference.translate(COMPLEMENT)[::-1],
                       location_link(name, start, end), "", "", relation, "hs1", COORDINATES]
                if previous is not None:
                    gap = start - previous[3] - 1
                    previous[8] = gap
                    row[7] = gap
                    writer.writerow(previous)
                previous = row
                count += 1
                forward += sign == "+"
                inside += relation == "inside"
                boundary += relation == "overlaps_boundary"
                if count == 1:
                    first_start, first_end = start, end
                last_start, last_end = start, end
            processed += len(chunk)
            if processed > expected_length:
                raise ValueError("FASTA exceeds expected length for %s" % name)
            tail = sequence[-(len(MOTIF) - 1):]
            if _EVENTS is not None:
                _EVENTS.put((name, processed, count))
        if previous is not None:
            writer.writerow(previous)
    if processed != expected_length:
        raise ValueError("FASTA length mismatch for %s: read %d; expected %d" % (name, processed, expected_length))
    c_start, c_end = centromere if centromere is not None else ("", "")
    verdict = ("not_annotated" if centromere is None else
               "no_sites" if not count else "yes" if inside == count else "no")
    return dict(zip(SUMMARY_FIELDS, [
        name, processed, count, forward, count - forward, first_start, first_end,
        last_start, last_end, c_start, c_end, inside if centromere else "",
        count - inside if centromere else "", boundary if centromere else "", verdict,
        location_link(name, first_start, first_end) if count else "",
        location_link(name, last_start, last_end) if count else "", "hs1", COORDINATES,
    ]))


class Progress:
    def __init__(self, total, slots, disabled):
        self.bars, self.overall, self.last_log = [], None, 0.0
        if not disabled:
            try:
                from tqdm import tqdm
            except ImportError:
                print("tqdm is unavailable; using periodic progress logs.", flush=True)
            else:
                self.overall = tqdm(total=total, desc="Chromosomes", unit="chr", position=0)
                self.bars = [tqdm(total=0, desc="Waiting", unit="bp", unit_scale=True,
                                  position=index + 1, leave=False) for index in range(slots)]
        self.disabled = disabled

    def start(self, slot, name, size):
        if self.bars:
            bar = self.bars[slot]
            bar.reset(total=size)
            bar.set_description(name)
            bar.set_postfix(sites=0)

    def update(self, slot, name, scanned, count, total_count, completed, total):
        if self.bars:
            bar = self.bars[slot]
            bar.update(scanned - bar.n)
            bar.set_postfix(sites="%s" % format(count, ","), refresh=False)
            self.overall.set_postfix(sites=format(total_count, ","), refresh=False)
        elif not self.disabled and time.monotonic() - self.last_log >= 5:
            print("Progress: %d/%d chromosomes; %s bp scanned on %s; %s sites so far" %
                  (completed, total, format(scanned, ","), name, format(total_count, ",")), flush=True)
            self.last_log = time.monotonic()

    def completed(self, summary, total_count):
        if self.overall:
            self.overall.update(1)
            self.overall.set_postfix(sites=format(total_count, ","))
        else:
            print("Finished %s: %s bp; %s sites" % (summary["Chromosome_Name"],
                  format(summary["Chromosome_Length"], ","), format(summary["Total_Sites"], ",")), flush=True)

    def close(self):
        for bar in reversed(self.bars):
            bar.close()
        if self.overall is not None:
            self.overall.close()


def run_scans(names, files, sizes, centromeres, temporary, workers, chunk_bytes, strand, disabled):
    progress = Progress(len(names), workers, disabled)
    summaries, counts, active_slots = {}, {}, {}
    pending_names = iter(names)
    context = multiprocessing.get_context("spawn")
    try:
        with context.Manager() as manager:
            events = manager.Queue()
            with ProcessPoolExecutor(max_workers=workers, mp_context=context,
                                     initializer=worker_init, initargs=(events,)) as executor:
                futures = {}

                def submit_next(slot):
                    name = next(pending_names, None)
                    if name is None:
                        return
                    active_slots[name] = slot
                    progress.start(slot, name, sizes[name])
                    future = executor.submit(scan_chromosome, name, files[name], sizes[name],
                                             centromeres.get(name), chunk_bytes, strand,
                                             temporary / (name + ".csv"))
                    futures[future] = (name, slot)

                for slot in range(workers):
                    submit_next(slot)
                try:
                    while futures:
                        completed, _ = wait(futures, timeout=0.15, return_when=FIRST_COMPLETED)
                        while True:
                            try:
                                name, scanned, count = events.get_nowait()
                            except queue.Empty:
                                break
                            counts[name] = count
                            if name in active_slots:
                                progress.update(active_slots[name], name, scanned, count,
                                                sum(counts.values()), len(summaries), len(names))
                        for future in completed:
                            name, slot = futures.pop(future)
                            summary = future.result()
                            summaries[name] = summary
                            counts[name] = summary["Total_Sites"]
                            progress.update(slot, name, sizes[name], counts[name],
                                            sum(counts.values()), len(summaries), len(names))
                            progress.completed(summary, sum(counts.values()))
                            del active_slots[name]
                            submit_next(slot)
                except BaseException:
                    for future in futures:
                        future.cancel()
                    raise
    finally:
        progress.close()
    return [summaries[name] for name in names]


def write_outputs(summaries, temporary, output_dir):
    site_path = temporary / "CENP_B_sites.csv"
    summary_path = temporary / "CENP_B_chromosome_summary.csv"
    with site_path.open("w", encoding="utf-8", newline="") as handle:
        csv.writer(handle).writerow(SITE_FIELDS)
        for summary in summaries:
            with (temporary / (summary["Chromosome_Name"] + ".csv")).open(encoding="utf-8", newline="") as fragment:
                shutil.copyfileobj(fragment, handle)
    with summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=SUMMARY_FIELDS)
        writer.writeheader()
        writer.writerows(summaries)
    outputs = [output_dir / site_path.name, output_dir / summary_path.name]
    os.replace(site_path, outputs[0])
    os.replace(summary_path, outputs[1])
    return outputs


def print_plot(summaries, width=70):
    if not summaries:
        return
    maximum = max(row["Chromosome_Length"] for row in summaries)
    print("\nChromosome overview: shared scale 0 .. %.2f Mb" % (maximum / 1_000_000))
    print("C: # = acen; H: |===| = first-to-last site envelope (not continuous occupancy)")
    print("Each nuclear chromosome has separate C and H tracks; - = chromosome.\n")
    for row in summaries:
        length = row["Chromosome_Length"]
        columns = max(1, round(length / maximum * (width - 1)) + 1)

        def track(start, end, character, endpoints=False):
            cells = ["-"] * columns
            if start != "":
                left = min(columns - 1, int((start - 1) / maximum * (width - 1)))
                right = min(columns - 1, int((end - 1) / maximum * (width - 1)))
                cells[left:right + 1] = [character] * (right - left + 1)
                if endpoints:
                    cells[left] = cells[right] = "|"
            return "".join(cells).ljust(width)

        print("%5s C %s acen=%s" % (row["Chromosome_Name"],
              track(row["Centromere_Start_Locus"], row["Centromere_End_Locus"], "#"),
              "%s..%s" % (row["Centromere_Start_Locus"], row["Centromere_End_Locus"])
              if row["Centromere_Start_Locus"] != "" else "unavailable"))
        print("      H %s n=%s; within=%s; span=%s" % (
              track(row["First_Site_Start_Locus"], row["Last_Site_End_Locus"], "=", True),
              format(row["Total_Sites"], ","), row["All_Sites_Within_Centromere"],
              "%s..%s" % (row["First_Site_Start_Locus"], row["Last_Site_End_Locus"])
              if row["Total_Sites"] else "no sites"))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--genome-dir", type=Path, default=DEFAULT_GENOME,
                        help="Directory of single-chromosome FASTA files (.txt also accepted).")
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULTS,
                        help="Result set containing Genome_Details; defaults to the requested 2026-11-01 set.")
    parser.add_argument("--cytobands", type=Path, help="Override Human_Cytobands.tsv (UCSC 0-based half-open).")
    parser.add_argument("--sizes", type=Path, help="Override Human_Genome_Chromosomes_length.csv.")
    parser.add_argument("--workers", type=int, default=min(4, os.cpu_count() or 1),
                        help="Concurrent chromosome processes (default: up to 4). Use 1 to reduce CPU usage.")
    parser.add_argument("--chunk-mb", type=int, default=4, help="Input block size per worker in MiB (default: 4).")
    parser.add_argument("--strand", choices=["both", "forward", "reverse"], default="both")
    parser.add_argument("--chromosomes", nargs="+", help="Optional subset, e.g. chr1 chrX chrM. Default: all.")
    parser.add_argument("--no-progress", action="store_true", help="Print chromosome completion logs only.")
    parser.add_argument("--no-plot", action="store_true", help="Skip terminal chromosome overview.")
    args = parser.parse_args(argv)
    if args.workers < 1 or args.chunk_mb < 1:
        parser.error("--workers and --chunk-mb must be positive integers")
    return args


def main(argv=None):
    args = parse_args(argv)
    output_dir = args.result_dir.resolve() / "Genome_Details"
    sizes_path = args.sizes or output_dir / "Human_Genome_Chromosomes_length.csv"
    cytobands_path = args.cytobands or output_dir / "Human_Cytobands.tsv"
    started = time.monotonic()
    lock_path = output_dir / ".cenp_b_analysis.lock"
    locked = False
    try:
        sizes = read_sizes(sizes_path)
        centromeres = read_centromeres(cytobands_path, sizes)
        files = discover_chromosomes(args.genome_dir.resolve(), sizes)
        names = sorted({chromosome_name(name) for name in args.chromosomes} if args.chromosomes else files,
                       key=chromosome_key)
        if not names or set(names) - set(files):
            raise ValueError("Selected chromosomes are unavailable in the supplied genome")
        output_dir.mkdir(parents=True, exist_ok=True)
        try:
            descriptor = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
        except FileExistsError:
            raise ValueError("Another scan is using this output directory. If a previous process was killed, "
                             "verify it has stopped, then remove %s" % lock_path) from None
        locked = True
        with os.fdopen(descriptor, "w") as handle:
            handle.write("pid=%d\n" % os.getpid())
        workers = min(args.workers, len(names))
        print("CENP-B motif: %s; X = A/C/G/T; strands=%s" % (MOTIF, args.strand), flush=True)
        print("Assembly: hs1; coordinates: %s; %d chromosomes; %d workers" %
              (COORDINATES, len(names), workers), flush=True)
        print("Input: %s\nOutput: %s" % (args.genome_dir.resolve(), output_dir), flush=True)
        print("Scanning full chromosomes, including soft-masked sequence. Ambiguous bases do not match.", flush=True)
        with tempfile.TemporaryDirectory(prefix=".cenp_b_", dir=output_dir) as temporary_dir:
            temporary = Path(temporary_dir)
            summaries = run_scans(names, files, sizes, centromeres, temporary, workers,
                                  args.chunk_mb * 1024 * 1024, args.strand, args.no_progress)
            outputs = write_outputs(summaries, temporary, output_dir)
        if not args.no_plot:
            print_plot(summaries, max(20, min(80, shutil.get_terminal_size((120, 24)).columns - 50)))
        print("\nCompleted: %s sites on %d chromosomes in %.1f seconds." %
              (format(sum(row["Total_Sites"] for row in summaries), ","), len(names), time.monotonic() - started))
        print("Matches are sequence candidates; acen comparison does not establish protein binding.")
        print("Distances count intervening bases: adjacent=0; overlap=negative; chromosome edges=blank.")
        for path in outputs:
            print("Written: %s" % path)
        return 0
    except (OSError, ValueError, csv.Error) as error:
        print("ERROR: %s" % error, file=sys.stderr, flush=True)
        return 1
    except KeyboardInterrupt:
        print("\nInterrupted; scan did not finish.", file=sys.stderr, flush=True)
        return 130
    finally:
        if locked:
            lock_path.unlink(missing_ok=True)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    sys.exit(main())
