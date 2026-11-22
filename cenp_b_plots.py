"""Plot CENP-B site density and aligned chromosome cytobands from saved CSVs.

Normally called by cenp_b_analysis.py. To redraw without scanning FASTA:
    python cenp_b_plots.py --result-dir <Results_Set_directory>

Dependencies: matplotlib (the same plotting dependency used by Circos).
PNG files are written to CENP_B_Plots beside Circos_Plots. No GUI/display or
network access is needed. One site is assigned to a bin by its Start_Locus;
strands are combined exactly as present in CENP_B_sites.csv. Empty bins are
retained. Density is sites / actual bin length * 1,000,000, including the
last, possibly shorter bin. Default bin size: 1,000,000 bp.

CSV coordinates are 1-based inclusive. Plot positions measure distance from
the chromosome start: bins and UCSC cytobands use 0-based half-open intervals,
so CSV start positions are shifted by -1 before binning. The red acen region
and both boundary lines use those same continuous plot coordinates.
"""

import argparse
import csv
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path


DEFAULT_BIN_SIZE = 1_000_000
DEFAULT_RESULTS = (Path(__file__).resolve().parent / "input" / "Dataset_31" /
                   "Results" / "Results_Set_2026_11_01_23_42_29")
CYTOBAND_COLORS = {
    "gneg": "#ffffff", "gpos25": "#c8c8c8", "gpos50": "#969696",
    "gpos75": "#646464", "gpos100": "#202020", "gvar": "#e3e3e3",
    "stalk": "#9ecfd4", "acen": "#d73027",
}


@dataclass(frozen=True)
class Chromosome:
    name: str
    length: int
    total_sites: int
    centromere: object  # None or a 1-based inclusive (start, end) pair.


@dataclass(frozen=True)
class Cytoband:
    start: int  # 0-based, inclusive.
    end: int    # 0-based, exclusive.
    name: str
    stain: str


@dataclass
class Histogram:
    edges: list
    counts: list
    density_per_mb: list


def _chromosome_name(name):
    suffix = name.strip()
    if suffix.lower().startswith("chr"):
        suffix = suffix[3:]
    if suffix.isdigit() and 1 <= int(suffix) <= 22:
        return "chr%02d" % int(suffix)
    if suffix.upper() in {"X", "Y", "M", "MT"}:
        return "chr" + ("M" if suffix.upper() == "MT" else suffix.upper())
    raise ValueError("Unrecognized chromosome name: %r" % name)


def _chromosome_key(name):
    suffix = name[3:]
    return int(suffix) if suffix.isdigit() else {"X": 23, "Y": 24, "M": 25}[suffix]


def check_plotting_dependencies():
    try:
        from matplotlib.backends.backend_agg import FigureCanvasAgg  # noqa: F401
    except ImportError as error:
        raise RuntimeError("Plotting requires matplotlib in the Python environment running this script. "
                           "Install matplotlib, or run the analysis with --no-figures.") from error


def read_chromosomes(path):
    result = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"Chromosome_Name", "Chromosome_Length", "Total_Sites",
                    "Centromere_Start_Locus", "Centromere_End_Locus"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Missing required columns in CENP-B chromosome summary")
        for row in reader:
            name = _chromosome_name(row["Chromosome_Name"])
            length, count = int(row["Chromosome_Length"]), int(row["Total_Sites"])
            if name in result or length <= 0 or count < 0:
                raise ValueError("Invalid or duplicate chromosome summary: %s" % name)
            if row.get("Assembly", "hs1") != "hs1":
                raise ValueError("Expected hs1 assembly in chromosome summary: %s" % name)
            if row.get("Coordinate_System", "1-based inclusive") != "1-based inclusive":
                raise ValueError("Expected 1-based inclusive summary coordinates")
            start, end = row["Centromere_Start_Locus"].strip(), row["Centromere_End_Locus"].strip()
            if bool(start) != bool(end):
                raise ValueError("Incomplete centromere bounds for %s" % name)
            centromere = (int(start), int(end)) if start else None
            if centromere is not None and not 1 <= centromere[0] <= centromere[1] <= length:
                raise ValueError("Centromere outside chromosome: %s" % name)
            result[name] = Chromosome(name, length, count, centromere)
    if not result:
        raise ValueError("Empty CENP-B chromosome summary")
    return result


def read_cytobands(path, chromosomes):
    result = {name: [] for name in chromosomes}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        for line_number, row in enumerate(csv.reader(handle, delimiter="\t"), 1):
            if not row or not row[0].strip() or row[0].startswith("#"):
                continue
            if row[0].lower() in {"chrom", "chromosome", "chromosome_name"}:
                continue
            if len(row) != 5:
                raise ValueError("Expected five cytoband columns at line %d" % line_number)
            name = _chromosome_name(row[0])
            if name not in chromosomes:
                continue  # A subset analysis still uses the complete cytoband file.
            start, end, stain = int(row[1]), int(row[2]), row[4].strip().lower()
            if not 0 <= start < end <= chromosomes[name].length or stain not in CYTOBAND_COLORS:
                raise ValueError("Invalid cytoband coordinates/stain at line %d" % line_number)
            result[name].append(Cytoband(start, end, row[3], stain))
    for name, chromosome in chromosomes.items():
        bands = sorted(result[name], key=lambda band: band.start)
        result[name] = bands
        if not bands:
            if name == "chrM" and chromosome.centromere is None:
                continue
            raise ValueError("Missing cytobands for %s" % name)
        previous_end = 0
        for band in bands:
            if band.start != previous_end:
                raise ValueError("Gap or overlap between cytobands on %s" % name)
            previous_end = band.end
        if previous_end != chromosome.length:
            raise ValueError("Cytobands and chromosome length differ for %s" % name)
        acen = [band for band in bands if band.stain == "acen"]
        for left, right in zip(acen, acen[1:]):
            if left.end != right.start:
                raise ValueError("Disjoint acen annotations on %s" % name)
        bounds = (acen[0].start + 1, acen[-1].end) if acen else None
        if bounds != chromosome.centromere:
            raise ValueError("Summary and cytoband acen boundaries differ for %s" % name)
    return result


def read_histograms(path, chromosomes, bin_size=DEFAULT_BIN_SIZE):
    """Stream the occurrence CSV once; retain only per-bin counts in memory."""
    if not isinstance(bin_size, int) or bin_size < 1:
        raise ValueError("Histogram bin size must be a positive integer in bp")
    histograms = {}
    for name, chromosome in chromosomes.items():
        edges = list(range(0, chromosome.length, bin_size)) + [chromosome.length]
        histograms[name] = Histogram(edges, [0] * (len(edges) - 1), [])
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not {"Chromosome_Name", "Start_Locus", "End_Locus"}.issubset(reader.fieldnames or []):
            raise ValueError("Missing required columns in CENP-B sites CSV")
        for row in reader:
            name = _chromosome_name(row["Chromosome_Name"])
            if name not in chromosomes:
                raise ValueError("Sites CSV contains a chromosome absent from summary: %s" % name)
            start, end = int(row["Start_Locus"]), int(row["End_Locus"])
            if not 1 <= start <= end <= chromosomes[name].length or end - start + 1 != 17:
                raise ValueError("Invalid 17 bp motif coordinates on %s" % name)
            if row.get("Coordinate_System", "1-based inclusive") != "1-based inclusive":
                raise ValueError("Expected 1-based inclusive site coordinates")
            histograms[name].counts[(start - 1) // bin_size] += 1
    for name, histogram in histograms.items():
        if sum(histogram.counts) != chromosomes[name].total_sites:
            raise ValueError("Site count differs between occurrence CSV and summary for %s" % name)
        histogram.density_per_mb = [count * 1_000_000 / (right - left)
                                    for count, left, right in zip(histogram.counts, histogram.edges,
                                                                  histogram.edges[1:])]
    return histograms


def create_chromosome_figure(chromosome, bands, histogram, bin_size=DEFAULT_BIN_SIZE, dpi=200):
    """Return a headless Matplotlib figure with two aligned genomic tracks."""
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.ticker import MaxNLocator, ScalarFormatter

    figure = Figure(figsize=(13.5, 5.5), dpi=dpi, constrained_layout=True, facecolor="white")
    FigureCanvasAgg(figure)
    grid = figure.add_gridspec(2, 1, height_ratios=(0.8, 4), hspace=0.06)
    ideogram = figure.add_subplot(grid[0])
    distribution = figure.add_subplot(grid[1], sharex=ideogram)
    figure.suptitle("CENP-B motif distribution | %s | T2T / hs1" % chromosome.name,
                   fontsize=16, fontweight="bold", color="#22354a")
    ideogram.set_title("Chromosome cytobands", loc="left", fontsize=9, color="#526174", pad=10)
    if bands:
        for band in bands:
            patch = Rectangle((band.start / 1_000_000, 0.25), (band.end - band.start) / 1_000_000,
                              0.5, facecolor=CYTOBAND_COLORS[band.stain], edgecolor="#707070",
                              linewidth=0.25, hatch="///" if band.stain == "gvar" else None)
            ideogram.add_patch(patch)
    else:
        ideogram.add_patch(Rectangle((0, 0.25), chromosome.length / 1_000_000, 0.5,
                                    facecolor="#eeeeee", edgecolor="#707070", linewidth=0.5))
        ideogram.text(0.5, 0.5, "Cytobands unavailable", transform=ideogram.transAxes,
                      ha="center", va="center", fontsize=9, color="#526174")
    ideogram.add_patch(Rectangle((0, 0.25), chromosome.length / 1_000_000, 0.5,
                                facecolor="none", edgecolor="#27384a", linewidth=0.8))
    legend = [Patch(facecolor=CYTOBAND_COLORS[key], edgecolor="#707070", label=label)
              for key, label in [("gneg", "gneg"), ("gpos75", "gpos (gray levels)"), ("acen", "acen")]]
    if any(band.stain == "stalk" for band in bands):
        legend.append(Patch(facecolor=CYTOBAND_COLORS["stalk"], label="stalk"))
    if any(band.stain == "gvar" for band in bands):
        legend.append(Patch(facecolor=CYTOBAND_COLORS["gvar"], hatch="///", label="gvar"))
    if bands:
        ideogram.legend(handles=legend, loc="lower right", bbox_to_anchor=(1, 1.02),
                        ncol=len(legend), frameon=False, fontsize=8)
    ideogram.set_ylim(0, 1)
    ideogram.set_yticks([])
    ideogram.tick_params(axis="x", which="both", bottom=False, labelbottom=False)
    for spine in ideogram.spines.values():
        spine.set_visible(False)

    left_edges = [edge / 1_000_000 for edge in histogram.edges[:-1]]
    widths = [(right - left) / 1_000_000 for left, right in zip(histogram.edges, histogram.edges[1:])]
    distribution.bar(left_edges, histogram.density_per_mb, width=widths, align="edge",
                     color="#356d9b", edgecolor="none", linewidth=0, zorder=2)
    distribution.set_xlim(0, chromosome.length / 1_000_000)
    distribution.set_ylim(bottom=0)
    distribution.set_xlabel("Position from chromosome start (Mb)", fontsize=11)
    distribution.set_ylabel("CENP-B motif sites / Mb", fontsize=11)
    distribution.set_title("%s sites | bin size: %s bp | chromosome length: %.3f Mb" %
                           (format(chromosome.total_sites, ","), format(bin_size, ","),
                            chromosome.length / 1_000_000), loc="left", fontsize=9, pad=10,
                           color="#526174")
    distribution.xaxis.set_major_locator(MaxNLocator(nbins=10, min_n_ticks=3))
    distribution.yaxis.set_major_locator(MaxNLocator(nbins=6))
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    distribution.yaxis.set_major_formatter(formatter)
    distribution.grid(axis="y", color="#e0e6ec", linewidth=0.7, zorder=0)
    distribution.spines["top"].set_visible(False)
    distribution.spines["right"].set_visible(False)
    if chromosome.centromere is not None:
        start, end = chromosome.centromere
        left, right = (start - 1) / 1_000_000, end / 1_000_000
        distribution.axvspan(left, right, color=CYTOBAND_COLORS["acen"], alpha=0.07, zorder=1)
        distribution.axvline(left, color="#d73027", linestyle="--", linewidth=1.3,
                             label="acen start: %.3f Mb" % left, zorder=3)
        distribution.axvline(right, color="#8e181c", linestyle=":", linewidth=1.5,
                             label="acen end: %.3f Mb" % right, zorder=3)
        distribution.legend(loc="upper right", fontsize=9, framealpha=0.95)
    else:
        distribution.text(0.99, 0.97, "No acen annotation", transform=distribution.transAxes,
                           ha="right", va="top", fontsize=9, color="#526174")
    if chromosome.total_sites == 0:
        distribution.set_ylim(0, 1)
        distribution.text(0.5, 0.45, "No motif sites detected", transform=distribution.transAxes,
                           ha="center", fontsize=12, color="#526174")
    return figure


def plot_chromosome_histograms(result_dir, cytobands_path=None, bin_size=DEFAULT_BIN_SIZE, dpi=200):
    """Render one PNG per summary chromosome, without reading genome FASTA."""
    check_plotting_dependencies()
    if not isinstance(bin_size, int) or bin_size < 1 or dpi < 1:
        raise ValueError("Bin size and DPI must be positive integers")
    result_dir = Path(result_dir).resolve()
    details = result_dir / "Genome_Details"
    chromosomes = read_chromosomes(details / "CENP_B_chromosome_summary.csv")
    bands = read_cytobands(cytobands_path or details / "Human_Cytobands.tsv", chromosomes)
    print("Building CENP-B histograms from saved CSVs (bin size: %s bp)..." % format(bin_size, ","), flush=True)
    histograms = read_histograms(details / "CENP_B_sites.csv", chromosomes, bin_size)
    output_dir = result_dir / "CENP_B_Plots"
    paths = []
    # Stage all figures before replacing existing PNGs; close every figure.
    with tempfile.TemporaryDirectory(prefix=".cenp_b_plots_", dir=result_dir) as temporary_dir:
        for index, name in enumerate(sorted(chromosomes, key=_chromosome_key), 1):
            figure = create_chromosome_figure(chromosomes[name], bands[name], histograms[name], bin_size, dpi)
            filename = "%s_CENP_B_distribution.png" % name
            try:
                figure.savefig(Path(temporary_dir) / filename, format="png", dpi=dpi)
            finally:
                figure.clear()
            paths.append(output_dir / filename)
            print("Plot %d/%d: %s (%s sites)" %
                  (index, len(chromosomes), name, format(chromosomes[name].total_sites, ",")), flush=True)
        output_dir.mkdir(parents=True, exist_ok=True)
        for path in paths:
            os.replace(Path(temporary_dir) / path.name, path)
    print("Written %d chromosome plots: %s" % (len(paths), output_dir), flush=True)
    return paths


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--cytobands", type=Path, help="Override Human_Cytobands.tsv.")
    parser.add_argument("--bin-size", type=int, default=DEFAULT_BIN_SIZE, help="Histogram window size in bp (default: 1000000).")
    parser.add_argument("--dpi", type=int, default=200)
    args = parser.parse_args(argv)
    if args.bin_size < 1 or args.dpi < 1:
        parser.error("--bin-size and --dpi must be positive integers")
    try:
        plot_chromosome_histograms(args.result_dir, args.cytobands, args.bin_size, args.dpi)
    except (OSError, ValueError, csv.Error, RuntimeError) as error:
        parser.exit(1, "ERROR: %s\n" % error)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
