"""Meaningful regression checks: python -m unittest test_cenp_b_analysis -v."""

import contextlib
import csv
import gzip
import io
import random
import tempfile
import unittest
from pathlib import Path
from unittest import mock
from urllib.parse import parse_qs, urlparse

import cenp_b_analysis as analysis


HIT = "ATTCGTTGGAACCGGGA"
REVERSE_HIT = HIT.translate(analysis.COMPLEMENT)[::-1]


def brute_force(sequence, strand="both"):
    sequence = sequence.upper()
    results = []
    for index in range(len(sequence) - 16):
        window = sequence[index:index + 17]
        for sign in ("+", "-"):
            if strand == "forward" and sign == "-" or strand == "reverse" and sign == "+":
                continue
            oriented = window if sign == "+" else window.translate(analysis.COMPLEMENT)[::-1]
            if all(base in "ACGT" and (expected == "X" or base == expected)
                   for base, expected in zip(oriented, analysis.MOTIF)):
                results.append((index + 1, sign, window))
    return results


class MotifTests(unittest.TestCase):
    def test_worker_defaults_reserve_one_logical_cpu(self):
        for detected, expected_cores, expected_workers in [(8, 8, 7), (2, 2, 1), (1, 1, 1), (None, 1, 1)]:
            with self.subTest(detected=detected), mock.patch.object(analysis.os, "cpu_count", return_value=detected):
                args = analysis.parse_args([])
                self.assertEqual(args.cpu_count, expected_cores)
                self.assertEqual(args.workers, expected_workers)

    def test_worker_limit_can_be_overridden(self):
        with mock.patch.object(analysis.os, "cpu_count", return_value=8):
            self.assertEqual(analysis.parse_args(["--workers", "3"]).workers, 3)

    def test_missing_saved_results_do_not_trigger_a_search(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", False), \
                    mock.patch.object(analysis, "discover_chromosomes") as discover, \
                    mock.patch.object(analysis, "run_scans") as scan, \
                    contextlib.redirect_stderr(io.StringIO()) as error:
                self.assertEqual(analysis.main(["--result-dir", directory, "--no-figures"]), 1)
            self.assertIn("Saved CENP-B CSVs are missing", error.getvalue())
            self.assertIn("RUN_SITE_ANALYSIS = True", error.getvalue())
            discover.assert_not_called()
            scan.assert_not_called()
            self.assertFalse((Path(directory) / "Genome_Details").exists())

    def test_regex_agrees_with_independent_basewise_oracle(self):
        rng = random.Random(73)
        # The 2 bp overlap also meets both boxes' fixed-base requirements.
        overlap = HIT[:-2] + "GTTCGTTGGAACCGGGA"
        sequence = "N" + overlap + "N" + REVERSE_HIT + "N" + HIT[:-1] + "N"
        for _ in range(500):
            sequence += "".join(rng.choice("ACGTNRY") for _ in range(30)) + rng.choice([HIT, REVERSE_HIT])
        expected = brute_force(sequence)
        self.assertGreater(len(expected), 500)
        self.assertEqual(list(analysis.find_sites(sequence.encode(), 0, analysis.compile_patterns())), expected)
        starts = [start for start, _, _ in expected]
        self.assertTrue(any(b - a < 17 for a, b in zip(starts, starts[1:])))

    def test_streaming_all_chunk_boundaries_and_strands(self):
        sequence = HIT + "N" + REVERSE_HIT.lower() + HIT.lower() + "R" + HIT[:-1] + "N"
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta = root / "input.fa.gz"
            with gzip.open(fasta, "wt") as handle:
                handle.write(">chr1 description\n" + "\n".join(sequence[i:i + 7] for i in range(0, len(sequence), 7)) + "\n")
            for chunk_bytes in [1, 2, 7, 16, 17, 18, 31, 4096]:
                for strand in ["both", "forward", "reverse"]:
                    with self.subTest(chunk_bytes=chunk_bytes, strand=strand):
                        output = root / "sites.csv"
                        summary = analysis.scan_chromosome("chr01", fasta, len(sequence), (1, 17),
                                                           chunk_bytes, strand, output)
                        with output.open(newline="") as handle:
                            rows = list(csv.DictReader(handle, fieldnames=analysis.SITE_FIELDS))
                        observed = [(int(row["Start_Locus"]), row["Strand"], row["Sequence"]) for row in rows]
                        self.assertEqual(observed, brute_force(sequence, strand))
                        self.assertEqual(summary["Total_Sites"], len(rows))
                        if rows:
                            self.assertEqual(rows[0]["Distance_from_previous"], "")
                            self.assertEqual(rows[-1]["Distance_to_next"], "")
                        for index, row in enumerate(rows):
                            self.assertEqual(int(row["End_Locus"]) - int(row["Start_Locus"]) + 1, 17)
                            self.assertEqual(row["Motif_Oriented_Sequence"], HIT)
                            query = parse_qs(urlparse(row["Link_to_Location"]).query)
                            self.assertEqual(query["db"], ["hs1"])
                            self.assertEqual(query["position"], ["chr1:%s-%s" % (row["Start_Locus"], row["End_Locus"])])
                            if index:
                                expected_gap = int(row["Start_Locus"]) - int(rows[index - 1]["End_Locus"]) - 1
                                self.assertEqual(int(row["Distance_from_previous"]), expected_gap)
                                self.assertEqual(int(rows[index - 1]["Distance_to_next"]), expected_gap)
                        if strand == "both":
                            self.assertEqual(summary["Sites_Inside_Centromere"], 1)
                            self.assertEqual(summary["Sites_Outside_Centromere"], 2)
                            self.assertEqual(summary["All_Sites_Within_Centromere"], "no")

    def test_centromere_boundaries_and_conversion(self):
        self.assertEqual(analysis.centromere_relation(11, 27, (11, 27)), "inside")
        self.assertEqual(analysis.centromere_relation(1, 17, (11, 27)), "overlaps_boundary")
        self.assertEqual(analysis.centromere_relation(27, 43, (11, 27)), "overlaps_boundary")
        self.assertEqual(analysis.centromere_relation(28, 44, (11, 27)), "outside")
        self.assertEqual(analysis.centromere_relation(1, 17, None), "not_annotated")
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cytobands.tsv"
            path.write_text("chr1\t0\t10\tp\tgneg\nchr1\t10\t20\tp11\tacen\n"
                            "chr1\t20\t27\tq11\tacen\nchr1\t27\t50\tq\tgneg\n")
            self.assertEqual(analysis.read_centromeres(path, {"chr01": 50}), {"chr01": (11, 27)})
            with self.assertRaises(ValueError):
                analysis.read_centromeres(path, {"chr01": 51})
            path.write_text("chr1\t0\t10\tp\tacen\nchr1\t20\t50\tq\tacen\n")
            with self.assertRaisesRegex(ValueError, "Disjoint"):
                analysis.read_centromeres(path, {"chr01": 50})

    def test_empty_chromosome_and_input_errors(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, output = root / "input.fa", root / "sites.csv"
            fasta.write_text(">chrM\nACGTNNN\n")
            summary = analysis.scan_chromosome("chrM", fasta, 7, None, 1, "both", output)
            self.assertEqual(summary["Total_Sites"], 0)
            self.assertEqual(summary["First_Site_Start_Locus"], "")
            self.assertEqual(summary["Sites_Inside_Centromere"], "")
            self.assertEqual(summary["All_Sites_Within_Centromere"], "not_annotated")
            self.assertEqual(output.read_text(), "")
            with self.assertRaisesRegex(ValueError, "length mismatch"):
                analysis.scan_chromosome("chrM", fasta, 8, None, 2, "both", output)
            fasta.write_text(">chrM\nACGT\n>chr2\nACGT\n")
            with self.assertRaisesRegex(ValueError, "Invalid sequence"):
                analysis.scan_chromosome("chrM", fasta, 8, None, 10, "both", output)

    def test_summary_centromere_distances(self):
        sequence = "N" * 30 + HIT + "N" * 13 + REVERSE_HIT + "N" * 30
        # First site spans 31..47; last site spans 61..77, independently of strand.
        cases = [
            ("inside", sequence, (11, 97), 20, 20),
            ("on_boundaries", sequence, (31, 77), 0, 0),
            ("left_outside", sequence, (41, 97), -10, 20),
            ("right_outside", sequence, (11, 67), 20, -10),
            ("both_outside", sequence, (41, 67), -10, -10),
            ("no_centromere", sequence, None, "", ""),
            ("no_sites", "N" * 107, (11, 97), "", ""),
            ("single_site", "N" * 30 + HIT + "N" * 30, (11, 67), 20, 20),
        ]
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fasta, output = root / "input.fa", root / "sites.csv"
            for case, dna, centromere, first_distance, last_distance in cases:
                with self.subTest(case=case):
                    fasta.write_text(">chr1\n" + dna + "\n")
                    summary = analysis.scan_chromosome("chr01", fasta, len(dna), centromere,
                                                       7, "both", output)
                    self.assertEqual(summary["First_Site_Distance_From_Centromere_Start_bp"], first_distance)
                    self.assertEqual(summary["Last_Site_Distance_To_Centromere_End_bp"], last_distance)

    def test_multiprocess_outputs_and_failed_rerun_preserves_results(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            genome = root / "genome"
            details = root / "results" / "Genome_Details"
            genome.mkdir()
            details.mkdir(parents=True)
            (genome / "chr1.fa").write_text(">chr1\n" + HIT + REVERSE_HIT + "\n")
            (genome / "chrM.fa").write_text(">chrM\nACGT\n")
            (details / "Human_Genome_Chromosomes_length.csv").write_text("Chromosome,Size\nchr01,34\nchrM,4\n")
            (details / "Human_Cytobands.tsv").write_text("chr1\t0\t34\tp11\tacen\n")
            args = ["--no-figures", "--genome-dir", str(genome), "--result-dir", str(details.parent),
                    "--workers", "2", "--no-progress"]
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", True), \
                    mock.patch.object(analysis.os, "cpu_count", return_value=8), \
                    contextlib.redirect_stdout(io.StringIO()) as log, contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(analysis.main(args[:-3]), 0)
            self.assertIn("Chromosome overview", log.getvalue())
            self.assertIn("Written:", log.getvalue())
            self.assertIn("CPU cores (logical): 8; worker limit: 7; active workers: 2", log.getvalue())
            sites = details / "CENP_B_sites.csv"
            summaries = details / "CENP_B_chromosome_summary.csv"
            with sites.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            with summaries.open(newline="") as handle:
                summary_rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual([row["Chromosome_Name"] for row in summary_rows], ["chr01", "chrM"])
            self.assertEqual(summary_rows[0]["All_Sites_Within_Centromere"], "yes")
            self.assertEqual(summary_rows[0]["First_Site_Distance_From_Centromere_Start_bp"], "0")
            self.assertEqual(summary_rows[0]["Last_Site_Distance_To_Centromere_End_bp"], "0")
            self.assertEqual(summary_rows[1]["First_Site_Distance_From_Centromere_Start_bp"], "")
            self.assertEqual(summary_rows[1]["Last_Site_Distance_To_Centromere_End_bp"], "")
            self.assertEqual(rows[0]["Distance_to_next"], "0")
            before = sites.read_bytes(), summaries.read_bytes()
            (genome / "chr1.fa").write_text(">chr1\nACGT\n")
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", True), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
                self.assertEqual(analysis.main(args), 1)
            self.assertEqual((sites.read_bytes(), summaries.read_bytes()), before)
            self.assertFalse((details / ".cenp_b_analysis.lock").exists())
            self.assertFalse(list(details.glob(".cenp_b_*")))


if __name__ == "__main__":
    unittest.main()
