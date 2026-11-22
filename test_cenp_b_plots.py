"""Synthetic-data regression tests for CENP-B plots; requires matplotlib."""

import contextlib
import csv
import io
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import cenp_b_analysis as analysis
import cenp_b_plots as plots


def make_fixture(root):
    result = root / "results"
    details = result / "Genome_Details"
    details.mkdir(parents=True)
    with (details / "CENP_B_chromosome_summary.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Chromosome_Name", "Chromosome_Length", "Total_Sites",
                         "Centromere_Start_Locus", "Centromere_End_Locus",
                         "First_Site_Start_Locus", "Last_Site_End_Locus", "All_Sites_Within_Centromere"])
        writer.writerows([("chr01", 250, 6, 101, 150, 1, 250, "no"),
                          ("chrM", 20, 0, "", "", "", "", "not_annotated")])
    with (details / "CENP_B_sites.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["Chromosome_Name", "Start_Locus", "End_Locus"])
        writer.writerows(("chr1", start, start + 16) for start in [1, 100, 101, 200, 201, 234])
    (details / "Human_Cytobands.tsv").write_text(
        "chr1\t0\t100\tp12\tgneg\nchr1\t100\t125\tp11\tacen\n"
        "chr1\t125\t150\tq11\tacen\nchr1\t150\t250\tq12\tgpos50\n",
        encoding="utf-8")
    return result, details


class PlotTests(unittest.TestCase):
    def test_histogram_boundaries_empty_bins_and_short_final_bin(self):
        with tempfile.TemporaryDirectory() as directory:
            _, details = make_fixture(Path(directory))
            chromosomes = plots.read_chromosomes(details / "CENP_B_chromosome_summary.csv")
            histograms = plots.read_histograms(details / "CENP_B_sites.csv", chromosomes, 100)
            histogram = histograms["chr01"]
            self.assertEqual(histogram.edges, [0, 100, 200, 250])
            self.assertEqual(histogram.counts, [2, 2, 2])
            self.assertEqual(histogram.density_per_mb, [20_000, 20_000, 40_000])
            self.assertEqual(histograms["chrM"].counts, [0])
            self.assertEqual(histograms["chrM"].density_per_mb, [0])
            smaller_bins = plots.read_histograms(details / "CENP_B_sites.csv", chromosomes, 25)
            self.assertEqual(smaller_bins["chr01"].counts, [1, 0, 0, 1, 1, 0, 0, 1, 1, 1])
            self.assertEqual(sum(smaller_bins["chr01"].counts), 6)

    def test_cytobands_and_boundary_lines_align_with_histogram(self):
        from matplotlib.colors import to_rgba

        with tempfile.TemporaryDirectory() as directory:
            _, details = make_fixture(Path(directory))
            chromosomes = plots.read_chromosomes(details / "CENP_B_chromosome_summary.csv")
            bands = plots.read_cytobands(details / "Human_Cytobands.tsv", chromosomes)
            histograms = plots.read_histograms(details / "CENP_B_sites.csv", chromosomes, 100)
            figure = plots.create_chromosome_figure(chromosomes["chr01"], bands["chr01"],
                                                    histograms["chr01"], 100)
            try:
                figure.canvas.draw()
                ideogram, histogram = figure.axes
                self.assertEqual(ideogram.get_xlim(), histogram.get_xlim())
                self.assertEqual(histogram.get_xlim(), (0, 250 / 1_000_000))
                self.assertEqual([bar.get_height() for bar in histogram.patches[:3]],
                                 [20_000, 20_000, 40_000])
                self.assertEqual([line.get_xdata()[0] for line in histogram.lines],
                                 [100 / 1_000_000, 150 / 1_000_000])
                self.assertEqual(ideogram.patches[1].get_facecolor(), to_rgba(plots.CYTOBAND_COLORS["acen"]))
                self.assertEqual(ideogram.patches[2].get_facecolor(), to_rgba(plots.CYTOBAND_COLORS["acen"]))
                self.assertEqual(bands["chrM"], [])
            finally:
                figure.clear()

    def test_headless_pngs_and_empty_mitochondrial_chromosome(self):
        from PIL import Image

        with tempfile.TemporaryDirectory() as directory:
            result, _ = make_fixture(Path(directory))
            with contextlib.redirect_stdout(io.StringIO()) as log:
                paths = plots.plot_chromosome_histograms(result, bin_size=100, dpi=100)
            self.assertEqual([path.name for path in paths],
                             ["chr01_CENP_B_distribution.png", "chrM_CENP_B_distribution.png"])
            self.assertTrue(all(path.parent == result / "CENP_B_Plots" for path in paths))
            for path in paths:
                with Image.open(path) as image:
                    self.assertEqual(image.size, (1350, 550))
                    self.assertEqual(image.format, "PNG")
                    image.verify()
            self.assertIn("Written 2 chromosome plots:", log.getvalue())
            self.assertFalse(list(result.glob(".cenp_b_plots_*")))

    def test_inconsistent_tables_do_not_replace_existing_plots(self):
        with tempfile.TemporaryDirectory() as directory:
            result, details = make_fixture(Path(directory))
            output = result / "CENP_B_Plots"
            output.mkdir()
            previous = output / "chr01_CENP_B_distribution.png"
            previous.write_bytes(b"existing plot")
            with (details / "CENP_B_sites.csv").open("a", newline="") as handle:
                csv.writer(handle).writerow(("chr1", 10, 26))
            with contextlib.redirect_stdout(io.StringIO()):
                with self.assertRaisesRegex(ValueError, "Site count differs"):
                    plots.plot_chromosome_histograms(result, bin_size=100)
            self.assertEqual(previous.read_bytes(), b"existing plot")

    def test_centromere_coordinate_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            _, details = make_fixture(Path(directory))
            chromosomes = plots.read_chromosomes(details / "CENP_B_chromosome_summary.csv")
            chromosomes["chr01"] = plots.Chromosome("chr01", 250, 6, (100, 150))
            with self.assertRaisesRegex(ValueError, "acen boundaries differ"):
                plots.read_cytobands(details / "Human_Cytobands.tsv", chromosomes)

    def test_analysis_redraws_changed_bins_without_searching_or_rewriting_csvs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, details = make_fixture(root)
            csv_paths = [details / "CENP_B_sites.csv", details / "CENP_B_chromosome_summary.csv"]
            before = [(path.read_bytes(), path.stat().st_mtime_ns) for path in csv_paths]
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", False), \
                    mock.patch.object(analysis, "read_sizes") as sizes, \
                    mock.patch.object(analysis, "read_centromeres") as centromeres, \
                    mock.patch.object(analysis, "discover_chromosomes") as discover, \
                    mock.patch.object(analysis, "run_scans") as scan, \
                    mock.patch.object(analysis, "write_outputs") as write, \
                    mock.patch.object(analysis, "plot_chromosome_histograms",
                                      wraps=plots.plot_chromosome_histograms) as draw:
                for bin_size in (100, 50):
                    with self.subTest(bin_size=bin_size), contextlib.redirect_stdout(io.StringIO()) as log:
                        self.assertEqual(analysis.main(["--result-dir", str(result), "--genome-dir",
                                                        str(root / "nonexistent_genome"),
                                                        "--plot-bin-size", str(bin_size)]), 0)
                    self.assertIn("Mode: plots only", log.getvalue())
                    self.assertIn("Chromosome overview", log.getvalue())
                    self.assertIn("active search workers: 0", log.getvalue())
                    self.assertEqual(draw.call_args.args[2], bin_size)
                for unused in (sizes, centromeres, discover, scan, write):
                    unused.assert_not_called()
                self.assertEqual(draw.call_count, 2)
            self.assertEqual([(path.read_bytes(), path.stat().st_mtime_ns) for path in csv_paths], before)
            self.assertEqual(len(list((result / "CENP_B_Plots").glob("*.png"))), 2)
            self.assertFalse((details / ".cenp_b_analysis.lock").exists())

    def test_analysis_calls_plotting_and_reports_plot_failure_after_csv_save(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            result, details = make_fixture(root)
            genome = root / "genome"
            genome.mkdir()
            motif = "ATTCGTTGGAACCGGGA"
            (genome / "chr1.fa").write_text(">chr1\n" + motif + "N" * 216 + motif + "\n")
            (genome / "chrM.fa").write_text(">chrM\n" + "N" * 20 + "\n")
            (details / "Human_Genome_Chromosomes_length.csv").write_text(
                "Chromosome,Size\nchr01,250\nchrM,20\n")
            args = ["--genome-dir", str(genome), "--result-dir", str(result), "--workers", "1",
                    "--no-progress", "--no-plot", "--plot-bin-size", "100"]
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", True), \
                    contextlib.redirect_stdout(io.StringIO()) as log:
                self.assertEqual(analysis.main(args), 0)
            self.assertIn("Written 2 chromosome plots:", log.getvalue())
            self.assertEqual(len(list((result / "CENP_B_Plots").glob("*.png"))), 2)
            with mock.patch.object(analysis, "RUN_SITE_ANALYSIS", True), \
                    mock.patch.object(analysis, "plot_chromosome_histograms", side_effect=ValueError("test failure")), \
                    contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()) as error:
                self.assertEqual(analysis.main(args), 1)
            self.assertIn("CSV outputs were saved, but chromosome plotting failed", error.getvalue())
            self.assertTrue((details / "CENP_B_sites.csv").is_file())
            self.assertTrue((details / "CENP_B_chromosome_summary.csv").is_file())
            self.assertFalse((details / ".cenp_b_analysis.lock").exists())


if __name__ == "__main__":
    unittest.main()
