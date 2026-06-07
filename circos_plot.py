import os
import io
import pickle
import re
from datetime import datetime
from urllib.error import URLError
from pathlib import Path
import pandas as pd
import numpy as np
from pandas import Series
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import utilities
from pycirclize import Circos
from pycirclize.utils import ColorCycler, load_eukaryote_example_dataset
from PIL import Image

def chr_sort_key(chr_name: str):
    # Natural sort key for human chromosomes: chr1..chr22, chrX, chrY, chrM (optional)
    m = re.fullmatch(r"chr(\d+)", chr_name)
    if m:
        return (0, int(m.group(1)))
    if chr_name == "chrX":
        return (1, 23)
    if chr_name == "chrY":
        return (1, 24)
    if chr_name == "chrM":
        return (2, 25)
    return (3, chr_name)

def normalize_chr(name: str) -> str:
    # Convert chr01 -> chr1, chr002 -> chr2; keep chrX/chrY/chrM unchanged
    name = str(name)
    m = re.fullmatch(r"(chr)0*([0-9]+)$", name)
    if m:
        return f"Chr{int(m.group(2))}"
    return name

def safe_filename(text: str) -> str:
    text = str(text)
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = text.replace("-", "")
    text = text.strip("_")
    return text

def plot(genome_folder_address, primary_patterns_file_address, final_results_data_path, base_str):
    target_folder_results = "Results_Set_2026_02_16_14_56_49"
    folder_of_current_results = genome_folder_address + os.sep + "Results" + os.sep + target_folder_results
    folder_of_colonies_location = folder_of_current_results + os.sep + "Colonies_members"
    folder_of_genome_details = folder_of_current_results + os.sep + "Genome_Details"
    folder_of_circos_plots = folder_of_current_results + os.sep + "Circos_Plots"

    Path(folder_of_circos_plots).mkdir(parents=True, exist_ok=True)

    # if not final_results_data_path:
    #     folder_of_current_data_structures = folder_of_current_results + os.sep + "Data_Structures"
    #     final_results_data_path = os.path.join(folder_of_current_data_structures, "total_colonies_data_results.pkl")
    #
    # total_colonies_results = None
    # # read data structures
    # with open(final_results_data_path, "rb") as f:  # rb = read binary
    #     total_colonies_results = pickle.load(f)

    # =================================================================================================================================================
    # -----------------------
    # Load chromosome sizes (hg38)
    # -----------------------
    sectors = {}
    # try:
    #     chromosome_sizes_url = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes"
    #     ucsc_chromosome_size_df = pd.read_csv(chromosome_sizes_url, sep="\t", header=None, names=["Chromosome", "Size"])
    #
    #     ucsc_chromosome_name = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    #     ucsc_chromosome_size_df = ucsc_chromosome_size_df[ucsc_chromosome_size_df["Chromosome"].isin(ucsc_chromosome_name)].copy()
    #     ucsc_chromosome_size_df["Chromosome"] = ucsc_chromosome_size_df["Chromosome"].map(utilities.standard_name)
    #
    #     # dict: {"chr1": 248956422, ...}
    #     sectors = dict(zip(ucsc_chromosome_size_df["Chromosome"], ucsc_chromosome_size_df["Size"]))
    # except URLError:
    #     # primary_chromosome_dict = utilities.read_chromosomes(genome_folder_address)
    #     # chromosome_length_dict = {}
    #     # for primary_chromosome_name, dna_sequence in primary_chromosome_dict.items():
    #     #     chromosome_name = utilities.standard_name(primary_chromosome_name.split("_")[0])
    #     #     chromosome_length_dict[chromosome_name] = len(dna_sequence)
    #     # df = pd.DataFrame(list(chromosome_length_dict.items()), columns=["Chromosome", "Size"])
    #     # df.sort_values(by="Chromosome", ascending=True, inplace=True)
    #     # df.to_csv(os.path.join(folder_of_genome_details, "Human_Genome_Chromosomes_length.csv"), index=False, encoding="utf-8-sig")
    #     sectors = pd.read_csv(os.path.join(folder_of_genome_details, "Human_Genome_Chromosomes_length.csv"), encoding="utf-8-sig")
    #     sectors = sectors.drop(sectors[sectors["Chromosome"] == "chrM"].index)
    #     sectors = sectors.sort_values(by="Chromosome", ascending=True, inplace=True)
    #     sectors = dict(zip(sectors_df["Chromosome"], sectors_df["Size"]))

    sectors_df = pd.read_csv(os.path.join(folder_of_genome_details, "Human_Genome_Chromosomes_length.csv"), encoding="utf-8-sig")
    sectors_df = sectors_df.drop(sectors_df[sectors_df["Chromosome"] == "chrM"].index)
    sectors_df = sectors_df.sort_values(by="Chromosome", ascending=True)
    sectors = dict(zip(sectors_df["Chromosome"], sectors_df["Size"]))

    # -----------------------
    # Load STR occurrences (STR_Map file)
    # -----------------------
    # base_str = "(GT)3-(TG)3"
    # base_str = "(GT)3"
    # base_str = "(TG)3"
    # base_str = "(AC)3-(CA)3"
    # base_str = "(AC)3"
    # base_str = "(CA)3"
    # base_str = "(AG)3-(GA)3"
    # base_str = "(AT)3-(TA)3"
    # base_str = "(CG)3-(GC)3"
    # base_str = "(CT)3-(TC)3"
    chromosome_str_map_dict = pd.read_csv(os.path.join(folder_of_colonies_location, F"Human_Integrated_D3Us_{base_str}_Map_Result.csv"))

    # -----------------------
    # Bin events -> density per bin
    # -----------------------
    circos_bin = 1_000_000  # 1Mb
    median = ((chromosome_str_map_dict["Start_Locus"].to_numpy() + chromosome_str_map_dict["End_Locus"].to_numpy()) // 2).astype(np.int64)
    chromosome_str_map_dict["Bin_Median"] = median
    chromosome_str_map_dict["Bin_Index"] = (chromosome_str_map_dict["Bin_Median"] // circos_bin).astype(np.int64)

    # count per (chrom, bin)
    binned_df = chromosome_str_map_dict.groupby(["Chromosome_Name", "Bin_Index"]).size().reset_index(name="Bin_STR_Count")

    # x-position for drawing bar chart by using of bins center
    binned_df["Bin_Center"] = binned_df["Bin_Index"] * circos_bin + circos_bin // 2

    # correcting last bin center in each track1
    for chromosome_name, chromosome_length in sectors.items():
        last_bin_index = binned_df.loc[binned_df["Chromosome_Name"] == chromosome_name, "Bin_Index"].max()
        if (last_bin_index + 1) * circos_bin <= chromosome_length :
            continue
        last_bin_start = last_bin_index * circos_bin
        length_of_last_bin = chromosome_length - last_bin_start
        if length_of_last_bin == circos_bin :
            raise IndexError("index out of range")
        mask = ((binned_df["Chromosome_Name"] == chromosome_name) & (binned_df["Bin_Index"] == last_bin_index))
        binned_df.loc[mask, "Bin_Center"] = last_bin_start + length_of_last_bin // 2

    binned_df["STR_Density_Per_MB"] = binned_df["Bin_STR_Count"] * (1_000_000 / circos_bin)

    binned_df["STR_Log1e_Density_Per_MB"] = np.log1p(binned_df["STR_Density_Per_MB"])
    binned_df["STR_Log2_Density_Per_MB"] = np.log2(binned_df["STR_Density_Per_MB"] + 1)

    median_baseline = np.median(binned_df["STR_Density_Per_MB"])
    mean_baseline = np.mean(binned_df["STR_Density_Per_MB"])
    binned_df["STR_Log2_Enrich_Median_Density_Per_MB"] = np.log2((binned_df["STR_Density_Per_MB"] + 1)/(median_baseline + 1))
    binned_df["STR_Log2_Enrich_Mean_Density_Per_MB"] = np.log2((binned_df["STR_Density_Per_MB"] + 1)/(mean_baseline + 1))

    # -----------------------
    # Make colony map
    # -----------------------

    grouped_map = chromosome_str_map_dict.groupby(["Chromosome_Name", "Colony_Index", "Colony_Member_Count", "Colony_Type"], as_index=False)
    colonies_map = (grouped_map.agg(member_count_check=("Core", "size"), Colony_Start=("Start_Locus", "min"), Colony_End=("End_Locus", "max")))
    colonies_map["Colony_length"] = (colonies_map["Colony_End"].to_numpy() - colonies_map["Colony_Start"].to_numpy())
    colonies_map["Colony_Bin_Median"] = ((colonies_map["Colony_Start"].to_numpy() + colonies_map["Colony_End"].to_numpy()) // 2).astype(np.int64)

    filtered_colony_map = colonies_map[colonies_map["Colony_Member_Count"] > 2].copy()
    filtered_colony_map["Colony_Bin_Index"] = (filtered_colony_map["Colony_Bin_Median"] // circos_bin).astype(np.int64)

    # count per (chrom, bin)
    colonies_binned_df = filtered_colony_map.groupby(["Chromosome_Name", "Colony_Bin_Index"]).size().reset_index(name="Bin_Colony_Count")

    # x-position for drawing bar chart by using of bins center
    colonies_binned_df["Colony_Bin_Center"] = colonies_binned_df["Colony_Bin_Index"] * circos_bin + circos_bin // 2

    # correcting last bin center in each track1
    for chromosome_name, chromosome_length in sectors.items():
        last_bin_index = colonies_binned_df.loc[colonies_binned_df["Chromosome_Name"] == chromosome_name, "Colony_Bin_Index"].max()
        if (last_bin_index + 1) * circos_bin <= chromosome_length:
            continue
        last_bin_start = last_bin_index * circos_bin
        length_of_last_bin = chromosome_length - last_bin_start
        if length_of_last_bin == circos_bin:
            raise IndexError("index out of range")
        mask = ((colonies_binned_df["Chromosome_Name"] == chromosome_name) & (colonies_binned_df["Colony_Bin_Index"] == last_bin_index))
        colonies_binned_df.loc[mask, "Colony_Bin_Center"] = last_bin_start + length_of_last_bin // 2

    colonies_binned_df["Colony_Density_Per_MB"] = colonies_binned_df["Bin_Colony_Count"] * (1_000_000 / circos_bin)

    colonies_binned_df["Colony_Log1e_Density_Per_MB"] = np.log1p(colonies_binned_df["Colony_Density_Per_MB"])
    colonies_binned_df["Colony_Log2_Density_Per_MB"] = np.log2(colonies_binned_df["Colony_Density_Per_MB"] + 1)

    median_baseline = np.median(colonies_binned_df["Colony_Density_Per_MB"])
    mean_baseline = np.mean(colonies_binned_df["Colony_Density_Per_MB"])
    colonies_binned_df["Colony_Log2_Enrich_Median_Density_Per_MB"] = np.log2((colonies_binned_df["Colony_Density_Per_MB"] + 1) / (median_baseline + 1))
    colonies_binned_df["Colony_Log2_Enrich_Mean_Density_Per_MB"] = np.log2((colonies_binned_df["Colony_Density_Per_MB"] + 1) / (mean_baseline + 1))

    # -----------------------
    # Build Circos with pyCirclize
    # -----------------------

    if not Path(os.path.join(folder_of_genome_details, "Human_Cytobands.tsv")).exists():
        # Load hg38 example cytoband file from pyCirclize dataset
        chr_bed_file, cytoband_file, chr_links = load_eukaryote_example_dataset("hg38")
        cytobands_df = pd.read_csv(cytoband_file, sep="\t", header=None, names=["Chromosome", "ChromosomeStart", "ChromosomeEnd", "name", "gieStain"])
        cytobands_df["Chromosome"] = cytobands_df["Chromosome"].map(utilities.standard_name)
        cytobands_df = cytobands_df[cytobands_df["Chromosome"] != "chrM"]
        cytobands_df.to_csv(os.path.join(folder_of_genome_details, "Human_Cytobands.tsv"), sep="\t", header=False, index=False)

    circos = Circos(sectors, space=2)
    # circos.text("Homo sapiens (hg38)", deg=360, r=150, size=25)
    # circos.text(F"Homo sapiens (hg38)\n\n{base_str}\nD3U Density", deg=360, r=0, size=15)
    circos.text(F"{base_str}\nD3U Density", deg=360, r=0, size=15)

    # Add cytoband track1 (this draws chromosome band pattern)
    circos.add_cytoband_tracks((97, 100), os.path.join(folder_of_genome_details, "Human_Cytobands.tsv"), track_name="cytoband")

    for sector in circos.sectors:
        sector.axis(lw=0.8, ec="none", fc="none")
        sector.text(normalize_chr(sector.name), r=120, size=10)

        sector.get_track("cytoband").xticks_by_interval(
            40000000,
            label_size=8,
            label_orientation="vertical",
            label_formatter=lambda v: f"{v / 1000000:.0f} Mb",
        )

        # track1 for histogram
        track1 = sector.add_track((74, 94), r_pad_ratio=0.05)
        track1.axis(lw=0.6)
        # global_sector_y_max = binned_df["STR_Log1e_Density_Per_MB"].max()
        # df_s = binned_df[binned_df["Chromosome_Name"] == sector.name]
        global_sector_y_max = colonies_binned_df["Colony_Density_Per_MB"].max()
        df_s = colonies_binned_df[colonies_binned_df["Chromosome_Name"] == sector.name]
        if len(df_s) == 0:
            continue
        # x = df_s["Bin_Center"].to_numpy()
        # y = df_s["STR_Log1e_Density_Per_MB"].to_numpy()
        x = df_s["Colony_Bin_Center"].to_numpy()
        y = df_s["Colony_Density_Per_MB"].to_numpy()

        # track1.bar(x, y, vmax=global_sector_y_max + 2, width=circos_bin, color="blue")
        track1.bar(x, y, vmax=global_sector_y_max + 2, width=circos_bin, color="green")

        # track2 for histogram
        track2 = sector.add_track((51, 71), r_pad_ratio=0.05)
        track2.axis(lw=0.6)
        global_sector_y_max = binned_df["STR_Density_Per_MB"].max()
        df_s = binned_df[binned_df["Chromosome_Name"] == sector.name]
        if len(df_s) == 0:
            continue
        x = df_s["Bin_Center"].to_numpy()
        y = df_s["STR_Density_Per_MB"].to_numpy()

        track2.bar(x, y, vmax=global_sector_y_max + 2, width=circos_bin, color="red")

        # # track3 for histogram
        # track3 = sector.add_track((52, 64), r_pad_ratio=0.05)
        # track3.axis(lw=0.6)
        # global_sector_y_max = colonies_binned_df["Colony_Density_Per_MB"].max()
        # df_s = colonies_binned_df[colonies_binned_df["Chromosome_Name"] == sector.name]
        # if len(df_s) == 0:
        #     continue
        # x = df_s["Colony_Bin_Center"].to_numpy()
        # y = df_s["Colony_Density_Per_MB"].to_numpy()
        #
        # track3.bar(x, y, vmax=global_sector_y_max + 2, width=circos_bin, color="green")

        # # track4 for histogram
        # track4 = sector.add_track((37, 49), r_pad_ratio=0.05)
        # track4.axis(lw=0.6)
        # global_sector_y_max = binned_df["STR_Log2_Enrich_Median_Density_Per_MB"].max()
        # df_s = binned_df[binned_df["Chromosome_Name"] == sector.name]
        # if len(df_s) == 0:
        #     continue
        # x = df_s["Bin_Center"].to_numpy()
        # y = df_s["STR_Log2_Enrich_Median_Density_Per_MB"].to_numpy()
        #
        # track4.bar(x, y, width=circos_bin, color="green")


    fig = circos.plotfig()
    now = datetime.now()
    format_now = now.strftime("%Y_%m_%d_%H_%M_%S")

    safe_base_str = safe_filename(base_str)
    png_path = os.path.join(folder_of_circos_plots, f"genome_density_{format_now}_{safe_base_str}.png")
    jpg_path = os.path.join(folder_of_circos_plots, f"genome_density_{format_now}_{safe_base_str}.jpg")

    fig.savefig(png_path, dpi=600, bbox_inches="tight", pad_inches=0.05, facecolor="white")
    fig.savefig(jpg_path, dpi=600, bbox_inches="tight", pad_inches=0.05, facecolor="white", pil_kwargs={"quality": 100, "subsampling": 0})

    panel_image = fig_to_pil_image(fig, dpi=600, bbox_inches="tight", pad_inches=0.05)

    plt.show()
    plt.close(fig)
    return {"base_str": base_str, "png_path": png_path, "jpg_path": jpg_path, "panel_image": panel_image, "output_folder": folder_of_circos_plots}

    # =================================================================================================================================================


    # # Create chromosome color dict
    # ColorCycler.set_cmap("hsv")
    # chr_names = [s.name for s in circos.sectors]
    # colors = ColorCycler.get_color_list(len(chr_names))
    # chr_name2color = {name: color for name, color in zip(chr_names, colors)}
    #
    # # Plot chromosome name
    # for sector in circos.sectors:
    #     sector.text(sector.name, r=120, size=10, color=chr_name2color[sector.name])
    #     # sector.text(sector.name, size=10)
    #     sector.get_track("cytoband").xticks_by_interval(
    #         40000000,
    #         label_size=8,
    #         label_orientation="vertical",
    #         label_formatter=lambda v: f"{v / 1000000:.0f} Mb",
    #     )
    #
    #
    # fig = circos.plotfig()
    # # =============================
    #
    # primary_chromosome_dict = utilities.read_chromosomes(genome_folder_address)
    # patterns_list = utilities.read_patterns(primary_patterns_file_address)
    #
    # core_length = 2
    # base_str = "(GT)3-(TG)3"
    # chromosome_length_dict = {}
    # chromosome_ste_map_dict = {}
    # for primary_chromosome_name, dna_sequence in primary_chromosome_dict.items():
    #     chromosome_name = utilities.standard_name(primary_chromosome_name.split("_")[0])
    #     chromosome_length_dict[chromosome_name] = len(dna_sequence)
    #
    #     chromosome_ste_map_dict[chromosome_name] = (
    #         pd.read_csv(os.path.join(folder_of_colonies_location, F"Human_{chromosome_name}_D3Us_{base_str}_Map_Result.csv"), header = 0, index_col = 0))


def fig_to_pil_image(fig, dpi=600, bbox_inches="tight", pad_inches=0.05):
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=dpi, bbox_inches=bbox_inches, pad_inches=pad_inches, facecolor="white")
    buffer.seek(0)
    image = Image.open(buffer).convert("RGB")
    return image


def plot_batch(genome_folder_address, primary_patterns_file_address, final_results_data_path, base_str_list):
    if len(base_str_list) != 6:
        raise ValueError("Exactly 6 base str are required for a 3x2 A4 panel.")

    a4_width_inch = 8.27
    a4_height_inch = 11.69

    fig, axes = plt.subplots(3, 2, figsize=(a4_width_inch, a4_height_inch), dpi=600)
    axes = axes.ravel()

    image_records = []

    for base_str in base_str_list:
        record = plot(genome_folder_address, primary_patterns_file_address, final_results_data_path, base_str)
        image_records.append(record)

    output_folder = image_records[0]["output_folder"]

    for ax, record in zip(axes, image_records):
        img = record["panel_image"]

        ax.imshow(img)
        # ax.set_title(record["base_str"], fontsize=9, pad=3)
        ax.axis("off")

    fig.suptitle("Homo sapiens (hg38) - D3U Density", fontsize=13, y=0.99)

    fig.subplots_adjust(left=0.01, right=0.99, top=0.95, bottom=0.08, wspace=0.02, hspace=0.08)

    # -----------------------
    # Bottom legend
    # -----------------------
    # legend_ax = fig.add_axes([0.15, 0.015, 0.70, 0.045])
    legend_ax = fig.add_axes([0.18, 0.018, 0.64, 0.045])
    legend_ax.set_xticks([])
    legend_ax.set_yticks([])
    legend_ax.set_xlim(0, 1)
    legend_ax.set_ylim(0, 1)

    for spine in legend_ax.spines.values():
        spine.set_visible(True)
        spine.set_edgecolor("0.45")
        spine.set_linewidth(0.6)

    # Rectangle settings
    line_width = 3.0

    # Red legend item: STR density
    legend_ax.plot([0.08, 0.23], [0.5, 0.5], color="red", linewidth=line_width, solid_capstyle="round", transform=legend_ax.transAxes)
    legend_ax.text(0.27, 0.5, "STR Density", transform=legend_ax.transAxes, va="center", ha="left", fontsize=9)

    # Green legend item: Colony density
    legend_ax.plot([0.55, 0.70], [0.5, 0.5], color="green", linewidth=line_width, solid_capstyle="round", transform=legend_ax.transAxes)
    legend_ax.text(0.74, 0.5, "Colony Density", transform=legend_ax.transAxes, va="center", ha="left", fontsize=9)

    now = datetime.now()
    format_now = now.strftime("%Y_%m_%d_%H_%M_%S")

    panel_png_path = os.path.join(output_folder, f"genome_density_{format_now}_total.png")
    panel_jpg_path = os.path.join(output_folder, f"genome_density_{format_now}_total.jpg")

    fig.savefig(panel_png_path, dpi=600, facecolor="white")
    fig.savefig(panel_jpg_path, dpi=600, facecolor="white", pil_kwargs={"quality": 100, "subsampling": 0})

    plt.show()
    plt.close(fig)
