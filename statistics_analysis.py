import os
import re
import math
from pathlib import Path
from datetime import datetime
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

try:
    from scipy.stats import poisson as _scipy_poisson
except Exception:
    _scipy_poisson = None


DEFAULT_BASE_STR_LIST = [
    "(AC)3-(CA)3",
    "(AG)3-(GA)3",
    "(AT)3-(TA)3",
    "(CG)3-(GC)3",
    "(CT)3-(TC)3",
    "(GT)3-(TG)3",
]


def safe_filename(text: str) -> str:
    text = str(text)
    text = re.sub(r"[^A-Za-z0-9._-]+", "_", text)
    text = text.replace("-", "")
    return text.strip("_")


def standardize_chromosome_name(chromosome_name: str) -> str:
    name = str(chromosome_name).strip()
    lower_name = name.lower()

    match = re.fullmatch(r"chr0*(\d+)", lower_name)
    if match:
        number = int(match.group(1))
        if number < 10:
            return f"chr0{number}"
        return f"chr{number}"

    if lower_name in {"chrx", "x"}:
        return "chrX"
    if lower_name in {"chry", "y"}:
        return "chrY"
    if lower_name in {"chrm", "chrmt", "m", "mt"}:
        return "chrM"

    return name


def chromosome_sort_key(chromosome_name: str):
    name = standardize_chromosome_name(chromosome_name)
    lower_name = name.lower()
    match = re.fullmatch(r"chr0*(\d+)", lower_name)
    if match:
        return 0, int(match.group(1))
    if lower_name == "chrx":
        return 1, 23
    if lower_name == "chry":
        return 1, 24
    if lower_name in {"chrm", "chrmt"}:
        return 2, 25
    return 3, name


def find_latest_result_set(genome_folder_address: str) -> str:
    results_root = Path(genome_folder_address) / "Results"
    if not results_root.exists():
        raise FileNotFoundError(f"Results folder was not found: {results_root}")

    candidates = [p for p in results_root.iterdir() if p.is_dir() and p.name.startswith("Results_Set_")]
    if not candidates:
        raise FileNotFoundError(f"No Results_Set_* folder was found under: {results_root}")

    candidates.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return candidates[0].name


def resolve_result_paths(genome_folder_address: str, result_set_name: Optional[str] = None) -> Dict[str, Path]:
    if result_set_name is None:
        result_set_name = find_latest_result_set(genome_folder_address)

    result_root = Path(genome_folder_address) / "Results" / result_set_name
    paths = {
        "result_root": result_root,
        "colonies_members": result_root / "Colonies_members",
        "genome_details": result_root / "Genome_Details",
        "statistics": result_root / "Statistics_Analysis",
    }

    if not paths["result_root"].exists():
        raise FileNotFoundError(f"Result set folder was not found: {paths['result_root']}")
    if not paths["colonies_members"].exists():
        raise FileNotFoundError(f"Colonies_members folder was not found: {paths['colonies_members']}")
    if not paths["genome_details"].exists():
        raise FileNotFoundError(f"Genome_Details folder was not found: {paths['genome_details']}")

    paths["statistics"].mkdir(parents=True, exist_ok=True)
    return paths


def load_chromosome_sizes(genome_details_folder: Path, include_chr_m: bool = False) -> pd.DataFrame:
    chromosome_size_path = genome_details_folder / "Human_Genome_Chromosomes_length.csv"
    if not chromosome_size_path.exists():
        raise FileNotFoundError(f"Chromosome length file was not found: {chromosome_size_path}")

    df = pd.read_csv(chromosome_size_path, encoding="utf-8-sig")
    required_columns = {"Chromosome", "Size"}
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Chromosome length file is missing columns: {sorted(missing_columns)}")

    df = df.copy()
    df["Chromosome"] = df["Chromosome"].map(standardize_chromosome_name)
    df["Size"] = pd.to_numeric(df["Size"], errors="raise").astype(np.int64)

    if not include_chr_m:
        df = df[~df["Chromosome"].str.lower().isin(["chrm", "chrmt"])].copy()

    df = df.sort_values("Chromosome", key=lambda s: s.map(chromosome_sort_key)).reset_index(drop=True)
    return df


def build_complete_bin_table(chromosome_sizes_df: pd.DataFrame, bin_size: int = 1_000_000) -> pd.DataFrame:
    rows = []
    for row in chromosome_sizes_df.itertuples(index=False):
        chromosome = str(row.Chromosome)
        chromosome_size = int(row.Size)
        n_bins = int(math.ceil(chromosome_size / bin_size))
        for bin_index in range(n_bins):
            bin_start = bin_index * bin_size + 1
            bin_end = min((bin_index + 1) * bin_size, chromosome_size)
            rows.append(
                {
                    "Chromosome_Name": chromosome,
                    "Chromosome_Size": chromosome_size,
                    "Bin_Index": int(bin_index),
                    "Bin_Start": int(bin_start),
                    "Bin_End": int(bin_end),
                    "Bin_Length": int(bin_end - bin_start + 1),
                    "Bin_Center": int((bin_start + bin_end) // 2),
                }
            )
    return pd.DataFrame(rows)


def load_cytobands(genome_details_folder: Path) -> Optional[pd.DataFrame]:
    cytoband_path = genome_details_folder / "Human_Cytobands.tsv"
    if not cytoband_path.exists():
        return None

    cytobands = pd.read_csv(
        cytoband_path,
        sep="\t",
        header=None,
        names=["Chromosome", "ChromosomeStart", "ChromosomeEnd", "Band", "GieStain"],
        comment="#",
        dtype=str,
    )

    cytobands["Chromosome"] = cytobands["Chromosome"].astype(str).str.strip().map(standardize_chromosome_name)
    cytobands["ChromosomeStart"] = pd.to_numeric(cytobands["ChromosomeStart"], errors="coerce")
    cytobands["ChromosomeEnd"] = pd.to_numeric(cytobands["ChromosomeEnd"], errors="coerce")

    cytobands = cytobands.dropna(subset=["Chromosome", "ChromosomeStart", "ChromosomeEnd"]).copy()
    if cytobands.empty:
        return None

    cytobands["ChromosomeStart"] = cytobands["ChromosomeStart"].astype(np.int64) + 1
    cytobands["ChromosomeEnd"] = cytobands["ChromosomeEnd"].astype(np.int64)
    cytobands["Band"] = cytobands["Band"].astype(str).str.strip()
    cytobands["GieStain"] = cytobands["GieStain"].astype(str).str.strip()
    return cytobands


def get_centromere_bounds(
    chromosome: str,
    cytobands_df: Optional[pd.DataFrame],
) -> Optional[Tuple[int, int]]:
    """Return the inclusive centromere interval derived from all acen cytobands."""
    if cytobands_df is None or cytobands_df.empty:
        return None

    chromosome = standardize_chromosome_name(chromosome)
    chromosome_bands = cytobands_df[
        (cytobands_df["Chromosome"] == chromosome)
        & (cytobands_df["GieStain"].astype(str).str.lower() == "acen")
    ]
    if chromosome_bands.empty:
        return None

    centromere_start = int(chromosome_bands["ChromosomeStart"].min())
    centromere_end = int(chromosome_bands["ChromosomeEnd"].max())
    return centromere_start, centromere_end


def interval_distance(
    first_start: int,
    first_end: int,
    second_start: int,
    second_end: int,
) -> int:
    """Return zero for overlapping intervals, otherwise the edge-to-edge distance."""
    if first_end < second_start:
        return int(second_start - first_end)
    if second_end < first_start:
        return int(first_start - second_end)
    return 0


def get_peak_region_annotation(
    chromosome: str,
    bin_start: int,
    bin_end: int,
    chromosome_size: int,
    cytobands_df: Optional[pd.DataFrame],
    telomere_radius_bp: int = 1_000_000,
    pericentromere_radius_bp: int = 2_500_000,
) -> Dict[str, object]:
    """Annotate a peak interval relative to telomeres and the cytoband centromere."""
    chromosome = standardize_chromosome_name(chromosome)
    bin_start = int(bin_start)
    bin_end = int(bin_end)
    chromosome_size = int(chromosome_size)

    if bin_start < 1 or bin_end < bin_start or bin_end > chromosome_size:
        raise ValueError(
            f"Invalid peak interval for {chromosome}: {bin_start}-{bin_end} "
            f"with chromosome size {chromosome_size}"
        )
    if telomere_radius_bp < 0:
        raise ValueError("telomere_radius_bp must be greater than or equal to zero")
    if pericentromere_radius_bp < 0:
        raise ValueError("pericentromere_radius_bp must be greater than or equal to zero")

    peak_center = int((bin_start + bin_end) // 2)
    distance_from_chromosome_start = int(bin_start - 1)
    distance_from_chromosome_end = int(chromosome_size - bin_end)
    distance_to_p_telomere = distance_from_chromosome_start
    distance_to_q_telomere = distance_from_chromosome_end

    if distance_to_p_telomere <= distance_to_q_telomere:
        nearest_telomere = "p_telomere"
        distance_to_nearest_telomere = distance_to_p_telomere
    else:
        nearest_telomere = "q_telomere"
        distance_to_nearest_telomere = distance_to_q_telomere

    centromere_bounds = get_centromere_bounds(chromosome, cytobands_df)
    centromere_available = centromere_bounds is not None

    centromere_start = float("nan")
    centromere_end = float("nan")
    centromere_center = float("nan")
    distance_to_centromere = float("nan")
    distance_center_to_centromere_center = float("nan")
    signed_center_position_from_centromere = float("nan")
    nearest_centromere_boundary = "unavailable"
    chromosome_arm = "unknown"
    overlaps_centromere = False

    if centromere_bounds is not None:
        centromere_start, centromere_end = centromere_bounds
        centromere_center = int((centromere_start + centromere_end) // 2)
        distance_to_centromere = interval_distance(
            bin_start,
            bin_end,
            centromere_start,
            centromere_end,
        )
        distance_center_to_centromere_center = int(abs(peak_center - centromere_center))
        signed_center_position_from_centromere = int(peak_center - centromere_center)

        if bin_end < centromere_start:
            chromosome_arm = "p"
            nearest_centromere_boundary = "p_side_boundary"
        elif bin_start > centromere_end:
            chromosome_arm = "q"
            nearest_centromere_boundary = "q_side_boundary"
        else:
            chromosome_arm = "centromeric"
            nearest_centromere_boundary = "overlap"
            overlaps_centromere = True

    if overlaps_centromere:
        region_class = "centromeric"
    elif distance_to_p_telomere <= telomere_radius_bp:
        region_class = "p_terminal"
    elif distance_to_q_telomere <= telomere_radius_bp:
        region_class = "q_terminal"
    elif (
        centromere_available
        and chromosome_arm == "p"
        and distance_to_centromere <= pericentromere_radius_bp
    ):
        region_class = "p_pericentromeric"
    elif (
        centromere_available
        and chromosome_arm == "q"
        and distance_to_centromere <= pericentromere_radius_bp
    ):
        region_class = "q_pericentromeric"
    else:
        region_class = "internal"

    landmark_distances = {
        "p_telomere": distance_to_p_telomere,
        "q_telomere": distance_to_q_telomere,
    }
    if centromere_available:
        landmark_distances["centromere"] = int(distance_to_centromere)

    nearest_landmark = min(landmark_distances, key=landmark_distances.get)
    distance_to_nearest_landmark = int(landmark_distances[nearest_landmark])

    return {
        "Region_Class": region_class,
        "Peak_Chromosome_Arm": chromosome_arm,
        "Peak_Center_Position": peak_center,
        "Peak_Center_Percent_Chromosome": float((peak_center / chromosome_size) * 100.0),
        "Distance_From_Chromosome_Start_Bp": distance_from_chromosome_start,
        "Distance_From_Chromosome_End_Bp": distance_from_chromosome_end,
        "Distance_To_P_Telomere_Bp": distance_to_p_telomere,
        "Distance_To_Q_Telomere_Bp": distance_to_q_telomere,
        "Nearest_Telomere": nearest_telomere,
        "Distance_To_Nearest_Telomere_Bp": int(distance_to_nearest_telomere),
        "Centromere_Annotation_Available": bool(centromere_available),
        "Centromere_Start": centromere_start,
        "Centromere_End": centromere_end,
        "Centromere_Center": centromere_center,
        "Peak_Overlaps_Centromere": bool(overlaps_centromere),
        "Nearest_Centromere_Boundary": nearest_centromere_boundary,
        "Distance_To_Centromere_Bp": distance_to_centromere,
        "Distance_Peak_Center_To_Centromere_Center_Bp": distance_center_to_centromere_center,
        "Signed_Peak_Center_From_Centromere_Center_Bp": signed_center_position_from_centromere,
        "Nearest_Chromosomal_Landmark": nearest_landmark,
        "Distance_To_Nearest_Chromosomal_Landmark_Bp": distance_to_nearest_landmark,
        "Telomere_Radius_Bp": int(telomere_radius_bp),
        "Pericentromere_Radius_Bp": int(pericentromere_radius_bp),
    }


def annotate_region_class(
    chromosome: str,
    bin_start: int,
    bin_end: int,
    chromosome_size: int,
    cytobands_df: Optional[pd.DataFrame],
    telomere_radius_bp: int = 1_000_000,
    pericentromere_radius_bp: int = 5_000_000,
) -> str:
    """Return only the coarse region class; retained for backward compatibility."""
    annotation = get_peak_region_annotation(
        chromosome=chromosome,
        bin_start=bin_start,
        bin_end=bin_end,
        chromosome_size=chromosome_size,
        cytobands_df=cytobands_df,
        telomere_radius_bp=telomere_radius_bp,
        pericentromere_radius_bp=pericentromere_radius_bp,
    )
    return str(annotation["Region_Class"])


def load_integrated_d3u_map(colonies_members_folder: Path, base_str: str) -> pd.DataFrame:
    input_path = colonies_members_folder / f"Human_Integrated_D3Us_{base_str}_Map_Result.csv"
    if not input_path.exists():
        raise FileNotFoundError(f"Integrated D3U map file was not found: {input_path}")

    df = pd.read_csv(input_path)
    required_columns = {
        "Chromosome_Name",
        "Start_Locus",
        "End_Locus",
        "Core",
        "Colony_Index",
        "Colony_Member_Count",
        "Colony_Type",
    }
    missing_columns = required_columns.difference(df.columns)
    if missing_columns:
        raise ValueError(f"Input file for {base_str} is missing columns: {sorted(missing_columns)}")

    df = df.copy()
    df["Chromosome_Name"] = df["Chromosome_Name"].map(standardize_chromosome_name)
    df["Start_Locus"] = pd.to_numeric(df["Start_Locus"], errors="raise").astype(np.int64)
    df["End_Locus"] = pd.to_numeric(df["End_Locus"], errors="raise").astype(np.int64)
    df["Colony_Index"] = pd.to_numeric(df["Colony_Index"], errors="raise").astype(np.int64)
    df["Colony_Member_Count"] = pd.to_numeric(df["Colony_Member_Count"], errors="raise").astype(np.int64)
    return df


def build_colony_map(d3u_map_df: pd.DataFrame, min_colony_members: int = 3) -> pd.DataFrame:
    if d3u_map_df.empty:
        return pd.DataFrame(
            columns=[
                "Chromosome_Name",
                "Colony_Index",
                "Colony_Member_Count",
                "Colony_Type",
                "member_count_check",
                "Colony_Start",
                "Colony_End",
                "Colony_Length",
                "Colony_Bin_Median",
            ]
        )

    grouped = d3u_map_df.groupby(
        ["Chromosome_Name", "Colony_Index", "Colony_Member_Count", "Colony_Type"],
        as_index=False,
    )
    colony_map = grouped.agg(
        member_count_check=("Core", "size"),
        Colony_Start=("Start_Locus", "min"),
        Colony_End=("End_Locus", "max"),
    )
    colony_map["Colony_Length"] = colony_map["Colony_End"] - colony_map["Colony_Start"] + 1
    colony_map["Colony_Bin_Median"] = ((colony_map["Colony_Start"] + colony_map["Colony_End"]) // 2).astype(np.int64)
    colony_map = colony_map[colony_map["Colony_Member_Count"] >= min_colony_members].copy()
    return colony_map


def build_pair_bin_counts(
    d3u_map_df: pd.DataFrame,
    base_bins_df: pd.DataFrame,
    base_str: str,
    bin_size: int = 1_000_000,
    min_colony_members: int = 3,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    d3u_df = d3u_map_df.copy()
    d3u_df["D3U_Bin_Median"] = ((d3u_df["Start_Locus"] + d3u_df["End_Locus"]) // 2).astype(np.int64)
    d3u_df["Bin_Index"] = ((d3u_df["D3U_Bin_Median"] - 1) // bin_size).astype(np.int64)

    d3u_counts = (
        d3u_df.groupby(["Chromosome_Name", "Bin_Index"])
        .size()
        .reset_index(name="D3U_Count")
    )

    colony_map = build_colony_map(d3u_map_df, min_colony_members=min_colony_members)
    if colony_map.empty:
        colony_counts = pd.DataFrame(columns=["Chromosome_Name", "Bin_Index", "Colony_Count"])
    else:
        colony_map["Bin_Index"] = ((colony_map["Colony_Bin_Median"] - 1) // bin_size).astype(np.int64)
        colony_counts = (
            colony_map.groupby(["Chromosome_Name", "Bin_Index"])
            .size()
            .reset_index(name="Colony_Count")
        )

    pair_bins = base_bins_df.merge(d3u_counts, on=["Chromosome_Name", "Bin_Index"], how="left")
    pair_bins = pair_bins.merge(colony_counts, on=["Chromosome_Name", "Bin_Index"], how="left")
    pair_bins["D3U_Count"] = pair_bins["D3U_Count"].fillna(0).astype(np.int64)
    pair_bins["Colony_Count"] = pair_bins["Colony_Count"].fillna(0).astype(np.int64)
    pair_bins["D3U_Density_Per_Mb"] = pair_bins["D3U_Count"] * (1_000_000 / pair_bins["Bin_Length"])
    pair_bins["Colony_Density_Per_Mb"] = pair_bins["Colony_Count"] * (1_000_000 / pair_bins["Bin_Length"])
    pair_bins["D3U_Pair"] = base_str

    ordered_columns = [
        "D3U_Pair",
        "Chromosome_Name",
        "Chromosome_Size",
        "Bin_Index",
        "Bin_Start",
        "Bin_End",
        "Bin_Length",
        "Bin_Center",
        "D3U_Count",
        "D3U_Density_Per_Mb",
        "Colony_Count",
        "Colony_Density_Per_Mb",
    ]
    return pair_bins[ordered_columns].copy(), colony_map


def poisson_upper_tail(observed_count: int, expected_lambda: float) -> float:
    observed_count = int(observed_count)
    expected_lambda = float(expected_lambda)

    if observed_count <= 0:
        return 1.0
    if expected_lambda < 0:
        return float("nan")
    if expected_lambda == 0:
        return 0.0

    if _scipy_poisson is not None:
        return float(_scipy_poisson.sf(observed_count - 1, expected_lambda))

    max_i = observed_count - 1
    if max_i < 100000 and expected_lambda < 10000:
        term = math.exp(-expected_lambda)
        cdf = term
        for i in range(1, max_i + 1):
            term *= expected_lambda / i
            cdf += term
        return float(min(max(1.0 - cdf, 0.0), 1.0))

    z = (observed_count - 0.5 - expected_lambda) / math.sqrt(expected_lambda)
    return float(0.5 * math.erfc(z / math.sqrt(2.0)))


def empirical_upper_tail(observed_density: float, background_density: pd.Series) -> float:
    clean_bg = background_density.dropna()
    n = int(clean_bg.shape[0])
    if n == 0:
        return float("nan")
    return float((1 + (clean_bg >= observed_density).sum()) / (1 + n))


def benjamini_hochberg(p_values: Iterable[float]) -> np.ndarray:
    p = np.asarray(list(p_values), dtype=float)
    q = np.full(p.shape, np.nan, dtype=float)
    finite_mask = np.isfinite(p)
    if finite_mask.sum() == 0:
        return q

    finite_p = p[finite_mask]
    order = np.argsort(finite_p)
    ranked = finite_p[order]
    n = ranked.shape[0]
    raw_q = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(raw_q[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)

    finite_q = np.empty_like(adjusted)
    finite_q[order] = adjusted
    q[finite_mask] = finite_q
    return q


def summarize_pair(pair_bins_df: pd.DataFrame, colony_map_df: pd.DataFrame, base_str: str) -> Dict[str, object]:
    d3u_peak = pair_bins_df.sort_values(
        ["D3U_Density_Per_Mb", "D3U_Count", "Chromosome_Name", "Bin_Index"],
        ascending=[False, False, True, True],
    ).iloc[0]

    colony_peak = pair_bins_df.sort_values(
        ["Colony_Density_Per_Mb", "Colony_Count", "Chromosome_Name", "Bin_Index"],
        ascending=[False, False, True, True],
    ).iloc[0]

    return {
        "D3U_Pair": base_str,
        "Total_D3U_Count": int(pair_bins_df["D3U_Count"].sum()),
        "Total_Colony_Count": int(colony_map_df.shape[0]),
        "Median_D3U_Density_Per_Mb": float(pair_bins_df["D3U_Density_Per_Mb"].median()),
        "Mean_D3U_Density_Per_Mb": float(pair_bins_df["D3U_Density_Per_Mb"].mean()),
        "Max_D3U_Density_Per_Mb": float(d3u_peak["D3U_Density_Per_Mb"]),
        "Max_D3U_Count": int(d3u_peak["D3U_Count"]),
        "Max_D3U_Chromosome": d3u_peak["Chromosome_Name"],
        "Max_D3U_Bin_Index": int(d3u_peak["Bin_Index"]),
        "Max_D3U_Bin_Start": int(d3u_peak["Bin_Start"]),
        "Max_D3U_Bin_End": int(d3u_peak["Bin_End"]),
        "Median_Colony_Density_Per_Mb": float(pair_bins_df["Colony_Density_Per_Mb"].median()),
        "Mean_Colony_Density_Per_Mb": float(pair_bins_df["Colony_Density_Per_Mb"].mean()),
        "Max_Colony_Density_Per_Mb": float(colony_peak["Colony_Density_Per_Mb"]),
        "Max_Colony_Count": int(colony_peak["Colony_Count"]),
        "Max_Colony_Chromosome": colony_peak["Chromosome_Name"],
        "Max_Colony_Bin_Index": int(colony_peak["Bin_Index"]),
        "Max_Colony_Bin_Start": int(colony_peak["Bin_Start"]),
        "Max_Colony_Bin_End": int(colony_peak["Bin_End"]),
    }


def calculate_peak_statistics_for_signal(
    pair_bins_df: pd.DataFrame,
    base_str: str,
    signal_type: str,
    count_column: str,
    density_column: str,
    local_radius_bins: int = 5,
    cytobands_df: Optional[pd.DataFrame] = None,
    telomere_radius_bp: int = 1_000_000,
    pericentromere_radius_bp: int = 5_000_000,
) -> List[Dict[str, object]]:
    records = []
    chromosomes = sorted(pair_bins_df["Chromosome_Name"].unique().tolist(), key=chromosome_sort_key)

    for chromosome in chromosomes:
        chr_df = pair_bins_df[pair_bins_df["Chromosome_Name"] == chromosome].copy()
        chr_df = chr_df.sort_values("Bin_Index").reset_index(drop=True)
        if chr_df.empty:
            continue

        peak_row = chr_df.sort_values(
            [density_column, count_column, "Bin_Index"],
            ascending=[False, False, True],
        ).iloc[0]

        peak_bin_index = int(peak_row["Bin_Index"])
        background_mask = (
            (chr_df["Bin_Index"] >= peak_bin_index - local_radius_bins)
            & (chr_df["Bin_Index"] <= peak_bin_index + local_radius_bins)
            & (chr_df["Bin_Index"] != peak_bin_index)
        )
        background_df = chr_df[background_mask].copy()

        observed_count = int(peak_row[count_column])
        observed_density = float(peak_row[density_column])
        target_bin_length = int(peak_row["Bin_Length"])

        background_total_count = int(background_df[count_column].sum()) if not background_df.empty else 0
        background_total_bp = int(background_df["Bin_Length"].sum()) if not background_df.empty else 0
        background_n_bins = int(background_df.shape[0])

        if background_total_bp > 0:
            local_rate_per_bp = background_total_count / background_total_bp
            expected_lambda = local_rate_per_bp * target_bin_length
            local_mean_density_per_mb = local_rate_per_bp * 1_000_000
            local_mean_count_per_bin = background_total_count / background_n_bins if background_n_bins > 0 else float("nan")
            local_median_count_per_bin = float(background_df[count_column].median()) if background_n_bins > 0 else float("nan")
            local_median_density_per_mb = float(background_df[density_column].median()) if background_n_bins > 0 else float("nan")
            poisson_p = poisson_upper_tail(observed_count, expected_lambda)
            empirical_p = empirical_upper_tail(observed_density, background_df[density_column])
        else:
            expected_lambda = 0.0
            local_mean_density_per_mb = 0.0
            local_mean_count_per_bin = float("nan")
            local_median_count_per_bin = float("nan")
            local_median_density_per_mb = float("nan")
            poisson_p = poisson_upper_tail(observed_count, expected_lambda)
            empirical_p = float("nan")

        fold_over_expected = (observed_count + 1.0) / (expected_lambda + 1.0)
        log2_fold_over_expected = float(np.log2(fold_over_expected))
        density_fold_over_local_mean = (observed_density + 1.0) / (local_mean_density_per_mb + 1.0)
        log2_density_fold_over_local_mean = float(np.log2(density_fold_over_local_mean))

        chromosome_median_count = float(chr_df[count_column].median())
        chromosome_mean_count = float(chr_df[count_column].mean())
        chromosome_median_density = float(chr_df[density_column].median())
        chromosome_mean_density = float(chr_df[density_column].mean())

        region_annotation = get_peak_region_annotation(
            chromosome=chromosome,
            bin_start=int(peak_row["Bin_Start"]),
            bin_end=int(peak_row["Bin_End"]),
            chromosome_size=int(peak_row["Chromosome_Size"]),
            cytobands_df=cytobands_df,
            telomere_radius_bp=telomere_radius_bp,
            pericentromere_radius_bp=pericentromere_radius_bp,
        )

        records.append(
            {
                "D3U_Pair": base_str,
                "Signal_Type": signal_type,
                "Chromosome": chromosome,
                **region_annotation,
                "Peak_Bin_Index": peak_bin_index,
                "Peak_Bin_Start": int(peak_row["Bin_Start"]),
                "Peak_Bin_End": int(peak_row["Bin_End"]),
                "Peak_Bin_Length": target_bin_length,
                "Observed_Count": observed_count,
                "Observed_Density_Per_Mb": observed_density,
                "Expected_Count_Local_Lambda": float(expected_lambda),
                "Local_Background_Radius_Bins": int(local_radius_bins),
                "Local_Background_N_Bins": background_n_bins,
                "Local_Background_Total_Count": background_total_count,
                "Local_Background_Total_Bp": background_total_bp,
                "Local_Mean_Count_Per_Bin": float(local_mean_count_per_bin),
                "Local_Median_Count_Per_Bin": float(local_median_count_per_bin),
                "Local_Mean_Density_Per_Mb": float(local_mean_density_per_mb),
                "Local_Median_Density_Per_Mb": float(local_median_density_per_mb),
                "Fold_Over_Expected_Count_Plus1": float(fold_over_expected),
                "Log2_Fold_Over_Expected_Count_Plus1": log2_fold_over_expected,
                "Fold_Over_Local_Mean_Density_Plus1": float(density_fold_over_local_mean),
                "Log2_Fold_Over_Local_Mean_Density_Plus1": log2_density_fold_over_local_mean,
                "Chromosome_Median_Count_Per_Bin": chromosome_median_count,
                "Chromosome_Mean_Count_Per_Bin": chromosome_mean_count,
                "Chromosome_Median_Density_Per_Mb": chromosome_median_density,
                "Chromosome_Mean_Density_Per_Mb": chromosome_mean_density,
                "Poisson_Local_One_Sided_P": float(poisson_p),
                "Empirical_Local_One_Sided_P": float(empirical_p),
            }
        )

    return records


def add_fdr_columns(peak_stats_df: pd.DataFrame) -> pd.DataFrame:
    df = peak_stats_df.copy()
    df["Poisson_Local_Q_Global_BH"] = benjamini_hochberg(df["Poisson_Local_One_Sided_P"].to_numpy())
    df["Empirical_Local_Q_Global_BH"] = benjamini_hochberg(df["Empirical_Local_One_Sided_P"].to_numpy())

    df["Poisson_Local_Q_By_Signal_BH"] = np.nan
    df["Empirical_Local_Q_By_Signal_BH"] = np.nan
    for signal_type, idx in df.groupby("Signal_Type").groups.items():
        signal_idx = list(idx)
        df.loc[signal_idx, "Poisson_Local_Q_By_Signal_BH"] = benjamini_hochberg(
            df.loc[signal_idx, "Poisson_Local_One_Sided_P"].to_numpy()
        )
        df.loc[signal_idx, "Empirical_Local_Q_By_Signal_BH"] = benjamini_hochberg(
            df.loc[signal_idx, "Empirical_Local_One_Sided_P"].to_numpy()
        )

    return df


def write_readme(
    output_folder: Path,
    result_set_name: str,
    bin_size: int,
    local_radius_bins: int,
    min_colony_members: int,
    telomere_radius_bp: int,
    pericentromere_radius_bp: int,
    alpha: float,
):
    readme_path = output_folder / "README_Statistics_Analysis.txt"
    content = f"""Statistics analysis output
Generated at: {datetime.now().isoformat(timespec='seconds')}
Result set: {result_set_name}
Bin size: {bin_size} bp
Local background radius: +/- {local_radius_bins} bins, excluding the peak bin
Minimum colony members: {min_colony_members}
Telomere terminal radius for annotation: {telomere_radius_bp} bp
Pericentromeric radius for annotation: {pericentromere_radius_bp} bp
Significance threshold: Poisson_Local_Q_Global_BH <= {alpha}

Main files:
1. All_Pairs_Complete_1Mb_Bin_Counts.csv
   Complete bin-level table with zero-count bins included.

2. Pair_Summary.csv
   One-row summary per D3U pair.

3. Per_Chromosome_Peak_Statistics.csv
   For each D3U pair, chromosome, and signal type, this table reports the highest-density bin and its local enrichment statistics.

4. Significant_Peaks_Poisson_FDR_{str(alpha).replace('.', '_')}.csv
   Subset of peak statistics with Poisson_Local_Q_Global_BH <= alpha.

Notes:
- Poisson_Local_One_Sided_P is P(X >= observed count) using lambda estimated from local flanking bins.
- Empirical_Local_One_Sided_P is computed from the local background bins as (1 + number of background bins with density >= observed density) / (1 + number of background bins).
- FDR columns are Benjamini-Hochberg adjusted q-values.
- Region_Class uses the following mutually exclusive priority: centromeric, p_terminal/q_terminal, p_pericentromeric/q_pericentromeric, internal.
- p_pericentromeric and q_pericentromeric are assigned when a peak lies on the corresponding chromosome arm and its edge-to-edge distance from the acen centromere interval is within the configured pericentromeric radius.
- Distance_To_Centromere_Bp is zero for a centromere-overlapping peak; otherwise it is the shortest edge-to-edge distance from the peak bin to the acen centromere interval.
- Distance_To_P_Telomere_Bp and Distance_To_Q_Telomere_Bp are measured from the nearest peak-bin edge to the corresponding chromosome end.
- Nearest_Chromosomal_Landmark identifies the closest of the p telomere, q telomere, and centromere.
- Cytoband files with or without a header row are supported.
"""
    readme_path.write_text(content, encoding="utf-8")


def analyze(
    genome_folder_address: str,
    base_str_list: Optional[List[str]] = None,
    result_set_name: Optional[str] = None,
    bin_size: int = 1_000_000,
    local_radius_bins: int = 5,
    min_colony_members: int = 3,
    telomere_radius_bp: int = 1_000_000,
    pericentromere_radius_bp: int = 5_000_000,
    alpha: float = 0.05,
    include_chr_m: bool = False,
    write_per_pair_bin_tables: bool = True,
) -> Dict[str, str]:
    if base_str_list is None:
        base_str_list = DEFAULT_BASE_STR_LIST

    paths = resolve_result_paths(genome_folder_address, result_set_name=result_set_name)
    resolved_result_set_name = paths["result_root"].name

    chromosome_sizes_df = load_chromosome_sizes(paths["genome_details"], include_chr_m=include_chr_m)
    base_bins_df = build_complete_bin_table(chromosome_sizes_df, bin_size=bin_size)
    cytobands_df = load_cytobands(paths["genome_details"])

    all_pair_bins = []
    pair_summaries = []
    peak_records = []

    per_pair_folder = paths["statistics"] / "Per_Pair_Bin_Tables"
    if write_per_pair_bin_tables:
        per_pair_folder.mkdir(parents=True, exist_ok=True)

    for base_str in base_str_list:
        d3u_map_df = load_integrated_d3u_map(paths["colonies_members"], base_str=base_str)
        pair_bins_df, colony_map_df = build_pair_bin_counts(
            d3u_map_df=d3u_map_df,
            base_bins_df=base_bins_df,
            base_str=base_str,
            bin_size=bin_size,
            min_colony_members=min_colony_members,
        )

        pair_summaries.append(summarize_pair(pair_bins_df, colony_map_df, base_str=base_str))

        peak_records.extend(
            calculate_peak_statistics_for_signal(
                pair_bins_df=pair_bins_df,
                base_str=base_str,
                signal_type="D3U",
                count_column="D3U_Count",
                density_column="D3U_Density_Per_Mb",
                local_radius_bins=local_radius_bins,
                cytobands_df=cytobands_df,
                telomere_radius_bp=telomere_radius_bp,
                pericentromere_radius_bp=pericentromere_radius_bp,
            )
        )
        peak_records.extend(
            calculate_peak_statistics_for_signal(
                pair_bins_df=pair_bins_df,
                base_str=base_str,
                signal_type="Colony",
                count_column="Colony_Count",
                density_column="Colony_Density_Per_Mb",
                local_radius_bins=local_radius_bins,
                cytobands_df=cytobands_df,
                telomere_radius_bp=telomere_radius_bp,
                pericentromere_radius_bp=pericentromere_radius_bp,
            )
        )

        if write_per_pair_bin_tables:
            per_pair_path = per_pair_folder / f"Complete_1Mb_Bin_Counts_{safe_filename(base_str)}.csv"
            pair_bins_df.to_csv(per_pair_path, index=False, encoding="utf-8-sig")

        all_pair_bins.append(pair_bins_df)

    all_pair_bins_df = pd.concat(all_pair_bins, axis=0, ignore_index=True)
    pair_summary_df = pd.DataFrame(pair_summaries)
    peak_stats_df = pd.DataFrame(peak_records)
    peak_stats_df = add_fdr_columns(peak_stats_df)

    pair_summary_df = pair_summary_df.sort_values("D3U_Pair").reset_index(drop=True)
    peak_stats_df["_Chromosome_Sort_Key"] = peak_stats_df["Chromosome"].map(
        lambda x: f"{chromosome_sort_key(x)[0]:02d}_{chromosome_sort_key(x)[1]:02d}"
    )
    peak_stats_df = peak_stats_df.sort_values(
        ["D3U_Pair", "Signal_Type", "_Chromosome_Sort_Key"]
    ).drop(columns=["_Chromosome_Sort_Key"]).reset_index(drop=True)

    all_bins_path = paths["statistics"] / "All_Pairs_Complete_1Mb_Bin_Counts.csv"
    summary_path = paths["statistics"] / "Pair_Summary.csv"
    peak_stats_path = paths["statistics"] / "Per_Chromosome_Peak_Statistics.csv"
    significant_path = paths["statistics"] / f"Significant_Peaks_Poisson_FDR_{str(alpha).replace('.', '_')}.csv"

    significant_df = peak_stats_df[peak_stats_df["Poisson_Local_Q_Global_BH"] <= alpha].copy()
    significant_df = significant_df.sort_values(
        ["Poisson_Local_Q_Global_BH", "Poisson_Local_One_Sided_P", "D3U_Pair", "Signal_Type", "Chromosome"]
    ).reset_index(drop=True)

    all_pair_bins_df.to_csv(all_bins_path, index=False, encoding="utf-8-sig")
    pair_summary_df.to_csv(summary_path, index=False, encoding="utf-8-sig")
    peak_stats_df.to_csv(peak_stats_path, index=False, encoding="utf-8-sig")
    significant_df.to_csv(significant_path, index=False, encoding="utf-8-sig")

    write_readme(
        output_folder=paths["statistics"],
        result_set_name=resolved_result_set_name,
        bin_size=bin_size,
        local_radius_bins=local_radius_bins,
        min_colony_members=min_colony_members,
        telomere_radius_bp=telomere_radius_bp,
        pericentromere_radius_bp=pericentromere_radius_bp,
        alpha=alpha,
    )

    return {
        "result_set": resolved_result_set_name,
        "statistics_folder": str(paths["statistics"]),
        "all_bins_csv": str(all_bins_path),
        "pair_summary_csv": str(summary_path),
        "peak_statistics_csv": str(peak_stats_path),
        "significant_peaks_csv": str(significant_path),
    }
