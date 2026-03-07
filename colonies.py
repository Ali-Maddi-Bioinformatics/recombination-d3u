import os
import pandas as pd
from tqdm import tqdm
import logging
import re
from datetime import datetime
from pathlib import Path
import pickle

import utilities


def reverse_str(target_str):
    pattern = r"\((.*?)\)(\d+)"
    match = re.match(pattern, target_str)

    if match:
        core = match.group(1)
        repeat = int(match.group(2))
        return F"({core[::-1]}){repeat}"
    raise ValueError("Invalid STR string...")


def calculate_colonies(dataset, key, chromosome_name):
    colonies = {}
    if dataset.empty:
        return colonies

    dataset.sort_values(by="Start_Locus", ascending=True, inplace=True)
    dataset.reset_index(drop=True, inplace=True)
    dataset["Distance_from_previous"] = dataset["Start_Locus"] - dataset["End_Locus"].shift(1)
    dataset.loc[0, "Distance_from_previous"] = int(0)
    dataset["Distance_from_previous"] = dataset["Distance_from_previous"].astype("int64")
    dataset["Distance_to_next"] = dataset["Start_Locus"].shift(-1) - dataset["End_Locus"]
    dataset.at[dataset.index[-1], "Distance_to_next"] = int(0)
    dataset["Distance_to_next"] = dataset["Distance_to_next"].astype("int64")

    colony_index = 1
    colonies_member_count_dict = {}
    last_seen_colony_member_count = 0

    for index, str in dataset.iterrows():
        dataset.at[index, "Colony_Index"] = colony_index
        last_seen_colony_member_count += 1

        if str["Distance_to_next"] > 500:
            colonies[last_seen_colony_member_count] = colonies.get(last_seen_colony_member_count, 0) + 1
            colonies_member_count_dict[colony_index] = last_seen_colony_member_count
            colony_index += 1
            last_seen_colony_member_count = 0


    colonies[last_seen_colony_member_count] = colonies.get(last_seen_colony_member_count, 0) + 1
    colonies_member_count_dict[colony_index] = last_seen_colony_member_count

    for index in dataset.index:
        colony_index = dataset.at[index, "Colony_Index"]
        dataset.at[index, "Colony_Member_Count"] = colonies_member_count_dict[colony_index]
        dataset.at[index, "Colony_Type"] = f"C({colonies_member_count_dict[colony_index]})"

    colonies_members_total_count = 0
    for colony_type_key in colonies.keys():
        colonies_members_total_count += colonies.get(colony_type_key) * colony_type_key
    if colonies_members_total_count != dataset.shape[0]:
        print(F"###### Error: The counting member of colonies has problem! We count total {colonies_members_total_count} {key} STRs in "
              F"chromosome {chromosome_name} while there is {dataset.shape[0]} STRs in it ######")
        raise ValueError(F"###### Error: The counting member of colonies has problem! We count total {colonies_members_total_count} {key} STRs in "
              F"chromosome {chromosome_name} while there is {dataset.shape[0]} STRs in it ######")
    return colonies


def process(genome_folder_address, primary_patterns_file_address):
    logging.basicConfig(level=logging.INFO)

    # chromosomes = {}
    total_colonies_results = {}
    bar_format = (
        "{desc:45}: {percentage:3.0f}%|"
        "{bar:30}| "
        "{n}/{total} "
        "[{elapsed}<{remaining}, {rate_fmt}, {postfix}]"
    )

    print("Starting scan chromosomes files...")

    primary_chromosome_dict = utilities.read_chromosomes(genome_folder_address)
    patterns_list = utilities.read_patterns(primary_patterns_file_address)

    patterns_map_dict = {}
    progressbar_0 = tqdm(primary_chromosome_dict.items(), desc="Scanning the chromosomes", bar_format=bar_format, dynamic_ncols=False)
    for primary_chromosome_name, dna_sequence in progressbar_0:
        occurs_list = []
        for pattern in patterns_list:
            # scan_start_index = 0
            # occur_index = chromosome_dna.find(pattern["extend"], 0)
            # while occur_index != -1:
            #     scan_start_index = occur_index + len(pattern["extend"])
            #     occurs_list_1.append({
            #         "Core": pattern["core"],
            #         "Repeat": pattern["repeat"],
            #         "STR": pattern["pack"],
            #         "Sequence": F"{chromosome_dna[occur_index:scan_start_index]}",
            #         "Start_Locus": occur_index,
            #         "End_Locus": scan_start_index
            #     })
            #     occur_index = chromosome_dna.find(pattern["extend"], scan_start_index)
            # previous_location = 0
            chromosome_number = primary_chromosome_name.split("_")[0].split("chr")[1]
            matches = re.finditer(pattern["extend"], dna_sequence, re.IGNORECASE)
            for match in matches:
                occurs_list.append({
                    "Core": match.group(0)[0:len(pattern["core"])],
                    "Repeat": pattern["repeat"],
                    "STR": pattern["pack"],
                    "Sequence": match.group(0),
                    "Start_Locus": match.start() + 1,
                    "End_Locus": match.end() + 1,
                    "Chromosome_Name": utilities.standard_name(primary_chromosome_name.split("_")[0]),
                    "Link_to_Location":
                        F"https://asia.ensembl.org/Homo_sapiens/Location/View?r={chromosome_number}:{match.start() + 1}-{match.end() + 1};db=core",
                    "Distance_from_previous": "",
                    "Distance_to_next": "",
                    "Colony_Index": "",
                    "Colony_Member_Count": 0,
                    "Colony_Type": ""
                })
                # previous_location = match.start()

        chromosome_df = pd.DataFrame(occurs_list, columns=["Core", "Repeat", "STR", "Sequence", "Start_Locus", "End_Locus",
                                                           "Chromosome_Name", "Link_to_Location", "Distance_from_previous", "Distance_to_next",
                                                           "Colony_Index", "Colony_Member_Count", "Colony_Type"])
        chromosome_df.sort_values(by="Start_Locus", ascending=True, inplace=True)
        patterns_map_dict[primary_chromosome_name] = chromosome_df

    print("Scanning chromosomes files were finished.")
    now = datetime.now()
    format_now = now.strftime("%Y_%m_%d_%H_%M_%S")
    folder_of_current_results = F"Results_Set_{format_now}"
    folder_of_current_results = genome_folder_address + os.sep + "Results" + os.sep + folder_of_current_results
    folder_of_current_strs_map_result = folder_of_current_results + os.sep + "STRs_Map"
    folder_of_current_colonies_result = folder_of_current_results + os.sep + "Colonies"
    folder_of_current_colonies_members_result = folder_of_current_results + os.sep + "Colonies_members"
    folder_of_current_data_structures = folder_of_current_results + os.sep + "Data_Structures"
    folder_of_current_genome_details = folder_of_current_results + os.sep + "Genome_Details"
    Path(folder_of_current_strs_map_result).mkdir(parents=True, exist_ok=True)
    Path(folder_of_current_colonies_result).mkdir(parents=True, exist_ok=True)
    Path(folder_of_current_colonies_members_result).mkdir(parents=True, exist_ok=True)
    Path(folder_of_current_data_structures).mkdir(parents=True, exist_ok=True)
    Path(folder_of_current_genome_details).mkdir(parents=True, exist_ok=True)

    chromosome_unique_names = set()
    for primary_chromosome_name, df_str_map in patterns_map_dict.items():
        df_str_map.to_csv(folder_of_current_strs_map_result + os.sep + primary_chromosome_name + ".csv", index=False)
        chromosome_unique_names.add(utilities.standard_name(primary_chromosome_name.split("_")[0]))

    print("Writing STRs Map files were finished.")

    if len(chromosome_unique_names) != len(patterns_map_dict):
        print("Duplicate chromosome name was detected.")
        raise ValueError("Duplicate chromosome name was found ...")

    print("All chromosome names are unique.")
    progressbar_1 = tqdm(patterns_map_dict.items(), desc="processing the STRs Map files", bar_format=bar_format, dynamic_ncols=False)
    for primary_chromosome_name, df_str_map in progressbar_1:
        progressbar_1.set_postfix_str(F"processing file {primary_chromosome_name}")
        coreLength = 2
        total_colonies_results.setdefault(coreLength, {})
        chromosome_name = utilities.standard_name(primary_chromosome_name.split("_")[0])

        progressbar_1.clear()

        logging.info(F"\n\nReading STRs from: {primary_chromosome_name}")
        # str_map_dataframe = pd.read_csv(strs_file_path, sep="\t", header=None)
        # str_map_dataframe.columns = ["Core", "Repeat", "STR", "Start_Locus", "End_Locus"]

        logging.info(F"The STRs with core length {coreLength} were read from chromosome folder {primary_chromosome_name} ...")
        logging.info(F"The dataset includes {df_str_map.shape[0]} STRs (row) and {df_str_map.shape[1]} columns")
        logging.info(df_str_map.head(10))

        # ==============================================================================================================================================
        logging.info("************** First phase was started (Calculating datasets based on pure STRs **************")

        keys = df_str_map["STR"].str.upper().unique().tolist()
        keys.sort()

        logging.info(F"The list of difference STRs with core length {coreLength} in the chromosome {chromosome_name}:")
        logging.info(keys)
        logging.info("\n")

        progressbar_1.refresh()

        str_categories_total_count = 0
        progressbar_3 = tqdm(keys, desc="processing the STRs", leave=False, bar_format=bar_format, dynamic_ncols=False)
        for key in progressbar_3:
            progressbar_3.set_postfix_str(F"processing the STR {key}")
            filtered_df = df_str_map[df_str_map["STR"].str.upper() == key].copy()
            filtered_df.sort_values(by="Start_Locus", ascending=True, inplace=True)

            dict_value = total_colonies_results[coreLength].setdefault(key, {"dataset": {}, "colonies": {}})
            dict_value["dataset"][chromosome_name] = filtered_df
            str_categories_total_count += dict_value["dataset"][chromosome_name].shape[0]

            # progressbar

            dict_value["colonies"][chromosome_name] = calculate_colonies(dict_value["dataset"][chromosome_name], key, chromosome_name)
            dict_value["dataset"][chromosome_name].to_csv(folder_of_current_colonies_members_result + os.sep +
                                                          F"{chromosome_name}_D3Us_{key}.csv", index=False)

        if str_categories_total_count != df_str_map.shape[0]:
            print(F"###### Error: The separation motif categories has problem! We count total {str_categories_total_count} STRs in categories while "
                  F"there is {df_str_map.shape[0]} STRs ######")
            raise ValueError(F"###### Error: The separation motif categories has problem! We count total {str_categories_total_count} STRs in categories while "
                  F"there is {df_str_map.shape[0]} STRs ######")

        progressbar_1.clear()

        logging.info("*************** First phase was ended (Calculating datasets based on pure STRs ***************")
        # ==============================================================================================================================================
        # ==============================================================================================================================================
        logging.info("************ Second phase was started (Calculating datasets based on grouped STRs ************")

        core_group = set()
        for str_key in keys:
            reversed_str_key = reverse_str(str_key)
            if F"{str_key}-{reversed_str_key}" not in core_group and F"{reversed_str_key}-{str_key}" not in core_group:
                core_group.add(F"{str_key}-{reversed_str_key}")
        key_group_list = list(core_group)

        key_group_list.sort()

        logging.info(F"The list of difference STRs group with core length {coreLength} in the chromosome {chromosome_name}:")
        logging.info(key_group_list)
        logging.info("\n")

        progressbar_1.refresh()

        str_categories_total_count = 0
        progressbar_3 = tqdm(key_group_list, desc="processing the STRs group", leave=False, bar_format=bar_format, dynamic_ncols=False)
        for key in progressbar_3:
            progressbar_3.set_postfix_str(F"processing the STR {key}")
            key_1, sep, key_2 = key.partition("-")
            # dataset_df_1 = total_colonies_results[coreLength][key_1]["dataset"][chromosome_name]
            # dataset_df_2 = total_colonies_results[coreLength][key_2]["dataset"][chromosome_name]
            dataset_df_1 = total_colonies_results[coreLength].get(key_1, {}).get("dataset", {}).get(chromosome_name, pd.DataFrame())
            dataset_df_2 = total_colonies_results[coreLength].get(key_2, {}).get("dataset", {}).get(chromosome_name, pd.DataFrame())
            merged_df = pd.concat([dataset_df_1, dataset_df_2], axis=0, ignore_index=True)
            merged_df.sort_values(by="Start_Locus", ascending=True, inplace=True)
            dict_value = total_colonies_results[coreLength].setdefault(key, {"dataset": {}, "colonies": {}})
            dict_value["dataset"][chromosome_name] = merged_df
            str_categories_total_count += dict_value["dataset"][chromosome_name].shape[0]
            # merged_df.to_csv()

            # progressbar

            dict_value["colonies"][chromosome_name] = calculate_colonies(dict_value["dataset"][chromosome_name], key, chromosome_name)
            dict_value["dataset"][chromosome_name].to_csv(folder_of_current_colonies_members_result + os.sep +
                                                          F"Human_{chromosome_name}_D3Us_{key}_Map_Result.csv", index=False)

        if str_categories_total_count != df_str_map.shape[0]:
            print(F"###### Error: The separation motif categories has problem! We count total {str_categories_total_count} STRs in categories while "
                  F"there is {df_str_map.shape[0]} STRs ######")
            raise ValueError(F"###### Error: The separation motif categories has problem! We count total {str_categories_total_count} STRs in categories while "
                  F"there is {df_str_map.shape[0]} STRs ######")

        progressbar_1.clear()

        logging.info("************* Second phase was ended (Calculating datasets based on grouped STRs *************")
        # ==============================================================================================================================================
        progressbar_1.refresh()

    print("\n\n========================================================================================================")
    print("============================================= final result =============================================")
    print("========================================================================================================\n")

    for core_length_key in total_colonies_results.keys():
        dict_based_core_length = total_colonies_results.get(core_length_key)
        for str_key in dict_based_core_length.keys():
            dict_based_str = dict_based_core_length.get(str_key)
            all_columns = set()

            for chromosome_key in dict_based_str["colonies"].keys():
                all_columns.update(dict_based_str["colonies"].get(chromosome_key).keys())
            table_columns_raw = sorted(list(all_columns))
            table_columns_final = []
            for col_name in table_columns_raw:
                table_columns_final.append(F"C({col_name})")

            table_indexes = sorted(list(dict_based_str["colonies"].keys()))
            result_table = pd.DataFrame(0, index=table_indexes, columns=table_columns_final)

            for chromosome_name in dict_based_str["colonies"].keys():
                chromosome_dict = dict_based_str["colonies"].get(chromosome_name)
                for colony_type_key, member_count in chromosome_dict.items():
                    result_table.at[chromosome_name, F"C({colony_type_key})"] = member_count
            print(F"\n****** Result of colonies count based on {str_key} str ******\n")
            print(result_table)
            print("\n")
            print(104 * "-")
            result_table.to_csv(os.path.join(folder_of_current_colonies_result, F"Human_{str_key}_Colonies_Map_Result.csv"), index=True)

            merged_str_map: pd.DataFrame = pd.concat(list(dict_based_str["dataset"].values()), axis=0, ignore_index=True)
            merged_str_map.to_csv(os.path.join(folder_of_current_colonies_members_result, F"Human_Integrated_D3Us_{str_key}_Map_Result.csv"), index=False)


    # save data structures
    final_results_data_path = os.path.join(folder_of_current_data_structures, "total_colonies_data_results.pkl")
    with open(final_results_data_path, "wb") as f:  # wb = write binary
        pickle.dump(total_colonies_results, f, protocol=pickle.HIGHEST_PROTOCOL)

    return final_results_data_path
