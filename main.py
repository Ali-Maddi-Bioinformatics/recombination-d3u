import sys

import circos_plot
import statistics_analysis
import colonies

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print(F"Number of argument is {len(sys.argv)}")
        print("Usage: python ColonyDetection.py <The genome folder address> <Input file address>")
        sys.exit(1)

    genome_folder_address = sys.argv[1]
    primary_patterns_file_address = sys.argv[2]

    # str_map_file_address = sys.argv[1]
    # target_str_core = sys.argv[2]
    # target_str_repeat = int(sys.argv[3])

    # final_results_data_path = colonies.process(genome_folder_address, primary_patterns_file_address)
    # circos_plot.plot(genome_folder_address, primary_patterns_file_address, final_results_data_path)

    # circos_plot.plot(genome_folder_address, primary_patterns_file_address, "", "(CT)3-(TC)3")

    BASE_STR_LIST = [
        "(AC)3-(CA)3",
        "(AG)3-(GA)3",
        "(AT)3-(TA)3",
        "(CG)3-(GC)3",
        "(CT)3-(TC)3",
        "(GT)3-(TG)3",
    ]
    circos_plot.plot_batch(genome_folder_address, primary_patterns_file_address, "", BASE_STR_LIST)
    # TARGET_RESULT_SET = "Results_Set_2026_02_16_14_56_49"
    #
    # stats_outputs = statistics_analysis.analyze(
    #     genome_folder_address=genome_folder_address,
    #     base_str_list=BASE_STR_LIST,
    #     result_set_name=TARGET_RESULT_SET,
    #     bin_size=1_000_000,
    #     local_radius_bins=5,
    #     min_colony_members=3,
    #     telomere_radius_bp=1_000_000,
    #     alpha=0.05,
    # )
    #
    # print(stats_outputs)