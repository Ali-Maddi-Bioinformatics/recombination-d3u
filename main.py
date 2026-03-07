import sys

import circos_plot
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

    circos_plot.plot(genome_folder_address, primary_patterns_file_address, "")