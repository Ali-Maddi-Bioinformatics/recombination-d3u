import os
import re

from tqdm import tqdm
from datetime import datetime

bar_format = (
        "{desc:45}: {percentage:3.0f}%|"
        "{bar:30}| "
        "{n}/{total} "
        "[{elapsed}<{remaining}, {rate_fmt}, {postfix}]"
    )

def read_chromosomes(genome_folder_address):
    primary_chromosome_dict = {}
    chromosomes_files_list = []
    for file in os.listdir(genome_folder_address):
        if os.path.isfile(os.path.join(genome_folder_address, file)):
            chromosomes_files_list.append(file)
    progressbar_0 = tqdm(chromosomes_files_list, desc="Reading the chromosomes files", bar_format=bar_format, dynamic_ncols=False)
    for chromosome_file in progressbar_0:
        dna_parts = []
        chromosome_name = ""
        with open(os.path.join(genome_folder_address, chromosome_file), 'r') as chromosome_file_content:
            # for line in tqdm(chromosome_file_content, desc="Reading content of chromosome line by line", bar_format=bar_format, dynamic_ncols=False):
            for line in chromosome_file_content:
                line = line.strip()
                if not line:
                    continue
                elif len(line) == 1 and line == '>':
                    raise ValueError("Invalid Chromosome Name...")
                elif line[0] == '>':
                    if chromosome_name and dna_parts:
                        primary_chromosome_dict[chromosome_name] = "".join(dna_parts)
                        chromosome_name = line[1:] + F"_{int(datetime.now().timestamp() * 1000)}"
                        dna_parts = []
                    else:
                        chromosome_name = line[1:] + F"_{int(datetime.now().timestamp() * 1000)}"
                else:
                    dna_parts.append(line)
            if chromosome_name and dna_parts:
                primary_chromosome_dict[chromosome_name] = "".join(dna_parts)
                chromosome_name = ""
                dna_parts = []

    return primary_chromosome_dict


def unpack_str(packed_str):
    packed_str = packed_str.upper()
    pattern = r"\((.*?)\)(\d+)"
    match = re.match(pattern, packed_str)

    if match:
        core = match.group(1)
        repeat = int(match.group(2))
        return {"core": core, "repeat": repeat, "extend": F"{core * repeat}", "pack": packed_str}
    raise ValueError("Invalid STR string...")


def read_patterns(primary_patterns_file_address):
    patterns_list = []
    with open(primary_patterns_file_address, 'r') as content:
        for line in content:
            line = line.strip()
            patterns_list.append(unpack_str(line))
    return patterns_list


def standard_name(chromosome_name):
    pattern = r"(.*?)(\d+)"
    match = re.match(pattern, chromosome_name)

    if match:
        number = len(match.group(2))
        if number < 2:
            return match.group(1) + "0" + match.group(2)
    return chromosome_name