import numpy as np
import pandas as pd
import os
import re
import math
import time
from subprocess import run, call
from Bio import SeqIO
from Bio.Seq import Seq



#USED
def extract_kmer_with_flanks(kmer,fasta_file, add_length: int = 50):
    kmer = kmer.lower()
    rev_kmer = str(Seq(kmer).reverse_complement())
    k = len(kmer)
    results = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq)
        seq_len = len(seq)
        for target in [kmer, rev_kmer]:
            i = 0
            while True:
                i = seq.find(target, i)
                if i == -1:
                    break
                start = max(0, i - add_length)
                end = min(seq_len, i + k + add_length)
                left = seq[start:i]
                middle = seq[i:i + k].upper()
                right = seq[i + k:end]
                region = left + middle + right
                results.append(region)
                i += 1  # move forward to allow overlapping matches
    return results
#USED
def filter_gene_description(info_from_GFF, desired_columns:list = ["gene", "gene_biotype", "Name"]):
    attrs = dict(field.split('=') for field in info_from_GFF.strip(';').split(';'))
    filtered = {key: attrs[key] for key in desired_columns if key in attrs}
    return filtered
#USED
def align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict):

    print("-------------Align To Genomes----------------- OUTPUT: alignments.csv, no_alignments.csv")

    if not os.path.exists("ref_seq_index.1.bt2"):
        run(f"bowtie2-build {ref_genome_file} ref_seq_index", shell=True)

    fasta_path = "all_queries.fa"
    if sequence_matrix != []:

        matrix = {}
        for i in range(len(sequence_matrix)):
            for j in range(len(sequence_matrix[0])):
                print(sequence_matrix[i][j])
                if sequence_matrix[i][j][0] != "NotPresent":
                    matrix[kmers[j], labelled_genomes[i][2:]] = sequence_matrix[i][j]

        with open(fasta_path, "w") as out:
            count = 0
            for (row, col), seqs in matrix.items():
                for seq in seqs:
                    out.write(f">query_{row}_{col}_{count}\n{seq}\n")
                    count += 1

    call([f"bowtie2 --very-sensitive-local -p 8 -x ref_seq_index -f {fasta_path} -S output.sam"], shell=True)

    alignments = []
    noalignments = []
    sam_path = "output.sam"
    new_lines = []
    print("opening SAM file: ", sam_path)
    with open("output.sam", "r") as sam:
        lines = sam.readlines()
        for line in lines:
            if line[0] == "@":
                continue
            line = line.strip().split("\t")
            new_lines.append(dict({
                "query_name":line[0].strip(),
                "is_reverse":line[1].strip() == "16",
                "reference_name":line[2].strip(),
                "reference_start":line[3].strip(),
                "map_quality":line[4].strip(),
                "cigar":line[5].strip(),
                "rnext":line[6].strip(),
                "pnext":line[7].strip(),
                "template_length":line[8].strip(),
                "query_sequence":line[9].strip(),
            }))

    for read in new_lines:
        # Extract kmer_id and genome_id from query name: query_row_col_count
        match = re.match(r"query_(.+?)_(.+?)_(\d+)", read["query_name"])
        if not match:
            continue
        kmer_id, genome_id, _ = match.groups()
        for labelled_genome in labelled_genomes:
            if genome_id == labelled_genome[2:]:
                genome_id = labelled_genome
            else:
                continue

        strand = "-" if read["is_reverse"] else "+"
        aln_start = read["reference_start"]  # 0-based leftmost position
        ref_name = read["reference_name"]

        if read["reference_name"] == "*":
            noalignments.append(
                {
                    "kmer_id": kmer_id,
                    "coefficient":kmers_coeffs_genomes_dict[kmer_id],
                    "genome_id": genome_id,
                    "ref": ref_name,
                    "pos": aln_start,
                    "seq": read["query_sequence"],
                    "strand":strand
                }
            )
        else:
            genes = find_in_GFF2(aln_start, str(ref_name), GFF_file, desired_columns = ["gene", "Name", "Note"])
            alignments.append({
                "kmer_id": kmer_id,
                "coefficient":kmers_coeffs_genomes_dict[kmer_id],
                "genome_id": genome_id,
                "ref": ref_name,
                "pos": aln_start,
                "seq": read["query_sequence"],
                "strand":strand,
                "genes": genes,
            })
    df = pd.DataFrame(alignments)
    df_n = pd.DataFrame(noalignments)
    
    df.to_csv("alignments.csv")
    df_n.to_csv("no_alignments.csv")





    return 0 ,0 ,0
#USED
def find_in_GFF2(location, chromosome_id ,GFF_file, padding:int =50, kmer_length:int = 13, desired_columns = []):

    with open(GFF_file, "r") as GFF:
        lines  = GFF.readlines()
        genes = ""
        for j in range(len(lines)):
            line = lines[j].strip().split("\t")
            if len(line) < 6:
                continue
            if line[2].strip() == "region" or "#" in str(line[0]) or str(line[0]).strip() != chromosome_id:
                continue
            
            nextLine = ""
            if j+1 < len(lines)-1 and "#" not in lines[j+1][0]:
                nextLine = lines[j+1].strip().split("\t")
            middle_of_align = padding + round(kmer_length/2)
            start = int(line[3])
            end = int(line[4])
            if int(location)+middle_of_align >= start and int(location)+middle_of_align <= end:
                add = f"{line[2].strip()}, {chromosome_id}, {start}, {location}, {end}, {filter_gene_description(line[8], desired_columns)}|"
                genes += add
            elif nextLine != "" and int(location)+middle_of_align >= end and int(location)+middle_of_align <= int(nextLine[3]):
                genes += f"intergenic, {chromosome_id}, {end}, {location}, {nextLine[3]}|"

        if genes == "":
            genes = "NoGenesInRef"
        return genes
#USED
def find_genes_alignment(significant_kmers, species, antibiotic, GFF_file, ref_genome_file = "Not added", reduce_genomes_to:int = 0, reduce_kmers_to:int = 0):
    '''
    Find genes associated with the kmers in the reference genome.
    The output files are genes.csv, intergenic.csv and Kmer_genome_loc_sequence.txt.
    The function returns 1 if the finding was successful.
    Parameters:
        significant_kmers: str - path to the significant kmers file
        species: str - name of the species
        antibiotic: str - name of the antibiotic
        ref_genome_file: str - path to the reference genome file
    Creates files:
        - genes.csv
        - intergenic.csv
        - Kmer_genome_loc_sequence.txt
    Returns:
        int: 1 if the finding was successful.
    '''

    #make file with genomes, kmers, locations in genome, larger DNA piece around kmer

    print("-------- Finding genes in genomes ----------- OUTPUT FILE: genes.csv, intergenic.csv, Kmer_genome_loc_sequence.txt")


    species = species.lower()
    antibiotic = antibiotic.lower()
    genome_folder = species.replace(' ', '_') + "__" + antibiotic


    data_frame = []
    kmers_coeffs_genomes_dict = {}
    kmers_genomes_dict = {}
    #read significant kmers from file with their original genomes and coefficients
    #       determine proportion of resistant and susceptible genomes

    genomes = ""
    print("         Processing significant k-mer file: ", significant_kmers)
    with open(significant_kmers, "r") as file:
        lines = file.readlines()
        for line in lines:
            line = line.strip().split("\t")
            data_frame.append([line[0].strip(), float(line[1].strip()), line[2].strip(), line[3].strip()])
            genomes += " " + line[3].strip()
            
            
        genomes = genomes.strip().split(" ")
        genomes = sorted(set(genomes), key=lambda x: int(x.split(".")[1]))
    data_frame = sorted(data_frame, key=lambda x: x[1], reverse=True)          #sort data frame by kmer
    print(f"            K-mers sorted for coefficients:  {len(data_frame)}")
    for line in data_frame:
        kmers_coeffs_genomes_dict[line[0]] = [line[1], line[2]]
        kmers_genomes_dict[line[0]] = line[3].strip().split(" ")
    
    labelled_genomes = []
    genome_folder_genome_names = os.listdir(f"Genomes/" + genome_folder)
    for genome in genomes:
        for genome_file in genome_folder_genome_names:
            if genome in genome_file:
                labelled_genomes.append(genome_file[:-4])
                break
    labelled_genomes = sorted(labelled_genomes, key = lambda x: x[0])
    print("             Unique genomes across all kmers: ", len(labelled_genomes))
    print("             Resistant genomes: ", len([elem for elem in labelled_genomes if elem[0] == "1"]))
    #K-mers with higher coefficients are first
    
    #print(labelled_genomes)
    # print("         Writing k-mers and source genomes to file: kmers_genomes.txt")
    # with open("kmers_genomes.txt", "w") as file:
    #     for kmer in kmers_genomes_dict:
    #         file.write(f"{kmer}\t{kmers_genomes_dict[kmer]}\n")
    
    kmers = list(kmers_genomes_dict.keys())
    kmer_count_res_sus_dict = {}
    print("         Making kmer, genome count, resistant, susceptible dictionary")
    for i in range(len(kmers)):
        genome_count = 0
        resistant_count = 0
        susceptible_count = 0
        for j in range(len(labelled_genomes)):
            if labelled_genomes[j][2:] in kmers_genomes_dict[kmers[i]]:
                genome_count += 1
                if int(labelled_genomes[j][0]) == 1:
                    resistant_count += 1
                else:
                    susceptible_count += 1
        kmer_count_res_sus_dict[kmers[i]] = f"T:{genome_count}; R:{resistant_count}/S:{susceptible_count}({round(resistant_count / genome_count*100, 2)}%)"
    longer_sequences = []
    sequence_matrix = []

    if reduce_genomes_to != 0 and reduce_genomes_to < len(labelled_genomes):
        labelled_genomes = labelled_genomes[:math.floor(reduce_genomes_to/2)] + labelled_genomes[-math.ceil(reduce_genomes_to/2):]
    if reduce_kmers_to != 0 and reduce_kmers_to < len(kmers_genomes_dict.keys()):
        kmers = list(kmers_genomes_dict.keys())[:reduce_kmers_to]

    if os.path.exists("all_queries.fa"):
        align_to_genome2([], ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict)
        return 0
    print(f"         Expanding {len(kmers)} k-mers from their original {len(labelled_genomes)} genomes. {len(labelled_genomes) * len(kmers)}")
    #make sequence_matrix
    kmer_count = 0
    times = []
    kmers_left = len(kmers)
    for i in range(len(kmers)):
        columns = []
        start = time.time()
        for j in range(len(labelled_genomes)):
            if labelled_genomes[j][2:] in kmers_genomes_dict[kmers[i]]:
                genome_file = "Genomes/" + genome_folder + "/" + labelled_genomes[j] + ".fna"
                #print("Genome file: ", genome_file)
                longer_sequences = extract_kmer_with_flanks(kmers[i], genome_file)
                columns.append(longer_sequences)
            else:
                columns.append(["NotPresent"])
        kmer_count += 1
        kmers_left -= 1
        duration = time.time() - start
        times.append(duration)
        print(f"\r         K-mers expanded:  {kmer_count} / {len(kmers)}\tTime per kmer:  {round(np.median(np.array(times)), 2)} s\tEstimated time: {round(kmers_left*np.median(times)/60,2)} min", end="", flush=True)
        sequence_matrix.append(columns)
    
    print("\n")    
    
    sequence_matrix = np.transpose(np.array(sequence_matrix))
    sequences = ""

    ref_genome_row = []
    ref_genome_row_string = ""

    #open resulting table file for writing
    with open("kmers_genomes_sequences_table.csv", "w") as file:
        file.write("Table")
        for i in range(len(kmers)):
            #add first row with kmers and their parameters in genomes.
            file.write(f", {kmers[i]} Counts(g; R; S): {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {round(kmers_coeffs_genomes_dict[kmers[i]][0],7)}")
            #make row for ref genome separately
            ref_genome_row.append(extract_kmer_with_flanks(kmers[i], ref_genome_file))
        
        #format genome rows and write to file
        shape = sequence_matrix.shape
        for i in range(shape[0]):
            file.write("\n" + labelled_genomes[i])
            for j in range(shape[1]):
            
                sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "")
                if len(sequence_matrix[i][j]) > 1:
                    sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
                file.write("," + sequences)
    
    align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict)          #align the sequences to the reference genome
