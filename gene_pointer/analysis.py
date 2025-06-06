import numpy as np
import pandas as pd
import os
import re
import math
import time
from subprocess import run, call
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path

#USED MINREQ
def find_significant_kmers(data_pheno_path: str, classifier: str = "log", kmer_length: int = 13, POPULATION_CORRECTION: bool = True):
    '''
        Find kmers associated with the resistance phenotype in input genomed using phenotypeseeker.
        The kmer length is set to 13 by default. The classifier can be either "log" or "RF".
        The output file is filtered_kmers_and_coeffs.txt.
        The function returns 1 if the kmers were found, 0 if not.

        Also creates a ML model for predicting the phenotype in unseen genomes.

        Parameters:
            data_pheno_path: str - path to the data.pheno file
            classifier: str - classifier to use, either "log" or "RF"
            kmer_length: int - length of the kmers to find, default is 13
        Creates files:
            - k-mers_and_coefficients_in_<classifier>_model_<antibiotic>.txt
            - and more which are note used in the following pipeline
        Returns:
            int: 1 if the kmers were found, 0 if not
    '''

    print("---------- Looking for significant kmers ---------- OUTPUT FILE: k-mers_and_coefficients_in_<classifier>_model_<antibiotic>.txt")
    antibiotic = ""
    try:
        with open(data_pheno_path, "r") as dataphenofile:           #open data.pheno file
            lines = dataphenofile.readlines()
            if len(lines) < 2:
                print("    No data in file")
                return 0
            #print("DataPheno Header: ", lines[0].strip().split("\t"))
            antibiotic = lines[0].strip().split()[2]           #get antibiotic name from first line of file
    except FileNotFoundError:
        print("    Data file not found")
        return 0
    
    # print("    Looking for file ./k-mers_and_coefficients_in_" + (classifier if classifier != "log" else "log_reg") + "_model" + "_" + antibiotic.capitalize() + ".txt")
    # if not os.path.exists("./k-mers_and_coefficients_in_" + (classifier if classifier != "log" else "log_reg") + "_model" + "_" + antibiotic.capitalize() + ".txt"):
    #     print("        File not found")

    absolute_path = str(Path(data_pheno_path).resolve())
    if POPULATION_CORRECTION:
        popcorr = "-w"
    else:
        popcorr = ""
    print(">>>>>>>>>>>>>>>>> RUNNING: phenotypeseeker modeling " + absolute_path + " -w -bc " + classifier + " -l " + str(kmer_length) + " <<<<<<<<<<<<<<<<<")             #if the k-mers and coefficients file does not exist, call phenotypeseeker to find significant kmers
    call(["phenotypeseeker modeling " + absolute_path + " " + popcorr + " -bc " + classifier + " -l " + str(kmer_length)], shell=True)
    print("    Finished executing phenotypeseeker")
    try:
        #add dictionary for classifier names if needed
        if classifier == "log":
            classifier = "log_reg"
        #check if the k-mers and coefficients file exists. redundant if it is checked in previous function
        with open("k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt", "r") as kmers_coeffs:      #open k-mers and coefficients file with significant kmers found by phenotypeseeker
            lines = kmers_coeffs.readlines()
            if len(lines) == 0:
                print("    No significant kmers found")
                return 0

            kmers = []
            coeffs = []
            whole_lines = []
            for line in lines[1:]:                  
                kmers.append(line.split("\t")[0])
                coeffs.append(line.split("\t")[1])
                whole_lines.append(line.strip().split("\t"))
            print("    Number of significant kmers from PhenotypeSeeker: ", len(whole_lines))

            whole_lines = sorted(whole_lines, key=lambda x: x[1], reverse=True)

            with open("significant_lines_from_kmer_coeff_file.txt", "w") as file:
                for line in whole_lines:
                    line = "\t".join(line).replace("|", "")
                    file.write(f"{line}\n")
     
    except FileNotFoundError:
        print("    k-mers and coefficients file not found: No kmers found by phenotypeseeker")
        return 0
    return 1
#USED MINREQ
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
#USED MINREQ
def filter_gene_description(info_from_GFF, desired_columns:list = ["gene", "gene_biotype", "Name"]):
    attrs = dict(field.split('=') for field in info_from_GFF.strip(';').split(';'))
    filtered = {key: attrs[key] for key in desired_columns if key in attrs}
    return filtered
#USED MINREQ
def align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict):

    refseq_ID = ref_genome_file.strip().split("/")[-1][:-4].replace(".","_")
    fasta_path = f"all_queries_for_{refseq_ID}.fa"
    print(f"---------- Align To Genome {refseq_ID} ---------- OUTPUT FILES: Aligned_kmer_results.csv, Un_aligned_kmer_results.csv")

    print(f"    Reference genome id: {refseq_ID}")
    refseq_index_file = f"refseq_{refseq_ID}_index"

    if not os.path.exists(f"{refseq_index_file}.1.bt2"):
        run(f"bowtie2-build {ref_genome_file} {refseq_index_file}", shell=True)

    
    if sequence_matrix != []:

        matrix = {}
        for i in range(len(sequence_matrix)):
            for j in range(len(sequence_matrix[0])):
                #print(sequence_matrix[i][j])
                if sequence_matrix[i][j][0] != "NotPresent":
                    matrix[kmers[j], labelled_genomes[i][2:]] = sequence_matrix[i][j]

        with open(fasta_path, "w") as out:
            count = 0
            for (row, col), seqs in matrix.items():
                for seq in seqs:
                    out.write(f">query_{row}_{col}_{count}\n{seq}\n")
                    count += 1

    sam_path = f"aligned_to_{refseq_ID}.sam"

    call([f"bowtie2 --very-sensitive-local -p 8 -x {refseq_index_file} -f {fasta_path} -S {sam_path}"], shell=True)

    alignments = []
    noalignments = []
    
    new_lines = []
    with open(sam_path, "r") as sam:
        print(f"---------- Extracting alignment results from SAM file: {sam_path} ----------")
        lines = sam.readlines()
        for line in lines:
            print(f"\r    Lines read: {len(new_lines)} / {len(lines)}", end="\r", flush=True)
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
    print("\n    Lines processed: ", len(new_lines), " / ", len(lines), " (metadata lines not proccessed)")

    
    times = [0]
    for read in new_lines:
        start = time.time()
        # Extract kmer_id and genome_id from query name: query_row_col_count
        print(f"\r    SAM lines processed: {len(alignments) + len(noalignments)} / {len(new_lines)}   ET:{round((len(new_lines)-len(alignments) + len(noalignments))*np.median(times)/60,2)} min", end="\r", flush=True)
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
                    "K-mer": kmer_id,
                    "Coefficient in ML model":kmers_coeffs_genomes_dict[kmer_id][0],
                    "Genome_ID": genome_id,
                    "AlignementBaseSeq_ID": ref_name,
                    "AlignmentLeftSide_Pos": aln_start,
                    "AlignementSeq": read["query_sequence"],
                    "Strand":strand
                }
            )
        else:
            genes = find_in_GFF2(aln_start, str(ref_name), GFF_file, desired_columns = ["gene", "Name", "Note"])
            alignments.append({
                #includes a mix of reverse and forward kmers
                "K-mer": kmer_id,
                "Coefficient in ML model":kmers_coeffs_genomes_dict[kmer_id][0],
                "Genome_ID": genome_id,
                "AlignementBaseSeq_ID": ref_name,
                "AlignmentLeftSide_Pos": aln_start,
                "AlignementSeq": read["query_sequence"],
                "Strand":strand,
                "Genes from base sequence": genes,
            })
        times.append(time.time() - start)

    print("    Making alignments.csv and no_alignments.csv files....")
    df = pd.DataFrame(alignments)
    df_n = pd.DataFrame(noalignments)
    
    df.to_csv("Aligned_kmer_results.csv")
    df_n.to_csv("Un_aligned_kmer_results.csv")
#USED MINREQ
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
#USED MINREQ
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

    print("---------- Finding genes in genomes ---------- OUTPUT FILE: kmers_genomes_sequences_table.csv, alignments.csv, no_alignments.csv")


    species = species.lower()
    antibiotic = antibiotic.lower()
    genome_folder = species.replace(' ', '_') + "__" + antibiotic


    data_frame = []
    kmers_coeffs_genomes_dict = {}
    kmers_genomes_dict = {}
    #read significant kmers from file with their original genomes and coefficients
    #       determine proportion of resistant and susceptible genomes

    genomes = ""
    print("    Processing significant k-mer file: ", significant_kmers)
    with open(significant_kmers, "r") as file:
        lines = file.readlines()
        for line in lines[1:]:
            line = line.replace("|", "")
            line = line.strip().split("\t")
            data_frame.append([line[0].strip(), float(line[1].strip()), line[2].strip(), line[3].strip()])
            genomes += " " + line[3].strip()
            
            
        genomes = genomes.strip().split(" ")
        genomes = list(set(genomes))  #remove duplicates
        genomes = sorted(genomes, key=lambda x: int(x.split(".")[1]))
    data_frame = sorted(data_frame, key=lambda x: x[1], reverse=True)          #sort data frame by kmer
    print(f"        K-mers sorted for coefficients:  {len(data_frame)}")
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
    print("        Unique genomes across all kmers: ", len(labelled_genomes))
    print("        Resistant genomes: ", len([elem for elem in labelled_genomes if elem[0] == "1"]))
    #K-mers with higher coefficients are first
    
    #print(labelled_genomes)
    # print("         Writing k-mers and source genomes to file: kmers_genomes.txt")
    # with open("kmers_genomes.txt", "w") as file:
    #     for kmer in kmers_genomes_dict:
    #         file.write(f"{kmer}\t{kmers_genomes_dict[kmer]}\n")
    
    kmers = list(kmers_genomes_dict.keys())
    kmer_count_res_sus_dict = {}
    print("    Making kmer, genome count, resistant, susceptible dictionary")
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

    #DO THIS ON MORE THAN 1 REFERENCE/ANNOTATED GENOME
    refseq_ID = ref_genome_file.strip().split("/")[-1][:-4].replace(".","_")
    fasta_path = f"all_queries_for_{refseq_ID}.fa"
    if os.path.exists(fasta_path) and os.path.exists("kmers_genomes_sequences_table.csv"):
        align_to_genome2([], ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict)
        return 0
    
    
    print(f"    Expanding {len(kmers)} k-mers from their original {len(labelled_genomes)} genomes. {len(labelled_genomes) * len(kmers)}")
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
        print(f"\r        K-mers expanded:  {kmer_count} / {len(kmers)}\tTime per kmer:  {round(np.median(np.array(times)), 2)} s\tEstimated time: {round(kmers_left*np.median(times)/60,2)} min", end="", flush=True)
        sequence_matrix.append(columns)
    
    print("\n")    
    
    sequence_matrix = np.transpose(np.array(sequence_matrix))
    sequences = ""

    # ref_genome_row = []
    # ref_genome_row_string = ""
    
    #open resulting table file for writing
    with open("kmers_genomes_sequences_table.csv", "w") as file:
        file.write("Table")
        for i in range(len(kmers)):
            #add first row with kmers and their parameters in genomes.
            file.write(f", {kmers[i]} Counts(g; R; S): {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {round(kmers_coeffs_genomes_dict[kmers[i]][0],7)}")
            
            #make row for ref genome separately
            #ref_genome_row.append(extract_kmer_with_flanks(kmers[i], ref_genome_file)
        
        #format genome rows and write to file
        shape = sequence_matrix.shape
        for i in range(shape[0]):
            file.write("\n" + labelled_genomes[i])
            for j in range(shape[1]):
            
                sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "")
                if len(sequence_matrix[i][j]) > 1:
                    sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
                file.write("," + sequences)
    
    #DO THIS FOR MORE THAN 1 REFERENCE/ANNOTATED GENOME
    align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict)          #align the sequences to the reference genome


#WIP
def extract_kmer_with_flanks2(kmer, fasta_file, add_length: int = 50):
    rev_kmer = str(Seq(kmer).reverse_complement())
    k = len(kmer)
    results = []

#WIP
def combine_GFFs(GFF_files_folder, output_file):
    '''
    Combine multiple GFF files into one.
    Parameters:
        GFF_files: list of str - paths to the GFF files to combine
        output_file: str - path to the output file
    Creates:
        - output_file: combined GFF file
    '''
    with open(output_file, "w") as out_file:
        for gff_file in GFF_files_folder.glob("*.gff"):
            with open(gff_file, "r") as in_file:
                lines = in_file.readlines()
                for line in lines:
                    if not line.startswith("#"):
                        out_file.write(line)
    print(f"Combined GFF files into {output_file}")

#MAYBE USED, would represent better the common elements between all input genomes that are associated with the phenotype?
def combine_kmers_by_locations(significant_kmers):
    #make file with combined kmers
    #only combine kmers that are close enough to eachother
    #fill in gaps with "_"
    pass