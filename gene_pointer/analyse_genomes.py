from subprocess import call, run
import numpy as np
import pandas as pd
import os
import time
import random
import math
from Bio import SeqIO
from Bio.Seq import Seq
#import pysam
import re

def find_significant_kmers(data_pheno_path: str, classifier: str = "log", kmer_length: int = 13):
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
            - filtered_kmers_and_coeffs.txt
            - and more which are note used in the following pipeline
        Returns:
            int: 1 if the kmers were found, 0 if not
    '''

    print("----- Looking for significant kmers -------------- OUTPUT FILE: filtered_kmers_and_coeffs.txt")
    antibiotic = ""
    try:
        with open(data_pheno_path, "r") as dataphenofile:           #open data.pheno file
            lines = dataphenofile.readlines()
            if len(lines) < 2:
                print("        No data in file")
                return 0
            #print("DataPheno Header: ", lines[0].strip().split("\t"))
            antibiotic = lines[0].strip().split()[2]           #get antibiotic name from first line of file
    except FileNotFoundError:
        print("        Data file not found")
        return 0
    print("         Looking for file ./k-mers_and_coefficients_in_" + (classifier if classifier != "log" else "log_reg") + "_model" + "_" + antibiotic.capitalize() + ".txt")
    if not os.path.exists("./k-mers_and_coefficients_in_" + (classifier if classifier != "log" else "log_reg") + "_model" + "_" + antibiotic.capitalize() + ".txt"):
        print("         EXECUTING: phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier + " -l " + str(kmer_length))             #if the k-mers and coefficients file does not exist, call phenotypeseeker to find significant kmers
    
    # try:
    #     call(["mkdir " + antibiotic], shell=True)           #create a directory for the antibiotic
    # except FileExistsError:
    #     print("         Directory already exists")
    # os.chdir(antibiotic)           #change directory to the antibiotic folder
    # print("         Changing directory to ", antibiotic)
        call(["phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier + " -l " + str(kmer_length)], shell=True)
    # os.chdir("..")           #change directory back to the original folder
        print("         Finished executing phenotypeseeker")
    # else:
    #     print("         Found")
    
    try:
        #add dictionary for classifier names if needed
        if classifier == "log":
            classifier = "log_reg"
        #check if the k-mers and coefficients file exists. redundant if it is checked in previous function
        with open("k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt", "r") as kmers_coeffs:      #open k-mers and coefficients file with significant kmers found by phenotypeseeker
            lines = kmers_coeffs.readlines()
            if len(lines) == 0:
                print("        No significant kmers found")
                return 0

            kmers = []
            coeffs = []
            whole_lines = []
            for line in lines[1:]:                  #take only kmers with non-zero coefficients from the file maybe separate kmer present in all phenotypes and in just resistant phenotypes
                if float(line.split("\t")[1]) != 0:
                    kmers.append(line.split("\t")[0])
                    coeffs.append(line.split("\t")[1])
                    whole_lines.append(line.strip().split("\t"))
            print("         Number of significant kmers from PhenotypeSeeker: ", len(whole_lines))

            whole_lines = sorted(whole_lines, key=lambda x: x[1], reverse=True)

            with open("significant_lines_from_kmer_coeff_file.txt", "w") as file:
                for line in whole_lines:
                    line = "\t".join(line).replace("|", "")
                    file.write(f"{line}\n")
     
    except FileNotFoundError:
        print("        k-mers and coefficients file not found: No kmers found by phenotypeseeker")
        return 0
    return 1

def index_genome(genome_path, output_file = "ref_gnome_kmer_locs.txt", kmer_length = 13):
    '''
        Index the reference genome using glistmaker and glistquery.
        The output file is ref_gnome_kmer_locs.txt.
        The function returns 0 if the index file already exists, 1 if the indexing was successful.
        If the reference genome file does not exist, the function returns 0.

        Parameters:
            reference_genome_path: str - path to the reference genome file
            output_file: str - path to the output file
            kmer_length: int - length of the kmers to find, default is 13
        Creates files:
            - out_<kmer_length>.index
            - ref_gnome_kmer_locs.txt
        Returns:
            int: 0 if the index file already exists, 1 if the indexing was successful
    '''
    print("----- Indexing genome ----------- OUTPUT FILE: ", output_file)
    if os.path.exists(output_file):
        print("         Index file already exists. Stopping indexing")
        return 0
    if not os.path.exists(genome_path):
        print("         Genome file not found")
        return 0
    if not os.path.exists(f"out_{kmer_length}.index"):
        print("         Index file not found")
        print("         Indexing genome")
        call([f"glistmaker {genome_path} -w {kmer_length} --index"], shell=True)           #index reference genome, locations of each kmer in the ref genome
    
    call([f"glistquery out_{kmer_length}.index --locations >> {output_file}"], shell=True)    #write the index file into a text file for reading later
    print("----- Finished indexing genome -----")

def filter_GFF_file(GFF_file, output_file="filtered_gff.txt"):
    
    '''
        Filter the GFF file to only include the elements of interest: element, start position, end position, description.
        The output file is filtered_gff.txt.
        The function returns 0 if the output file already exists, 1 if the filtering was successful.
        If the GFF file does not exist, the function returns 0.
        Parameters:
            GFF_file: str - path to the GFF file
            output_file: str - path to the output file
        Creates files:
            - filtered_gff.txt
        Returns:
            int: 0 if the output file already exists, 1 if the filtering was successful.
    '''

    # if os.path.exists(output_file):
    #     print("        Output file already exists. Stopping filtering")
    #     return 0
    print("----- Reducing GFF file ------------------ OUTPUT FILE: ", output_file)

    try:
        with open(GFF_file, "r") as file:
            lines = file.readlines()
            lines = lines[7:]           #skip the first 7 lines of the GFF file, meta data
            elements = []
            start_end = []
            ids = []
            #columns of GFF file
            # 1. seqid, 2. source, 3. type, 4. start, 5. end, 6. score, 7. strand, 8. phase, 9. description
            previous_start = 0
            id_line = ""
            beginning = 0
            #print(lines[4].strip().split("\t"))
            for i in range(len(lines)):              #taking only the lines with the elements of interest: element type, start, end, id
                if "##sequence" in lines[i].strip() or "###" in lines[i]:
                    continue
                if "##species" in lines[i] and i < 15:
                    continue
                if "##species" in lines[i].strip() and i > 10:
                    beginning = beginning + int(lines[i-2].strip().split("\t")[4])
                    continue
                if "region" == lines[i].strip().split("\t")[2].strip():
                    continue
                line = lines[i].strip().split("\t")
                start = lines[i].strip().split("\t")[3]
                if start != previous_start:
                    elements.append(line[2])
                    start_end.append((beginning + int(line[3]),beginning + int(line[4])))
                    id_line += line[8].strip() + ";"
                    if lines[i+1].strip().split("\t")[2] == "CDS":
                        id_line += lines[i+1].strip().split("\t")[8].strip()
                    ids.append(id_line)
                    id_line = ""
                    previous_start = start
                
        print("         Elements: ", len(elements))
    except FileNotFoundError:
        print("         GFF file not found")
        return 0
    try:
        with open("filtered_gff.txt", "w") as file:          #write the filtered GFF file to a new file
            for i in range(len(elements)):
                file.write(f"{elements[i]}\t{start_end[i][0]}\t{start_end[i][1]}\t{ids[i]}\n")
        print("----- Finished reducing GFF file -----")
    except FileNotFoundError:
        print("         Output file not found")
        return 0

def find_ref_genome_kmers(indexed_kmer_file: str,kmers_coeffs_file: str, max_mismatches: int = 0, output_file: str="pheno_kmers_ref_genome.txt"):
    '''
        Find kmers in the reference genome using the indexed kmer file and the kmers and coefficients file.
        The output file is pheno_kmers_ref_genome.txt.
        Parameters:
            indexed_kmer_file: str - path to the indexed kmer file
            kmers_coeffs_file: str - path to the kmers and coefficients file
            output_file: str - path to the output file default: pheno_kmers_ref_genome.txt
        Creates files:
            - pheno_kmers_ref_genome.txt
        Returns:
            int: 0 if the output file already exists or no kmers were found, 1 if the filtering was successful.
    '''

    kmers_coeffs = []

    try: 
        with open(kmers_coeffs_file, "r") as file:
            print("         Reading kmers and coefficients file: ", kmers_coeffs_file)
            lines = file.readlines()
            for line in lines:
                line = line.strip().split("\t")
                kmers_coeffs.append((line[0], line[1]))
    except FileNotFoundError:
        print(f"        {kmers_coeffs_file} not found")
        return 0

    print("----- Searching for kmers in ref genome---------- OUTPUT FILE: ", output_file)        #Zipping kmers and coefficients together
    if len(kmers_coeffs) == []:
        print("        No kmers found")
        return 0   

    print("         Sorting kmers by coefficients")     # Sort kmers by coefficients, greatest to lowest
    kmers_coeffs = sorted(kmers_coeffs, key=lambda x: x[1])
    print("         Making reverse complement kmer list")  
    kmers_coeffs_revcomp = []
    for line in kmers_coeffs:
        kmer  = line[0]
        kmers_coeffs_revcomp.append((reverse_complement(kmer)))
    #make dictionary instead
    print("         Finding kmers in indexed file. May take a while...")       # Find kmers in indexed reference genome file
    shared_kmers = {}
    try:
        with open(indexed_kmer_file, "r") as file:          #open indexed ref genome
            print("         Reading indexed kmer file: ", indexed_kmer_file)
            lines = file.readlines()


            kmer_count = 0
            times = []
            lines_left = len(lines)
            for i in range(len(lines)):
                start = time.time()
                lines_left -= 1
                if i % 1000 == 0:
                    print(f"\r       Searching lines: {i} / {len(lines)} Time per line: {round(np.median(times), 4)} s Estimated time: {round(np.median(times) * lines_left / 60, 1)} min Kmers found: {kmer_count}", end="", flush=True)
                if lines[i][0] == 0:
                    continue
                #print(lines[i].split("\t"))
                line_start = str((lines[i].split("\t"))[0])
                index_count = int(lines[i].strip().split("\t")[1])
                
                for j in range(len(kmers_coeffs)):
                    kmer = str(kmers_coeffs[j][0])
                    kmer_revcomp = str(kmers_coeffs_revcomp[j])
                    #make a dictionary for the kmers and their locations             
                    if (kmer == line_start) or (kmer_revcomp == line_start):                    #check if a phenotype kmer (or its reverse complement) is in this line of indexed ref genome 
                        
                        index_lines = lines[i+1:i+index_count+1]         #add all locations of the kmer to the shared kmers list
                        for index in index_lines:
                            index = index.strip().split("\t")
                            if int(index[-1]) > max_mismatches:
                                continue
                            if not shared_kmers.get(kmer):
                                shared_kmers[kmer] = [kmers_coeffs[j][1]]
                                kmer_count += 1
                            shared_kmers[kmer].append(index[2])
                        break
                times.append(round(time.time()-start, 5))

    except FileNotFoundError: 
        print("        Indexed kmer file not found")
        return 0

    print(f"     \nNumber of kmers found in indexed file (max_mismatch: {max_mismatches}): ", len(list(shared_kmers.keys())))          # Write phenotype associated kmers found in reference genome to file
    with open(output_file, "w") as file:
        for kmer in shared_kmers.keys():
            locations = ""
            coeff = shared_kmers[kmer][0]
            for location in shared_kmers[kmer][1:]:
                locations += location + "\t"
            file.write(f"{coeff}\t{kmer}\t{locations}\n")
    print("----- Finished searching for kmers -----")
    return 1

def filter_id_line(line):
    '''
    Filter the ID line to only include the elements of interest: Parent, Name, gene, gene_biotype, product.
    The function returns the filtered line.
    Parameters:
        line: str - line to be filtered
    Returns:
        str: filtered line
    '''
    items = line.split(";")
    new_items = []
    for item in items:
        if "Parent=" in item:
            new_items.append(item)
        if "Name=" in item:
            new_items.append(item)
        if "gene" in item:
            new_items.append(item)
        if "gene_biotype=" in item:
            new_items.append(item)
        if "product=" in item:
            new_items.append(item)
        
    return "(" +" ".join(new_items) + ")".replace(",", " ")

def find_genes(pheno_kmers_in_ref_genome_file, filtered_GFF_file):
    '''
    Find genes associated with the kmers in the reference genome.
    The output files are genes.txt, intergenic.txt and unidentified.txt.

    Parameters:
        pheno_kmers_in_ref_genome_file: str - path to the pheno kmers in reference genome file
        filtered_GFF_file: str - path to the filtered GFF file
        min_mismatches: int - minimum number of mismatches to consider a kmer as significant
        output_file: str - path to the output file
    Creates files:
        - genes.txt
        - intergenic.txt
        - unidentified_kmers.txt
    Returns:
        int: 1 if the filtering was successful.

    '''
    print("----- Finding genes ------------------ OUTPUT FILES: genes.csv, intergenic.csv, ")

    kmers_with_locations = []
    coefficients = []
    kmers_without_locations = []
    print("         Extracting kmers and their locations")
    try:

        #open file with kmers associated with phenotype and present in reference genome
        with open(pheno_kmers_in_ref_genome_file, "r") as significant_kmer_file:
            lines = significant_kmer_file.readlines()
            for line in lines:
                line = line.strip().split("\t")
                kmer = line[1]
                locations = line[2:]
                if locations == []:
                    kmers_without_locations.append([kmer])
                else:
                    kmers_with_locations.append([kmer, locations])
                    coefficients.append(line[0])   

    except FileNotFoundError:
        print("         Significant kmers file not found")
        return 0
    
    
    print("         Number of kmers with locations in reference genome (direct map method, alignment based matching will be implemented): ", len(kmers_with_locations))
    print("         Identifying kmer parent elements in reference genome")


    genes, intergenic = find_in_GFF_no_align(kmers_with_locations, filtered_GFF_file)

    print("         UNIDENTIFIED kmers: ", len(kmers_without_locations))
    print("         ELEMENTS found: ", len(genes))
    print("         INTERGENIC found: ", len(intergenic))

    print("         Writing tables..")
    with open("genes.csv", "w") as gene_file:
        genes = sorted(genes, key=lambda x: x[3])
        for line in genes:
            line[5] = filter_id_line(line[5])
            gene_file.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]}\n")
    with open("intergenic.csv", "w") as intergenic_file:
        i = 1
        for line in intergenic:
            intergenic_file.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]}\n")
            i += 1

    if len(kmers_without_locations) == 0:
        print("         No unidentified k-mers")
    else:
        with open("unidentified_kmers.txt", "w") as unidentified_kmers_file:
            
            for line in kmers_without_locations:
                locations = ""
                for item in line[1:]:
                    locations += str(item) + " "
                unidentified_kmers_file.write(f"{line[0]},{locations}\n")
    print("----- Finished finding genes -----")
    return 1

def find_in_GFF_no_align(dnas_with_locations, filtered_GFF_file):
    
    genes =  []
    intergenic = []
    try:   
        with open(filtered_GFF_file, "r") as file:          #open filtered GFF file
            lines = file.readlines()
            for k in range(len(dnas_with_locations)):
                dnaloc = dnas_with_locations[k]
                for j in range(len(lines)):
                    line = lines[j].strip().split("\t")
                    nextLine = ""
                    if j < len(lines)-1:
                        nextLine = lines[j+1].strip().split("\t")
                    if line[0] == "region":
                        continue
                    else:
                        for i in range(len(dnaloc[1])):
                            if int(dnaloc[1][i]) >= int(line[1]) and int(dnaloc[1][i]) <= int(line[2]):
                                add = [dnaloc[0], line[0], line[1], dnaloc[1][i], line[2], line[3]]
                                genes.append(add)
                            elif nextLine != "" and int(dnaloc[1][i]) >= int(line[2]) and int(dnaloc[1][i]) <= int(nextLine[1]):
                                intergenic.append([dnaloc[0], "intergenic",int(line[2]), dnaloc[1][i], nextLine[1]])
        return genes, intergenic
    except FileNotFoundError:
        print("         Filtered GFF file not found")
        return 0,0

def find_in_GFF(dnas_with_locations, filtered_GFF_file):
    
    gene_column =  []
    intergenic_column = []
    try:   
        with open(filtered_GFF_file, "r") as file:          #open filtered GFF file
            lines = file.readlines()
            for k in range(len(dnas_with_locations)):
                dnaloc = dnas_with_locations[k]
                genes = ""
                intergenic = ""
                for j in range(len(lines)):
                    line = lines[j].strip().split("\t")
                    nextLine = ""
                    if j < len(lines)-1:
                        nextLine = lines[j+1].strip().split("\t")
                    if line[0] == "region":
                        continue
                    else:
                        for alignment in dnaloc:
                            middle_of_align = round(len(alignment[1]) / 2)
                            if int(alignment[0]) == 0:
                                continue
                            elif int(alignment[0])+middle_of_align >= int(line[1]) and int(alignment[0])+middle_of_align <= int(line[2]):
                                add = f"{line[0].strip()}, {line[1]}, {alignment[0]}, {line[2]}, {filter_gene_description(line[3])}|"
                                genes += add
                            elif nextLine != "" and int(alignment[0])+middle_of_align >= int(line[2]) and int(alignment[0])+middle_of_align <= int(nextLine[1]):
                                genes += f"intergenic, {int(line[2])}, {alignment[0]}, {nextLine[1]}|"
                                intergenic += f"intergenic,{int(line[2])}, {alignment[0]}, {nextLine[1]}|"
                if genes == '':
                    genes = "NoGenes"
                gene_column.append(genes)
                intergenic_column.append(intergenic)
    except FileNotFoundError:
        print("         Filtered GFF file not found")
        return 0
    #[['AAATCGAAAGTCACACCGATTACTACAAAGATTAACGTGACCACATTCAAAAACATTAAAAGATATGCTGATAATCTCATCATTTTTTTAGTTCACTTAC', ['0'], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']]]

    #print(genes)

    return gene_column, intergenic_column

def filter_gene_description(info_from_GFF, desired_columns:list = ["gene", "gene_biotype", "Name"]):
    attrs = dict(field.split('=') for field in info_from_GFF.strip(';').split(';'))
    filtered = {key: attrs[key] for key in desired_columns if key in attrs}
    return filtered

def align_to_genome(sequence_matrix, genome_file, kmers, labelled_genomes):
    '''
    Align the sequences to the reference genome.
    The output file is aligned_sequences.txt.
    Parameters:
        sequence_matrix: list - matrix of sequences
        genome_file: str - path to the reference genome file
    Returns:
        int: 1 if the filtering was successful.
    '''
    print("----- Aligning sequences to genome ----------- OUTPUT FILE: alignment_matrix.csv, gene_matrix.csv, intergenic_matrix.csv")

    alignment_results = []

    alignment_matrix = []
    dnas_locations = []
    genes_matrix = []
    intergenic_matrix = []
    sequences_aligned = 0
    if not os.path.exists("ref_seq_index.1.bt2"):
        run(f"bowtie2-build {genome_file} ref_seq_index", shell=True)
    print("         Shape of sequence matrix: ", sequence_matrix.shape)

    # #remove contig separators and take longest remaining segment of DNA
    # for i in range(len(sequence_matrix)):
    #     for j in range(len(sequence_matrix[0])):
    #         if sequence_matrix[i][j] != "NotPresent":
    #             sequences = []
    #             for sequence in sequence_matrix[i][j]:
    #                 try:
    #                     index = sequence.index("___")
    #                     if index <= 49:
    #                         sequence = sequence.split("___")[1]
    #                     elif index > 49:
    #                         sequence = sequence.split("___")[0]
    #                     sequences.append(sequence.replace("_", ""))
    #                 except:
    #                     sequences.append(sequence.replace("_", ""))
    #             sequence_matrix[i][j] = sequences


    print("         Aligning sequences to ref genome.")
    kmers_found = 0
    durations = []
    num_sequences = len(sequence_matrix)*len(sequence_matrix[0])
    sequences_left = num_sequences
    for i in range(len(sequence_matrix)):
        alignment_results = []
        dnas_locations = []
        for j in range(len(sequence_matrix[0])):
            start = time.time()
            if sequence_matrix[i][j][0] != "NotPresent":
                
                
                #only take first match of kmer to genome, other matches are not processed right now
                sequences = ",".join(sequence_matrix[i][j])
                result = run(f"bowtie2 --very-sensitive-local -x ref_seq_index -c {sequences}",shell=True, check=True, capture_output=True)

                #result format keeps changing between genomes???
                results = result.stdout.decode("utf-8").split("\n")
                result_line = []
                for row in results:
                    if row != '':
                        if row[0] not in "@":
                            row = row.split("\t")
                            if int(row[3]) != 0:
                                result_line.append([row[3], row[9]])
                            else:
                                result_line.append([0,"-----NoAlignment-----"])

                alignment_results.append(result_line)
                
                
                dnas_locations.append(result_line)
            else:
                alignment_results.append([[0,"-----NotPresent-----"]])
                dnas_locations.append([[0, "-----NotPresent-----"]])
            sequences_aligned += 1
            sequences_left -= 1
            duration = time.time() - start
            durations.append(duration)
            print(f"\r        Sequences aligned: {sequences_aligned} / {num_sequences} Time per sequence: {round(np.median(durations),3)} s Estimated time: {round((np.median(durations) * sequences_left) / 60, 2)} min", end="", flush=True)
        
        genes, intergenic = find_in_GFF(dnas_locations, "filtered_gff.txt")
        # if any(ele != ['None', '0', '0', '0', '', '', ''] for ele in genes):
        #     kmers_found += 1

        alignment_matrix.append(alignment_results)
        genes_matrix.append(genes)
        intergenic_matrix.append(intergenic)
    print()
    #print("         Kmers found via alignment to reference: ", kmers_found)
    print("\nAlignment Matrix done:")
    print("         Columns: ", len(alignment_matrix[0]))
    print("         Rows: ", len(alignment_matrix))

    genes_matrix = np.array(genes_matrix)
    intergenic_matrix = np.array(intergenic_matrix)
    alignment_matrix = np.array(alignment_matrix)
    print(f"          Matrix shapes genes{genes_matrix.shape}, intergenic{intergenic_matrix.shape}, alignments{alignment_matrix.shape}, sequences{sequence_matrix.shape}")



    return genes_matrix, intergenic_matrix, alignment_matrix

def align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict):

    print("-------------Align To Genomes----------------- OUTPUT: alignments.csv, no_alignments.csv")

    if not os.path.exists("ref_seq_index.1.bt2"):
        run(f"bowtie2-build {ref_genome_file} ref_seq_index", shell=True)

    if not os.path.exists("all_queries.fa"):
        matrix = {}
        fasta_path = "all_queries.fa"
        for i in range(len(sequence_matrix)):
            for j in range(len(sequence_matrix[0])):
                if sequence_matrix[i][j][0] != "NotPresent":
                    matrix[kmers[j], labelled_genomes[i][2:]] = sequence_matrix[i][j]

        with open(fasta_path, "w") as out:
            count = 0
            for (row, col), seqs in matrix.items():
                for seq in seqs:
                    out.write(f">query_{row}_{col}_{count}\n{seq}\n")
                    count += 1
    else:
        print("         all_queries file exists. Skipping query writing.")
    if not os.path.exists("output.sam"):
        call([f"bowtie2 --very-sensitive-local -p 8 -x ref_seq_index -f {fasta_path} -S output.sam"], shell=True)
    else:
        print("         SAM file exists!. Skipping alignment")

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
    print("         Writing k-mers and source genomes to file: kmers_genomes.txt")
    with open("kmers_genomes.txt", "w") as file:
        for kmer in kmers_genomes_dict:
            file.write(f"{kmer}\t{kmers_genomes_dict[kmer]}\n")
    
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

    # inp = input("Extend k-mers from their original genomes again? n to skip")

    # if os.path.exists("./kmers_genomes_sequences_table.csv") and inp == "n":
    #     with open("./kmers_genomes_sequences_table.csv", "r") as file:
    #         lines = file.readlines()
    #         col_names = lines[0]
    #         lines = lines[2:]
    #         for line in lines:
    #             sequence_row = []
    #             for item in line.strip().split(",")[1:]:
    #                 #print(item)
    #                 try:
    #                     sequence_row.append(item.split(";"))
    #                 except:
    #                     sequence_row.append([item])
    #         sequence_matrix.append(sequence_row)
    #     sequence_matrix = np.array(sequence_matrix)
    # else:
    #take subset of data
    if reduce_genomes_to != 0 and reduce_genomes_to < len(labelled_genomes):
        labelled_genomes = labelled_genomes[:math.floor(reduce_genomes_to/2)] + labelled_genomes[-math.ceil(reduce_genomes_to/2):]
    if reduce_kmers_to != 0 and reduce_kmers_to < len(kmers_genomes_dict.keys()):
        kmers = list(kmers_genomes_dict.keys())[:reduce_kmers_to]
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

        # print(ref_genome_row)

        # #format ref genome row and write to file
        # for sequence_set in ref_genome_row:
        #     if len(sequence_set) < 2:
        #         ref_genome_row_string += "," + str(sequence_set).replace("[", "").replace("]", "").replace("'", "")
        #     else:
        #         ref_genome_row_string += "," + str(sequence_set).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
        # file.write("\nREFERENCE " + ref_genome_row_string)
        
        #format genome rows and write to file
        shape = sequence_matrix.shape
        for i in range(shape[0]):
            file.write("\n" + labelled_genomes[i])
            for j in range(shape[1]):
            
                sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "")
                if len(sequence_matrix[i][j]) > 1:
                    sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
                file.write("," + sequences)


    genes_matrix, intergenic_matrix, alignment_matrix = align_to_genome2(sequence_matrix, ref_genome_file, kmers, labelled_genomes, GFF_file, kmers_coeffs_genomes_dict)          #align the sequences to the reference genome

    #print(alignment_matrix[0])

    #GENE_MATRIX[0]
    #[['TATTCGTTCAATGATCCGAGGGAAACTTGGGGATTGGATCTTAAGTATTTTGGAAAACAAATATGACTTAAATCACCTGGACGCGATGAAATTATATCAA', 'gene', '8671', '10231', '10377', 'ID=gene-PA0007;Name=PA0007;gbkey=Gene;gene_biotype=protein_coding;locus_tag=PA0007'], 
    # ['TATTCGTTCAATGATCCGAGGGAAACTTGGGGATTGGATCTTAAGTATTTTGGAAAACAAATATGACTTAAATCACCTGGACGCGATGAAATTATATCAA', 'CDS', '8671', '10231', '10377', 'ID=cds-AAG03397.1;Parent=gene-PA0007;Dbxref=NCBI_GP:AAG03397.1;Name=AAG03397.1;Note=Product name confidence: Class 4 (Homologs of previously reported genes of unknown function%2C or no similarity to any previously reported sequences);gbkey=CDS;locus_tag=PA0007;product=hypothetical protein;protein_id=AAG03397.1;transl_table=11']]


    #ALIGNMENT_MATRIX[0]
    #[['0', 'GATTCTACCTTCTAAACACGTAATAAAGTATCTAAAAAAATTAAAAGAGAAAACATTAAAAGAACAATTTTTAA', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'], 
    # ['10231', 'TATTCGTTCAATGATCCGAGGGAAACTTGGGGATTGGATCTTAAGTATTTTGGAAAACAAATATGACTTAAATCACCTGGACGCG', 'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII']

    # locations_matrix = []
    # original_sequence_matrix = []
    # aligned_sequence_matrix = []

    # with open("gene_matrix.csv", "w") as file:
    #     # sequence_matrix = np.transpose(sequence_matrix)
    #     # for col in sequence_matrix:
    #     #     col_header = str(set(col)).replace(",", ";")
    #     #     file.write(f"{col_header},")
    #     # file.write("\n")
    #     file.write("GENES")
    #     for i in range(len(kmers)):
    #         #add first row with kmers and their parameters in genomes.
    #         file.write(f", {kmers[i]}/{reverse_complement(kmers[i])} {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {kmers_coeffs_genomes_dict[kmers[i]][0]}")
    #     file.write("\n")
    #     for i in range(len(genes_matrix)):
    #         file.write(f"{labelled_genomes[i]},")
    #         for j in range(len(genes_matrix[i])):
    #             if len(genes_matrix[i][j]) > 1:
    #                 if genes_matrix[i][j] != '':
    #                     string = genes_matrix[i][j].replace(",", " ")
    #                     file.write(f"({string}),")
    #                 else:
    #                     file.write("---------,")
    #         file.write("\n")

    # for i in range(len(genes_matrix)):
    #     row1 = []
    #     row2 = []
    #     row3 = []
    #     for j in range(len(genes_matrix[i])):
    #         if genes_matrix[i][j] != "NoGenes":
    #             genes = genes_matrix[i][j].split("|")
    #             locations = ""
    #             if len(genes) > 2:
    #                 for gene in genes:
    #                     if gene != '':
    #                         locations += gene.split(",")[2] + "|"
    #             else:
    #                 locations = genes[0].split(",")[2]
    #             row1.append(locations)
    #         else:
    #             row1.append("NoGenes")
    #         row2.append(str(sequence_matrix[i][j]))
    #         alignments = "|".join([f"{elem[0]};{elem[1]}" for elem in alignment_matrix[i][j]])
    #         row3.append(alignments)
    #     locations_matrix.append(row1)
    #     original_sequence_matrix.append(row2)
    #     aligned_sequence_matrix.append(row3)

    # with open("alignment_matrix.csv", "w") as file:
    #     # sequence_matrix = np.transpose(sequence_matrix)
    #     # for col in sequence_matrix:
    #     #     col_header = str(np.unique(col)).replace(",",";")
    #     #     file.write(f"{col_header},")
    #     # file.write("\n")
    #     file.write("ALIGNMENTS")
    #     for i in range(len(kmers)):
    #         #add first row with kmers and their parameters in genomes.
    #         file.write(f", {kmers[i]}/{reverse_complement(kmers[i])} {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {kmers_coeffs_genomes_dict[kmers[i]][0]}")
    #     file.write("\n")
    #     for i in range(len(aligned_sequence_matrix)):
    #         file.write(f"{labelled_genomes[i]},")
    #         for j in range(len(aligned_sequence_matrix[i])):
    #             file.write(f"{aligned_sequence_matrix[i][j]},")
    #         file.write("\n")
    # # with open("expanded_sequence_matrix.csv", "w") as file:
    # #     # sequence_matrix = np.transpose(sequence_matrix)
    # #     # for col in sequence_matrix:
    # #     #     col_header = str(set(col)).replace(",", ";")
    # #     #     file.write(f"{col_header},")
    # #     # file.write("\n")
    # #     file.write("EXPANDED,")
    # #     for kmer in kmers:
    # #         file.write(f"{kmer},")
    # #     file.write("\n")
    # #     for i in range(len(original_sequence_matrix)):
    # #         file.write(f"{labelled_genomes[i]},")
    # #         for j in range(len(original_sequence_matrix[i])):
    # #             sequences = str(original_sequence_matrix[i][j]).replace(",", ";")
    # #             file.write(f"{sequences},")
    # #         file.write("\n")
   
    # with open("intergenic_matrix.csv", "w") as file:
    #     # sequence_matrix = np.transpose(sequence_matrix)
    #     # for col in sequence_matrix:
    #     #     col_header = str(set(col)).replace(",", ";")
    #     #     file.write(f"{col_header},")
    #     # file.write("\n")
    #     file.write("INTERGENIC")
    #     for i in range(len(kmers)):
    #         #add first row with kmers and their parameters in genomes.
    #         file.write(f", {kmers[i]}/{reverse_complement(kmers[i])} {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {kmers_coeffs_genomes_dict[kmers[i]][0]}")
    #     file.write("\n")
    #     for i in range(len(intergenic_matrix)):
    #         file.write(f"{labelled_genomes[i]},")
    #         for j in range(len(intergenic_matrix[i])):
    #             if len(intergenic_matrix[i][j]) > 1:
    #                 if intergenic_matrix[i][j][0] != "None":
    #                     file.write(f"{intergenic_matrix[i][j][0]};{intergenic_matrix[i][j][1]};{intergenic_matrix[i][j][2]},")
    #                 else:
    #                     file.write("-----------,")
    #         file.write("\n")
    # with open("locations_matrix.csv", "w") as file:
    #     # sequence_matrix = np.transpose(sequence_matrix)
    #     # for col in sequence_matrix:
    #     #     col_header = str(set(col)).replace(",", ";")
    #     #     file.write(f"{col_header},")
    #     # file.write("\n")
    #     file.write("LOCATIONS")
    #     for i in range(len(kmers)):
    #         #add first row with kmers and their parameters in genomes.
    #         file.write(f", {kmers[i]}/{reverse_complement(kmers[i])} {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {kmers_coeffs_genomes_dict[kmers[i]][0]}")
    #     file.write("\n")
    #     for i in range(len(locations_matrix)):
    #         file.write(f"{labelled_genomes[i]},")
    #         for j in range(len(locations_matrix[i])):
    #             if locations_matrix[i][j]:
    #                 file.write(f"{locations_matrix[i][j]},")
    #             else:
    #                 file.write(f"--------,")
    #         file.write("\n")

    

    #make readable excel file. rows as genomes, columns as kmers, sequences as values
    #add kmer coefficients to the file

import re

def find_all_indexes(text, substring):
    # Use re.finditer to find all matches of the substring
    return [match.start() for match in re.finditer(re.escape(substring), text)]

def reverse_complement(dna_sequence):
    """
    Returns the reverse complement of a DNA sequence.
    Parameters:
        dna_sequence (str): The input DNA sequence (e.g., "ATGC")
    Returns:
        str: The reverse complement of the DNA sequence.
    """
    # Define the complement mapping
    complement = {
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C'
    }
    # Generate the reverse complement
    reverse_comp = ''.join(complement[base] for base in reversed(dna_sequence.upper()))
    return reverse_comp

#UNUSED, TOO SLOW
def straighten_DNA(genome):
    straight_dna = ""
    with open(genome, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line[0]  == ">":
                straight_dna += "___"
            else:
                straight_dna += line.strip()
    with open("new_genome.txt", "w") as file:
        file.write(straight_dna)
    return straight_dna

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
#UNUSED, TOO SLOW
def find_kmers_in_genome(kmer, genome, add_length: int = 50):
    '''
    Find kmers in the genome using the indexed kmer file and the kmers and coefficients file.
    The output file is pheno_kmers_ref_genome.txt.
    Parameters:
        indexed_kmer_file: str - path to the indexed kmer file
        kmers_coeffs_file: str - path to the kmers and coefficients file
        output_file: str - path to the output file default: pheno_kmers_ref_genome.txt
    Creates files:
        - pheno_kmers_ref_genome.txt
    Returns:
        int: 0 if the output file already exists or no kmers were found, 1 if the filtering was successful.
    '''
    
    straight_dna = straighten_DNA(genome)
    
    locations = []
    indexes = find_all_indexes(straight_dna, kmer.lower())
    #print(indexes)
    for index in indexes:
        locations.append(index)
    
    indexes = find_all_indexes(straight_dna, reverse_complement(kmer).lower())
    for index in indexes:
        locations.append(index)
    #print(find_all_indexes(straight_dna, (kmer.lower())[::-1]))
    longer_sequences = []
    for location in locations:
        #print("Location: ", location)
        if location != []:
            if location > add_length:
                if location < len(straight_dna) - len(kmer) - add_length:
                    longer_sequences.append(straight_dna[location-add_length:location]+str(straight_dna[location:location+len(kmer)]).upper()+straight_dna[location+len(kmer):location+add_length])
                else:
                    longer_sequences.append(straight_dna[location-add_length:location+len(kmer)])
            else:
                longer_sequences.append(straight_dna[location:location+len(kmer)+add_length])
    if longer_sequences == []:
        #print("No longer sequences found")
        longer_sequences.append(genome)
    
    return longer_sequences

def combine_kmers_by_locations(significant_kmers):
    #make file with combined kmers
    #only combine kmers that are close enough to eachother
    #fill in gaps with "_"
    pass

# anda ka genoomid kus sellised kmeerid esinesid täpsemaks analüüsiks
    #joondada ka nendes genoomides pikemaid DNA juppe kmeeride ümber referents genoomi vastu (valideerimine ja täpsemate mutatsioonide leidmine)

# mis ei seostu millegiga
# kmeeride liitmine
# erki artikklist valideerida
# kmeerid asukoha järgi kokku panna




#võta genoom
#indekseeri
#võta kmeerid
#kui kmeer esineb genoomis
#leia kmeeri asukoht genoomis
#võta kmeeri ümber 1000bp
#kui kõik kmeerid on kontrollitud kustuta indekseeritud fail



