from subprocess import call, run
import numpy as np
import pandas as pd
import os
import time
import re




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

    genes: list[list] = []
    intergenic: list[list] = []

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
    
    genes: list[list] =  []
    intergenic: list[list] = []
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
        return [],[]

#USED
def filter_gene_description(info_from_GFF, desired_columns:list = ["gene", "gene_biotype", "Name"]):
    attrs = dict(field.split('=') for field in info_from_GFF.strip(';').split(';'))
    filtered = {key: attrs[key] for key in desired_columns if key in attrs}
    return filtered

def find_in_GFF(dnas_with_locations, filtered_GFF_file):
    
    gene_column: list[str] =  []
    intergenic_column: list[str] = []
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
        return [],[]
    #[['AAATCGAAAGTCACACCGATTACTACAAAGATTAACGTGACCACATTCAAAAACATTAAAAGATATGCTGATAATCTCATCATTTTTTTAGTTCACTTAC', ['0'], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']], ['NotFound', [0], ['None', '0', '0', '0', '', '', ''], ['None', '0', '0', '0', '', '', '']]]

    #print(genes)

    return gene_column, intergenic_column

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



