from subprocess import call, run
import numpy as np
import os


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
        #call(["phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier + " -l " + str(kmer_length)], shell=True)
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

def filter_GFF_file(GFF_file, output_file="filtered_gff.txt", repeat:bool = False):
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
            for line in lines:              #taking only the lines with the elements of interest: element type, start, end, id
                line = line.strip().split("\t")
                if len(line) < 9:
                    continue
                elements.append(line[2].strip())
                start_end.append((line[3].strip(), line[4].strip()))
                ids.append(line[8].strip())
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
    print("         Finding kmers in indexed file. May take a while...")       # Find kmers in indexed reference genome file
    #make dictionary instead
    shared_kmers = {}
    try:
        #take


        with open(indexed_kmer_file, "r") as file:          #open indexed ref genome
            print("         Reading indexed kmer file: ", indexed_kmer_file)
            lines = file.readlines()


            kmer_count = 0
            for i in range(len(lines)):
                print("\r         Searching lines: ", i, " / ", len(lines), "    Kmers found: ", kmer_count, end="", flush=True)
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
                        kmer_count += 1
                        index_lines = lines[i+1:i+index_count+1]         #add all locations of the kmer to the shared kmers list
                        for index in index_lines:
                            index = index.strip().split("\t")
                            if int(index[-1]) > max_mismatches:
                                continue
                            if not shared_kmers.get(kmer):
                                shared_kmers[kmer] = [kmers_coeffs[j][1]]
                            shared_kmers[kmer].append(index[2])
                        break
    except FileNotFoundError: 
        print("        Indexed kmer file not found")
        return 0

    print(f"     Number of kmers found in indexed file (max_mismatch: {max_mismatches}): ", len(list(shared_kmers.keys())))          # Write phenotype associated kmers found in reference genome to file
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
        
    return ";".join(new_items)

def find_genes(pheno_kmers_in_ref_genome_file, filtered_GFF_file,output_file="genes.txt"):
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
    print("----- Finding genes ------------------ OUTPUT FILE: ", output_file)

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


    genes, intergenic = find_in_GFF(kmers_with_locations, filtered_GFF_file)

    print("         UNIDENTIFIED kmers: ", len(kmers_without_locations))
    print("         ELEMENTS found: ", len(genes))
    print("         INTERGENIC found: ", len(intergenic))

    with open("genes.csv", "w") as gene_file:
        i = 1
        for line in genes:
            line[5] = filter_id_line(line[5])
            gene_file.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]}\n")
            i += 1
    with open("intergenic.csv", "w") as intergenic_file:
        i = 1
        for line in intergenic:
            intergenic_file.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]}\n")
            i += 1
    with open("unidentified_kmers.txt", "w") as unidentified_kmers_file:
        i = 1
        for line in kmers_without_locations:
            locations = ""
            for item in line[1:]:
                locations += str(item) + " "
            unidentified_kmers_file.write(f"{line[0]},{locations}\n")
            i += 1
    print("----- Finished finding genes -----")
    return 1

def find_in_GFF(dnas_with_locations, filtered_GFF_file):
    
    genes = []
    intergenic = []
    try:    
        with open(filtered_GFF_file, "r") as file:          #open filtered GFF file
            lines = file.readlines()
            
            for dnaloc in dnas_with_locations:
                for j in range(len(lines)):
                    line = lines[j].strip().split("\t")
                    nextLine = ""
                    if j < len(lines)-1:
                        nextLine = lines[j+1].strip().split("\t")
                    if line[0] == "region":
                        continue
                #print(line[:3])
                    else:
                        for i in range(len(dnaloc[1])):
                            #print("range ", range(1, len(kmerline)))
                            if int(dnaloc[1][i]) >= int(line[1]) and int(dnaloc[1][i]) <= int(line[2]):
                                #print("line ", line)
                                add = [dnaloc[0], line[0], line[1], dnaloc[1][i], line[2], line[3]]
                                #print(add)
                                if add not in genes:
                                    genes.append(add)
                            elif nextLine != "" and int(dnaloc[1][i]) >= int(line[2]) and int(dnaloc[1][i]) <= int(nextLine[1]):
                                intergenic.append([dnaloc[0],"intergenic",int(line[2]), dnaloc[1][i], nextLine[1]," "])
    except FileNotFoundError:
        print("         Filtered GFF file not found")
        return 0

    
    
    return genes, intergenic


def align_to_genome(sequence_matrix, genome_file):
    '''
    Align the sequences to the reference genome.
    The output file is aligned_sequences.txt.
    Parameters:
        sequence_matrix: list - matrix of sequences
        genome_file: str - path to the reference genome file
    Returns:
        int: 1 if the filtering was successful.
    '''
    print("----- Aligning sequences to genome ----------- OUTPUT FILE: aligned_sequences.txt")

    alignment_results = []

    alignment_matrix = []
    dnas_locations = []
    genes_matrix = []
    intergenic_matrix = []
    sequences_aligned = 0
    if not os.path.exists("ref_seq_index.1.bt2"):
        run(f"bowtie2-build {genome_file} ref_seq_index", shell=True)
    print("         Shape of sequence matrix: ", sequence_matrix.shape)
    for i in range(len(sequence_matrix)):
        alignment_results = []
        dnas_locations = []
        for j in range(len(sequence_matrix[0])):
            if sequence_matrix[i][j] != ["No longer sequences found"]:
                try:
                    sequence_matrix[i][j] = str(sequence_matrix[i][j][0]).split("___")[1]
                except:
                    pass
                if ";" in sequence_matrix[i][j]:
                    sequence_matrix[i][j] = sequence_matrix[i][j].split(";")[0]
                #only take first match of kmer to genome, other matches are not processed right now
                result = run(f"bowtie2 -x ref_seq_index -c {str(sequence_matrix[i][j][0]).lower()}",shell=True, check=True, capture_output=True)

                results = result.stdout.decode("utf-8").split("\t")
                alignment_results.append([results[2+10], results[8+10], results[9+10]])
                dnas_locations.append([results[8+10], [results[2+10]]])
                #print(dnas_locations[-1])
                #print(alignment_results[-1])
            else:
                alignment_results.append(["No sequence in reference genome with kmer"])
            sequences_aligned += 1
            print("\r        Sequences aligned: ", sequences_aligned, " / ", len(sequence_matrix)*len(sequence_matrix[i]), end="", flush=True)
        genes, intergenic = find_in_GFF(dnas_locations, "filtered_gff.txt")


        alignment_matrix.append(alignment_results)
        genes_matrix.append(genes)
        intergenic_matrix.append(intergenic)
    print("\nAlignment Matrix done:")
    print("     Columns: ", len(alignment_matrix[0]))
    print("     Rows: ", len(alignment_matrix))

    genes_matrix = np.transpose(np.array(genes_matrix))
    intergenic_matrix = np.transpose(np.array(intergenic_matrix))   

    with open("alignment_matrix.csv", "w") as file:
        for i in range(len(alignment_matrix)):
            for j in range(len(alignment_matrix[i])):
                if len(alignment_matrix[i][j]) > 1:
                    file.write(f"{alignment_matrix[i][j][0]};{alignment_matrix[i][j][1]};{alignment_matrix[i][j][2]},")
                else:
                    file.write(f"{alignment_matrix[i][j]},")
            file.write("\n")
    with open("gene_matrix.csv", "w") as file:
        for i in range(len(genes_matrix)):
            for j in range(len(genes_matrix[i])):
                if len(genes_matrix[i][j]) > 1:
                    file.write(f"({genes_matrix[i][j][0].strip()};{genes_matrix[i][j][1].strip()};{genes_matrix[i][j][3].strip()}),")
            file.write("\n")
    with open("intergenic_matrix.csv", "w") as file:
        for i in range(len(intergenic_matrix)):
            for j in range(len(intergenic_matrix[i])):
                if len(intergenic_matrix[i][j]) > 1:
                    file.write(f"{intergenic_matrix[i][j][0]};{intergenic_matrix[i][j][1]};{intergenic_matrix[i][j][2]},")
            file.write("\n")

    #return location_matrix


def find_genes_alignment(significant_kmers, species, antibiotic, ref_genome_file = "Not added"):
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

    print("         Processing significant k-mer file.")
    with open(significant_kmers, "r") as file:
        lines = file.readlines()
        genomes = ""
        for line in lines:
            line = line.strip().split("\t")
            data_frame.append([line[0].strip(), line[1].strip(), line[2].strip(), line[3].strip()])
            genomes += " " + line[3].strip()
            kmers_coeffs_genomes_dict[line[0]] = [line[1], line[2]]
            kmers_genomes_dict[line[0]] = line[3].strip().split(" ")
        genomes = genomes.strip().split(" ")
        genomes = set(genomes)
    labelled_genomes = []
    genome_folder_genome_names = os.listdir(f"Genomes/" + genome_folder)
    for genome in genomes:
        for genome_file in genome_folder_genome_names:
            if genome in genome_file:
                labelled_genomes.append(genome_file[:-4])
                break
    #K-mers with higher coefficients are first
    print("         K-mers sorted for coefficients.")
    data_frame = sorted(data_frame, key=lambda x: x[1])          #sort data frame by kmer

    print("         Writing k-mers and source genomes to file.")
    with open("kmers_genomes.txt", "w") as file:
        for kmer in kmers_genomes_dict:
            file.write(f"{kmer}\t{kmers_genomes_dict[kmer]}\n")

    #print(kmer_coeff_genomes_dict)
    longer_sequences = []
    kmers = list(kmers_genomes_dict.keys())
    sequence_matrix = []
    kmer_count = 0
    print("         Expanding k-mers from their original genomes")

    for i in range(len(kmers)):
        columns = []
        for j in range(len(labelled_genomes)):
            if labelled_genomes[j][2:] in kmers_genomes_dict[kmers[i]]:
                genome_file = "Genomes/" + genome_folder + "/" + labelled_genomes[j] + ".fna"
                #print("Genome file: ", genome_file)
                longer_sequences = find_kmers_in_genome(kmers[i], genome_file)
                columns.append(longer_sequences)
            else:
                columns.append(["No longer sequences found"])
        kmer_count += 1
        print(f"\r      K-mers searched:  {kmer_count} / {len(kmers)}", end="", flush=True)
        sequence_matrix.append(columns)
    
    kmer_count_res_sus_dict = {}
    for i in range(len(kmers)):
        genome_count = 0
        resistant_count = 0
        susceptible_count = 0
        for j in range(len(labelled_genomes)):
            if labelled_genomes[j][2:] in kmers_genomes_dict[kmers[i]]:
                genome_count += 1
                print(labelled_genomes[j][0])
                if int(labelled_genomes[j][0]) == 1:
                    resistant_count += 1
                else:
                    susceptible_count += 1
        kmer_count_res_sus_dict[kmers[i]] = f"T:{genome_count}; R:{resistant_count}/S:{susceptible_count}({round(resistant_count / genome_count, 2)*100}%)"

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
            file.write(f", {kmers[i]} Counts(g; R; S): {kmer_count_res_sus_dict[kmers[i]]} Coefficient: {kmers_coeffs_genomes_dict[kmers[i]][0]}")
            #make row for ref genome separately
            ref_genome_row.append(find_kmers_in_genome(kmers[i], ref_genome_file))
        
        #format ref genome row and write to file
        for sequence_set in ref_genome_row:
            if len(sequence_set) < 2:
                ref_genome_row_string += "," + str(sequence_set).replace("[", "").replace("]", "").replace("'", "")
            else:
                ref_genome_row_string += "," + str(sequence_set).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
        file.write("\nREFERENCE " + ref_genome_row_string)
        
        #format genome rows and write to file
        for i in range(len(labelled_genomes)):
            file.write("\n" + labelled_genomes[i])
            for j in range(len(kmers)):
                
                if sequence_matrix[i][j] != ["No longer sequences found"]:
                    sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "")
                    if len(sequence_matrix[i][j]) > 1:
                        sequences = str(sequence_matrix[i][j]).replace("[", "").replace("]", "").replace("'", "").replace(",", ";")
                    file.write("," + sequences)
                else:
                    file.write(",")

    align_to_genome(sequence_matrix, ref_genome_file)          #align the sequences to the reference genome

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
    locations = []
    indexes = find_all_indexes(straight_dna, kmer.lower())
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
        longer_sequences.append("No longer sequences found")
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



