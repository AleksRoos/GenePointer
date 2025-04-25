from subprocess import call
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
        call(["phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier + " -l " + str(kmer_length)], shell=True)
    # os.chdir("..")           #change directory back to the original folder
        print("         Finished executing phenotypeseeker")

    # else:
    #     print("         Found")
    
    try:
        if classifier == "log":
            classifier = "log_reg"
        with open("k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt", "r") as kmers_coeffs:      #open k-mers and coefficients file with significant kmers found by phenotypeseeker
            lines = kmers_coeffs.readlines()
            if len(lines) == 0:
                print("        No significant kmers found")
                return 0

            kmers = []
            coeffs = []
            whole_lines = []
            for line in lines[1:]:                  #take only kmers with non-zero coefficients from the file
                if float(line.split("\t")[1]) != 0:
                    kmers.append(line.split("\t")[0])
                    coeffs.append(line.split("\t")[1])
                    whole_lines.append(line.strip().split("\t"))
            print("         Number of significant kmers: ", len(kmers))
            with open("filtered_kmers_and_coeffs.txt", "w") as kmercoefffile:
                i = 0
                for kmer in kmers:
                    kmercoefffile.write(f"{kmer}\t{coeffs[i]}\n")
                    i += 1
            with open("significant_lines_from_kmer_coeff_file.txt", "w") as file:
                for line in whole_lines:
                    genomes = line[3] #add genomes cleanly to end of line or make clean list of them
                    genomes = genomes.split(" ")[2:]
                    string_genomes = ""
                    for item in genomes:
                        string_genomes += str(item) + ' '
                    #print(string_genomes)
                    file.write(f"{line[0]}\t{line[1]}\t{line[2]}\t{string_genomes}\n")
    except FileNotFoundError:
        print("        k-mers and coefficients file not found: No kmers found by phenotypeseeker")
        return 0
    return 1

def index_genome(genome_path, output_file = "ref_gnome_kmer_locs.txt", kmer_length = 13, repeat:bool = False):
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

def find_ref_genome_kmers(indexed_kmer_file: str,kmers_coeffs_file: str, output_file: str="pheno_kmers_ref_genome.txt", repeat:bool = False):
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
    kmers = []
    coeffs = []

    try:
        with open(kmers_coeffs_file, "r") as kmers_coeffs:
            lines = kmers_coeffs.readlines()
            for line in lines:
                line = line.strip().split("\t")
                kmers.append(line[0])
                coeffs.append(line[1])
    except FileNotFoundError:
        print("        Kmers and coefficients file not found")
        return 0
        
    print("----- Searching for kmers in ref genome---------- OUTPUT FILE: ", output_file)        #Zipping kmers and coefficients together
    if kmers == []:
        print("        No kmers found")
        return 0   
    
    
    kmers_with_coeffs = zip(kmers, coeffs)
    print("         Sorting kmers by coefficients")     # Sort kmers by coefficients, greatest to lowest
    sorted_kmers = []
    sorted_kmers_with_coeffs = sorted(kmers_with_coeffs, key=lambda x: x[1])
    for line in sorted_kmers_with_coeffs:
        sorted_kmers.append((line[0],line[1]))

    print("     Finding kmers in indexed file. May take a while...")       # Find kmers in indexed reference genome file
    shared_kmers = []
    matching_kmers = []
    try:
        with open(indexed_kmer_file, "r") as file:          #open indexed ref genome
            lines = file.readlines()
            for i in range(len(lines)):
                for line in sorted_kmers:
                    kmer = line[0]                   
                    if kmer in lines[i]:                    #check if a phenotype kmer is in this line of indexed ref genome 
                        matching_kmers.append(line)
                        index_lines = lines[i:i+int(lines[i].strip().split()[1])+1]         #add all locations of the kmer to the shared kmers list
                        for line in index_lines:
                            shared_kmers.append(line)
                        break
    except FileNotFoundError:
        print("        Indexed kmer file not found")
        return 0

    print("     Number of kmers found in indexed file: ", len(matching_kmers))          # Write phenotype associated kmers found in reference genome to file
    with open("pheno_kmers_ref_genome_with_coeffs.txt", "w") as file:
        for line in matching_kmers:
            file.write(f"{line[0]}\t{line[1]}\n")
    with open(output_file, "w") as file:
        for line in shared_kmers:
            file.write(line)
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

def find_genes(pheno_kmers_in_ref_genome_file, filtered_GFF_file, min_mismatches = 0,output_file="genes.txt", repeat:bool = False):
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
    print("         Extracting kmers and their locations")
    try:

        with open(pheno_kmers_in_ref_genome_file, "r") as significant_kmer_file:
            lines = significant_kmer_file.readlines()
            for i in range(len(lines)):
                kmer = []
                locations = []
                line = lines[i].strip().split("\t")
                if len(line) == 2:                          #start  parsing on lines that begin with the kmer. 2 items on the line
                    kmer = line[0]
                    j = 1                          #kmer is the first element
                    while j < int(line[1])+1:         #iterate over line[1] lines under kmer. there are lines with locations and mismatches under the kmer 
                        sub_line = lines[i+j].strip().split("\t")        #a line under kmer line
                        if int(sub_line[-1]) <= min_mismatches:              #if kmer mismatches is less then minimum save the location
                            locations.append(sub_line[2])
                        j += 1
                    if len(locations) > 0:                  #if there are locations, add the kmer and its locations to the list
                        kmers_with_locations.append([kmer, *locations])
                    else:
                        kmers_with_locations.append([kmer, "No locations found"])
                else:
                    continue
    except FileNotFoundError:
        print("         Significant kmers file not found")
        return 0
    
    
    print("         Number of kmers with locations: ", len(kmers_with_locations))
    print("         Identifying kmer parent elements in reference genome")

    genes = []
    intergenic = []
    used_kmers = []   
    kmers_without_locations = []         
    try:    
        with open(filtered_GFF_file, "r") as file:          #open filtered GFF file
            lines = file.readlines()
            
            for kmerline in kmers_with_locations:
                if kmerline[1] == "No locations found":          #if there are no locations for the kmer, add it to the kmers_without_locations list
                    kmers_without_locations.append(kmerline)
                    kmers_with_locations.pop(kmers_with_locations.index(kmerline))
                    continue
                else:
                    for j in range(len(lines)):
                        line = lines[j].strip().split("\t")
                        nextLine = ""
                        if j < len(lines)-1:
                            nextLine = lines[j+1].strip().split("\t")
                        if line[0] == "region":
                            continue
                    #print(line[:3])
                        else:
                            for i in range(1,len(kmerline)):
                                #print("range ", range(1, len(kmerline)))
                                if int(kmerline[i]) >= int(line[1]) and int(kmerline[i]) <= int(line[2]):
                                    #print("line ", line)
                                    add = [kmerline[0], line[0], line[1], kmerline[i], line[2], line[3]]
                                    if add not in genes:
                                        genes.append(add)
                                    used_kmers.append(kmerline[0])
                                elif nextLine != "" and int(kmerline[i]) >= int(line[2]) and int(kmerline[i]) <= int(nextLine[1]):
                                    intergenic.append([kmerline[0],"intergenic",int(line[2]), kmerline[i], nextLine[1]," "])
                                    used_kmers.append(kmerline[0])
    except FileNotFoundError:
        print("         Filtered GFF file not found")
        return 0
   

    for item in kmers_with_locations:
        if item[0] not in used_kmers:
            if len(item)>1:
                kmers_without_locations.append([item[0], item[1:]])               

    print(kmers_without_locations)

    print("         ELEMENTS found: ", len(genes))
    print("         INTERGENIC found: ", len(intergenic))
    print("         UNIDENTIFIED kmers: ", len(kmers_without_locations))


    with open(output_file, "w") as gene_file:
        i = 1
        for line in genes:
            line[5] = filter_id_line(line[5])
            gene_file.write(f"{line[0]},{line[1]},{line[2]},{line[3]},{line[4]},{line[5]}\n")
            i += 1
    with open("intergenic.txt", "w") as intergenic_file:
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


def show_source_genomes(significant_kmers):    
    #make file with genomes, kmers, locations in genome, larger DNA piece around kmer

    with open(significant_kmers, "r") as file:
        lines = file.readlines()
        data_frame = []
        genomes = ""
        for line in lines:
            line = line.strip().split("\t")
            data_frame.append([line[0].strip(), line[1].strip(), line[2].strip(), line[3].strip()])
            genomes += " " + line[3].strip()
        genomes = genomes.strip().split(" ")
        genomes = set(genomes)
        print(data_frame[0])
    pass

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