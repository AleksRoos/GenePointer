from subprocess import call
import os


def find_significant_kmers(data_pheno_path: str, classifier: str = "log", kmer_length: int = 13):
    print("----- Looking for significant kmers -----")
    antibiotic = ""
    with open(data_pheno_path, "r") as dataphenofile:           #open data.pheno file
        lines = dataphenofile.readlines()
        if len(lines) < 2:
            print("        No data in file")
            return 0
        #print("DataPheno Header: ", lines[0].strip().split("\t"))
        antibiotic = lines[0].strip().split()[2]           #get antibiotic name from first line of file

    print("         Looking for file ./k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt")
    if not os.path.exists("./k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt"):
        print("         EXECUTING: phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier+ " -l " + str(kmer_length))             #if the k-mers and coefficients file does not exist, call phenotypeseeker to find significant kmers
        
        # try:
        #     call(["mkdir " + antibiotic], shell=True)           #create a directory for the antibiotic
        # except FileExistsError:
        #     print("         Directory already exists")
        # os.chdir(antibiotic)           #change directory to the antibiotic folder
        # print("         Changing directory to ", antibiotic)
        call(["phenotypeseeker modeling " + data_pheno_path + " -w -bc " + classifier + " -l " + str(kmer_length)], shell=True)
        # os.chdir("..")           #change directory back to the original folder
        print("         Finished executing phenotypeseeker")

    else:
        print("         Found")
    try:
        with open("k-mers_and_coefficients_in_" + classifier + "_model" + "_" + antibiotic.capitalize() + ".txt", "r") as kmers_coeffs:      #open k-mers and coefficients file with significant kmers found by phenotypeseeker
            lines = kmers_coeffs.readlines()
            if len(lines) == 0:
                print("        No significant kmers found")
                return 0

            kmers = []
            coeffs = []
            for line in lines[1:]:                  #take only kmers and non-zero coefficients from the file
                if float(line.split("\t")[1]) != 0:
                    kmers.append(line.split("\t")[0])
                    coeffs.append(line.split("\t")[1])
            print("         Number of significant kmers: ", len(kmers))
            with open("filtered_kmers_and_coeffs.txt", "w") as kmercoefffile:
                i = 0
                for kmer in kmers:
                    kmercoefffile.write(f"{kmer}\t{coeffs[i]}\n")
                    i += 1
    except FileNotFoundError:
        print("        k-mers and coefficients file not found: No kmers found by phenotypeseeker")
        return 0
    return 1



def index_reference_genome(reference_genome_path, output_file = "ref_gnome_kmer_locs.txt", kmer_length = 13):
    print("----- Indexing reference genome -----")
    if os.path.exists(output_file):
        print("         Index file already exists. Stopping indexing")
        return 0
    if not os.path.exists(reference_genome_path):
        print("         Reference genome file not found")
        return 0
    if not os.path.exists(f"out_{kmer_length}.index"):
        print("         Index file not found")
        print("         Indexing reference genome")
        call([f"glistmaker {reference_genome_path} -w {kmer_length} --index"], shell=True)           #index reference genome, locations of each kmer in the ref genome
    
    call([f"glistquery out_{kmer_length}.index --locations >> {output_file}"], shell=True)    #write the index file into a text file for reading later
    print("----- Finished indexing reference genome -----")

def filter_GFF_file(GFF_file, output_file="filtered_gff.txt"):
    # if os.path.exists(output_file):
    #     print("        Output file already exists. Stopping filtering")
    #     return 0
    print("----- Reducing GFF file -----")
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
    with open("filtered_gff.txt", "w") as file:          #write the filtered GFF file to a new file
        for i in range(len(elements)):
            file.write(f"{elements[i]}\t{start_end[i][0]}\t{start_end[i][1]}\t{ids[i]}\n")
    print("----- Finished reducing GFF file -----")

def find_ref_genome_kmers(indexed_kmer_file: str,kmers_coeffs_file: str, output_file: str="pheno_kmers_ref_genome.txt"):
    
    kmers = []
    coeffs = []

    with open(kmers_coeffs_file, "r") as kmers_coeffs:
        lines = kmers_coeffs.readlines()
        for line in lines:
            line = line.strip().split("\t")
            kmers.append(line[0])
            coeffs.append(line[1])

        
    print("----- Searching for kmers in ref genome-----")        #Zipping kmers and coefficients together
    if kmers == []:
        print("        No kmers found")
        return 0   
    kmers_with_coeffs = zip(kmers, coeffs)
    
    
    print("         Sorting kmers by coefficients")     # Sort kmers by coefficients, greatest to lowest
    sorted_kmers = []
    sorted_kmers_with_coeffs = sorted(kmers_with_coeffs, key=lambda x: x[1])
    for line in sorted_kmers_with_coeffs:
        sorted_kmers.append(line[0])

    if os.path.exists(output_file):
        print("         Output file already exists. Stopping searching")
    else:
        print("     Finding kmers in indexed file. May take a while...")       # Find kmers in indexed reference genome file
        shared_kmers = []
        matching_kmers = []
        with open(indexed_kmer_file, "r") as file:          #open indexed ref genome
            lines = file.readlines()
            for i in range(len(lines)):
                for kmer in sorted_kmers:                   
                    if kmer in lines[i]:                    #check if a phenotype kmer is in this line of indexed ref genome 
                        matching_kmers.append(kmer)
                        index_lines = lines[i:i+int(lines[i].strip().split()[1])+1]         #add all locations of the kmer to the shared kmers list
                        for line in index_lines:
                            shared_kmers.append(line)
                        break


        print("     Number of kmers found in indexed file: ", len(matching_kmers))          # Write phenotype associated kmers found in reference genome to file
        with open(output_file, "w") as file:
            for line in shared_kmers:
                file.write(line)
    print("----- Finished searching for kmers -----")

def filter_id_line(line):

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

def find_genes(pheno_kmers_in_ref_genome_file, filtered_GFF_file, min_mismatches = 0,output_file="genes.txt"):
    print("----- Finding genes -----")

    kmers_with_locations = []
    print("         Extracting kmers and their locations")
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
                kmers_with_locations.append([kmer, *locations])
            else:
                continue
    
    print("         Number of kmers with locations: ", len(kmers_with_locations))


    print("         Identifying kmer parent elements")

    genes = []
    intergenic = []
    used_kmers = []   
    kmers_without_locations = []         

    with open(filtered_GFF_file, "r") as file:          #open filtered GFF file
        lines = file.readlines()
        
        for kmerline in kmers_with_locations:
            if len(kmerline) == 1:
                kmers_without_locations.append(kmerline[0])
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
                                genes.append(add)
                                used_kmers.append(kmerline[0])
                            elif nextLine != "" and int(kmerline[i]) >= int(line[2]) and int(kmerline[i]) <= int(nextLine[1]):
                                intergenic.append([kmerline[0],"intergenic"," ", kmerline[i], nextLine[1]," "])
                                used_kmers.append(kmerline[0])

   

    for item in kmers_with_locations:
        if item[0] not in used_kmers:
            kmers_without_locations.append(item)
            kmers_with_locations.pop(kmers_with_locations.index(item))

    print("         ELEMENTS found: ", len(genes))
    print("         INTERGENIC found: ", len(intergenic))
    print("         UNIDENTIFIED kmers: ", len(kmers_without_locations))


    with open(output_file, "w") as gene_file:
        i = 1
        for line in genes:
            line[5] = filter_id_line(line[5])
            gene_file.write(f"\n------------------------------------------ \n{i}\nKMER: {line[0]}\nELEMENT: {line[1]},START: {line[2]}, LOCUS: {line[3]}, END: {line[4]}\nID: {line[5]} \n------------------------------------------\n")
            i += 1
    with open("intergenic.txt", "w") as intergenic_file:
        i = 1
        for line in intergenic:
            intergenic_file.write(f"\n------------------------------------------\n{i}\nKMER: {line[0]}\nELEMENT: {line[1]},START: {line[2]}, LOCUS: {line[3]}, END: {line[4]}\nID: {line[5]}  \n------------------------------------------\n")
            i += 1
    with open("unidentified_kmers.txt", "w") as unidentified_kmers_file:
        i = 1
        for line in kmers_without_locations:
            unidentified_kmers_file.write(f"\n------------------------------------------\n{i}\nKMER: {line[0]} LOCATIONS: {line[:1]} \n------------------------------------------\n")
            i += 1


# kmerid väljaspool geene ei loeta praegu
# mis ei seostu millegiga

# kmeeride liitmine

# erki artikklist valideerida

# kmeerid asukoha järgi kokku panna
