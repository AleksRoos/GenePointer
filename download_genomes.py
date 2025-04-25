from gene_finder import data_collection


SPECIES = "mycobacterium tuberculosis" #Species to download genomes for
ANTIBIOTIC = "isoniazid" #Antibiotic to download genomes for

inp = input("Are you in the desired directory? Press Enter to continue... n to exit.")
if inp == "n":
    exit(0)

#download genomes and make appropriate data.pheno file in current folder
dataPhenofolder = data_collection.download(SPECIES, ANTIBIOTIC, 
                   max_genomes=150, #Number of genomes to download
                   patric_meta_file_path="/home/sass/Dev/PhenotypeSeeker/GeneFinder1/PATRIC_genomes_AMR.txt", #Currently required
                   intermediate_phenos=False,
                   num_threads=8,

                   output_dir="./Genomes") #Whether intermediate phenotype genomes are downloaded

#if genomes already downloaded, use this line instead of get_data.download
#dataPhenofolder = get_data.make_data_pheno_file(_species.replace(" ", "_") + "__" + _antibiotic, "isoniazid", min_resist=20, min_suscept=20, max_genomes=50, patric_meta_file_path="/home/sass/Dev/PhenotypeSeeker/GeneFinder1/PATRIC_genomes_AMR.txt", intermediate_phenos=False)

print("DataPheno folder: ", dataPhenofolder)
print(dataPhenofolder.split("_"))

DATAPHENO_PATH = dataPhenofolder + "/data_all.pheno" #Path to the data.pheno file. if downloaded genomes is larger than 99 the program makes different data.pheno files for different genome sizes. This is the one that is used in the analysis by default containing all available. example: data_all.pheno, data_100.pheno, data_200.pheno, data_30.pheno contain the respective number of genomes in a equal as possible distribution of R/S phenotypes.

print("DataPheno path: ", DATAPHENO_PATH)


