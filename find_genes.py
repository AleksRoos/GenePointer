from GenePointer.gene_pointer import other_functions
from gene_pointer import data_collection
from subprocess import call
import os
import time

GENEPOINTER_DIR = "/home/sass/Dev/PhenotypeSeeker"

bacteria_and_ref_genome_and_gff = {
    "mycobacterium tuberculosis": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/MycoTuber/GCA_000195955.2_ASM19595v2_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/MycoTuber/genomic.gff"],
    "klebsiella pneumoniae": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/KlebsPneum/GCF_023546055.1_ASM2354605v1_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/KlebsPneum/genomic.gff"],
    "escherichia coli": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/EColi/GCA_000005845.2_ASM584v2_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/EColi/genomic.gff"],
    "pseudomonas aeruginosa": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/PseudoAerugi/GCA_000006765.1_ASM676v1_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/PseudoAerugi/genomic.gff"],
    "enterococcus faecium": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/EnteroFaecium/GCF_003071425.1_ASM307142v1_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/EnteroFaecium/genomic.gff"],
    "staphylococcus aureus": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/StaphylAureus/GCA_000013425.1_ASM1342v1_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/StaphylAureus/genomic.gff"],
    "streptococcus pneumoniae": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/StreptoPneum/GCA_001457635.1_NCTC7465_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/StreptoPneum/genomic.gff"],
    "salmonella enterica": [GENEPOINTER_DIR + "/GenePointer/RefGenomes/SalmEnter/GCA_000006945.2_ASM694v2_genomic.fna", GENEPOINTER_DIR + "/GenePointer/RefGenomes/SalmEnter/genomic.gff"],
}

SPECIES = "mycobacterium tuberculosis" #Species name
ANTIBIOTIC = "isoniazid" #Antibiotic name
DATAPHENO_PATH = f"/home/sass/Dev/PhenotypeSeeker/MycobacTuberGenes/Isoniazid/GenomPhenoFiles___{SPECIES.replace(" ", "_")}_{ANTIBIOTIC}/data_all.pheno"

#For PhenotypeSeeker
#RandomForest - RF, logistic - log_reg  (choose from 'log', 'SVM', 'RF', 'NB', 'XGBC', 'DT')
REGRESSION_MODEL = "RF"
KMER_LENGTH = 13

#likely unnecessary
MIN_MISMATCHES = 0

def main():
    print("------ GenePointer -------")
    start = time.time()
    global SPECIES, ANTIBIOTIC, REF_GENOME_PATH, GFF_PATH, DATAPHENO_PATH, KMER_LENGTH, MIN_MISMATCHES,REGRESSION_MODEL
    inp = input("       Are you in the desired directory? Press Enter to continue... n to exit.")
    if inp == "n":
        exit(0)    

    if SPECIES == "":
        print("      No species defined!")
        exit()
    else:
        if SPECIES not in list(bacteria_and_ref_genome_and_gff.keys()):
            print("         Species name not correct or no reference genomes for this species!")
            exit()
    if ANTIBIOTIC == "":
        print("      No antibiotic defined!")
        exit()
    



    SPECIES = SPECIES.lower()
    ANTIBIOTIC = ANTIBIOTIC.lower() 

    #Use full paths everywhere
    REF_GENOME_PATH = ""
    GFF_PATH = ""
    if bacteria_and_ref_genome_and_gff != {}:
        REF_GENOME_PATH = bacteria_and_ref_genome_and_gff[SPECIES][0] #Path to the reference genome
        GFF_PATH = bacteria_and_ref_genome_and_gff[SPECIES][1] #Path to the GFF file

    if DATAPHENO_PATH == "":
        print("         no ID,Address,Phenotype file path defined")
        # print("         looking for data.pheno file")
        # print("GenomePhenoFiles___" + SPECIES.replace(" ", "_") + "_" + ANTIBIOTIC)
        # if  os.path.exists("./GenomePhenoFiles___" + SPECIES.replace(" ", "_") + "_" + ANTIBIOTIC):
        #     print(os.listdir("GenomePhenoFiles___" + SPECIES.replace(" ", "_") + "_" + ANTIBIOTIC))
        #     DATAPHENO_PATH = "GenomePhenoFiles___" + SPECIES.replace(" ", "_") + "_" + ANTIBIOTIC + "/data_all.pheno" #Path to the data.pheno file. if downloaded genomes is larger than 99 the program makes different data.pheno files for different genome sizes. This is the one that is used in the analysis by default containing all available. example: data_all.pheno, data_100.pheno, data_200.pheno, data_30.pheno contain the respective number of genomes in a equal as possible distribution of R/S phenotypes.
        exit()
  
    

    print("------ Running GenePointer -------")


    #----------------------FIND SIGNIFICANT KMERS----------------------
    #required folder for phenotypeseeker
    if not os.path.exists("./K-mer_lists"):
        call(["mkdir K-mer_lists"], shell=True)
    good_result = other_functions.find_significant_kmers(DATAPHENO_PATH, classifier = REGRESSION_MODEL, kmer_length = 13)
    if os.path.exists("K-mer_lists/"):
        print("------Emptying K-mer_lists directory-------")
        call(["rm -rf K-mer_lists/*"], shell=True)
    if not good_result:
        print("        No significant kmers found")
        print("-----------------ENDING ANALYSIS------------------")
        exit(0)



    #---------------- REF GENOME AND ANNOTATED GENE FILE PREPROCESSING
    #find all kmers in ref genome and their locations -> file = ref_genome_kmer_locs.txt
    other_functions.index_genome(REF_GENOME_PATH, output_file="ref_gnome_kmer_locs.txt", kmer_length=KMER_LENGTH)
    #make GFF file smaller, take only required columns
    other_functions.filter_GFF_file(GFF_PATH, output_file="filtered_gff.txt")

    #-----------------------FIND GENES IN REF GENOME
    #find kmers in common with the reference genome and the significant kmers -> file = pheno_kmers_ref_genome.txt

    if os.path.exists("./pheno_kmers_ref_genome.txt"):
        inp = input("           Do you want to look for kmers in reference index again? n to stop")
        if inp != "n":
            other_functions.find_ref_genome_kmers("ref_gnome_kmer_locs.txt", "significant_lines_from_kmer_coeff_file.txt",output_file="pheno_kmers_ref_genome.txt")
    #find genes
    other_functions.find_genes("pheno_kmers_ref_genome.txt", "filtered_gff.txt")

    

    ##CLEARING FILES I DONT NEED RIGHT NOW

    # call(["rm -rf ref_gnome_kmer_locs.txt"], shell=True)
    # call(["rm -rf filtered_kmers_and_coeffs.txt"], shell=True)
    # call(["rm -rf filtered_gff.txt"], shell=True)
    # call(["rm -rf pheno_kmers_ref_genome.txt"], shell=True)
    # call(["rm -rf genes.txt"], shell=True)
    # call(["rm -rf K-mer_lists/*"], shell=True)
    # call(["rm -rf Genomes/*"], shell=True)
    # call(["rm -rf Genomes"], shell=True)
    # call(["rm -rf data_all.pheno"], shell=True)


    print("DEV ------------REMOVING UNWANTED FILES-----------------")
    call(["rm -rf distances.mat"], shell=True)
    call(["rm -rf *.list"], shell=True)
    call(["rm -rf log.txt"], shell=True)
    call(["rm -rf tree_xml.txt"], shell=True)
    call(["rm -rf tree_newick.txt"], shell=True)
    call(["rm -rf mash_distances.mat"], shell=True)
    call(["rm -rf out_13.index"], shell=True)
    call(["rm -rf reference.msh"], shell=True)
    call(["rm -rf ref_gnome_kmer_locs"], shell=True)
    call(["rm -rf kmers_genomes.txt"], shell=True)



    #implement way of making this program command line firendly
    #give species and antibitoic on cmd line
    #how to automate data.pheno file usage so I dont have to give the path to it.
    #fix the issues with file paths
    #       the paths should be relative to the current working directory but they are absolute for now.
    print("-------------------find genes DONE------------------------\n")
    print("Time taken for find_genes.py program: ", round((time.time()-start) / 60, 3), " min")


if __name__ == "__main__":
    main()
