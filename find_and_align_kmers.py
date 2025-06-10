from GenePointer.gene_pointer import analysis, other_functions
from gene_pointer import summarise_data
from subprocess import call
import os
import time
from Bio import SeqIO
from Bio.Seq import Seq
from dotenv import get_key




def main():

    try:
        GENEPOINTER_DIR = get_key("/home/sass/Dev/PhenotypeSeeker/.env","GENEPOINTER_DIR")
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
        #FOR PHENOTYPESEEKER
        ML_CLASSIFIER_DICT = {
            "log": "log_reg",
            "log_reg": "log_reg",
            "logistic_regression": "log_reg",
            "logistic regression": "log_reg",
            "RF": "RF",
            "random_forest": "RF",
            "random forest": "RF",
            "RandomForest": "RF",
        }

        print(">>>>>>>>>>>>>>>>> RUNNING GENEPOINTER <<<<<<<<<<<<<<<<<")


        analysis.name_files_with_seqID("/home/sass/Dev/PhenotypeSeeker/GenePointer/RefGenomes/EnteroFaecium")
        
        #INPUT PARAMETERS
        SPECIES = "enterococcus faecium" #Species to download genomes for
        ANTIBIOTIC = "vancomycin"
        ML_CLASSIFIER = "log"
        POPULATION_CORRECTION = True

        KMER_LENGTH = 13
        DATA_PHENO_PATH = f"./GenomPhenoFiles___{SPECIES.replace(' ', '_')}_{ANTIBIOTIC}/data_all.pheno"

        REFSEQ = bacteria_and_ref_genome_and_gff[SPECIES][0]
        GFF = bacteria_and_ref_genome_and_gff[SPECIES][1]

        #summarise_data.summarise_all()

        print("    CHECK PARAMETERS    (s to skip all, enter to continue, or type to change)")
        inp = input(f"    GENEPOINTER_DIR: {GENEPOINTER_DIR} -> ")
        if inp.lower() != "s":
            GENEPOINTER_DIR = GENEPOINTER_DIR if inp == "" else inp
            inp = input(f"    SPECIES: {SPECIES} -> ")
            SPECIES = SPECIES if inp == "" else inp.lower()
            inp = input(f"    ANTIBIOTIC: {ANTIBIOTIC} -> ")
            ANTIBIOTIC = ANTIBIOTIC if inp == "" else inp.lower()
            inp = input(f"    ML_CLASSIFIER: {ML_CLASSIFIER} -> ")
            ML_CLASSIFIER = ML_CLASSIFIER if inp == "" else inp.lower()
            inp = input(f"    KMER_LENGTH: {KMER_LENGTH} -> ")
            KMER_LENGTH = KMER_LENGTH if inp == "" else int(inp)
            inp = input(f"    DATA_PHENO_PATH: {DATA_PHENO_PATH} -> ")
            DATA_PHENO_PATH = DATA_PHENO_PATH if inp == "" else inp
            # inp = input(f"    REFSEQ: {REFSEQ} -> ")
            # REFSEQ = REFSEQ if inp == "" else inp
            # inp = input(f"    GFF: {GFF} -> ")
            # GFF = GFF if inp == "" else inp



        print("    PARAMETERS:")
        print("        GENEPOINTER_DIR: ", GENEPOINTER_DIR)
        print("        SPECIES: ", SPECIES)
        print("        ANTIBIOTIC: ", ANTIBIOTIC)
        print("        ML_CLASSIFIER: ", ML_CLASSIFIER)
        print("        KMER_LENGTH: ", KMER_LENGTH)
        print("        DATA_PHENO_PATH: ", DATA_PHENO_PATH)
        print("        REFSEQ: ", REFSEQ)
        print("        GFF: ", GFF)

        inp = input("    Are these parameters correct? enter to continue: ")
        if inp.lower() != "":
            print("    Exiting the analysis. Please run the script again with corrected parameters.")
            return 1
        else:
            print("    Parameters confirmed. Proceeding with the analysis.")

        start = time.time()
        print(f"    Checking for: k-mers_and_coefficients_in_{ML_CLASSIFIER_DICT[ML_CLASSIFIER]}_model_{ANTIBIOTIC.capitalize()}.txt")
        if not os.path.exists(f"k-mers_and_coefficients_in_{ML_CLASSIFIER_DICT[ML_CLASSIFIER]}_model_{ANTIBIOTIC.capitalize()}.txt"):
            if not os.path.exists("./K-mer_lists"):
                call(["mkdir K-mer_lists"], shell=True)
            good_result = analysis.find_significant_kmers(DATA_PHENO_PATH, ML_CLASSIFIER, KMER_LENGTH, POPULATION_CORRECTION)
            if os.path.exists("K-mer_lists/"):
                print("    Emptying K-mer_lists directory")
                call(["rm -rf K-mer_lists/*"], shell=True)
            if not good_result:
                print("    No significant kmers found")
                print("------------------ENDED ANALYSIS------------------")
                exit(0)
        else:
            inp = input("    Find k-mers again? n to skip .... ")
            if inp != "n":
                if not os.path.exists("./K-mer_lists"):
                    call(["mkdir K-mer_lists"], shell=True)
                good_result = analysis.find_significant_kmers(DATA_PHENO_PATH, ML_CLASSIFIER, KMER_LENGTH, POPULATION_CORRECTION)
                if os.path.exists("K-mer_lists/"):
                    print("    Emptying K-mer_lists directory")
                    call(["rm -rf K-mer_lists/*"], shell=True)
                if not good_result:
                    print("    No significant kmers found")
                    print("------------------ENDED ANALYSIS------------------")
                    exit(0)
                
        
        REFSEQ_DIR = "/home/sass/Dev/PhenotypeSeeker/GenePointer/RefGenomes/EnteroFaecium"

        analysis.find_genes_alignment(f"k-mers_and_coefficients_in_{ML_CLASSIFIER_DICT[ML_CLASSIFIER]}_model_{ANTIBIOTIC.capitalize()}.txt", SPECIES , ANTIBIOTIC, ref_genome_dir = REFSEQ_DIR, reduce_genomes_to=10)
            
        print("    Time taken for k-mer alignment-based marker search: ", round((time.time()-start), 6), " s")

        print("")

        

        summarise_data.summarise_all()


        try:
            call(["mkdir Results"], shell=True)
            #call(["mv k-mers_and_coefficients_in_" + ML_CLASSIFIER_DICT[ML_CLASSIFIER] + "_model_" + ANTIBIOTIC.capitalize() + ".txt Results/"], shell=True)
            # call(["mv Aligned_kmer_results.csv Results/"], shell=True)
            # call(["mv Un_aligned_kmer_results.csv Results/"], shell=True)
            call(["mv Summary_aligned_kmers.csv Results/"], shell=True)
            call(["mv Summary_unaligned_kmers.csv Results/"], shell=True)
            call(["cp kmers_genomes_sequences_table.csv Results/"], shell=True)
        except FileExistsError:
            print("    Results directory already exists, skipping creation.")
        except Exception as e:
            print("    Error while moving to /Results directory:", str(e))
        

        print("------------------ENDED ANALYSIS------------------")
        return 1
    except KeyboardInterrupt:
        print("\nAnalysis interrupted by user.\n")
        return 1
    except Exception as e:
        print("\nERROR OCCURED DURING ANALYSIS", str(e), "\n")
        return 1
if __name__ == "__main__":
    main()