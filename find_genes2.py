from gene_pointer import aligning, analyse_genomes
from subprocess import call
import os
import time
from Bio import SeqIO
from Bio.Seq import Seq

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

def main():

    SPECIES = "enterococcus faecium"
    ANTIBIOTIC = "vancomycin"
    ML_CLASSIFIER = "RF"
    KMER_LENGTH = 13
    DATA_PHENO_PATH = f"./GenomPhenoFiles___{SPECIES.replace(' ', '_')}_{ANTIBIOTIC}/data_all.pheno"

    REFSEQ = bacteria_and_ref_genome_and_gff[SPECIES][0]
    GFF = bacteria_and_ref_genome_and_gff[SPECIES][1]

    print("SPECIES: ", SPECIES)
    print("ANTIBIOTIC: ", ANTIBIOTIC)
    print("ML_CLASSIFIER: ", ML_CLASSIFIER)
    print("KMER_LENGTH: ", KMER_LENGTH)
    print("DATA_PHENO_PATH: ", DATA_PHENO_PATH)
    print("REFSEQ: ", REFSEQ)
    print("GFF: ", GFF)
    print("GENEPOINTER_DIR: ", GENEPOINTER_DIR)

    start = time.time()

    if not os.path.exists(f"kmers_and_coefficients_{ML_CLASSIFIER}_model_{ANTIBIOTIC.capitalize()}.txt"):
        print("Finding significant k-mers...")	
        analyse_genomes.find_significant_kmers(DATA_PHENO_PATH, ML_CLASSIFIER, KMER_LENGTH)
    
    aligning.find_genes_alignment(f"k-mers_and_coefficients_in_{ML_CLASSIFIER}_model_{ANTIBIOTIC.capitalize()}.txt", SPECIES , ANTIBIOTIC, GFF, ref_genome_file = REFSEQ, reduce_kmers_to=10)
    
    print("Time taken for find_gene2.py: ", round((time.time()-start), 6), " s")


if __name__ == "__main__":
    main()