from gene_finder import analyse_genomes
from subprocess import call
import os


analyse_genomes.find_genes_alignment("significant_lines_from_kmer_coeff_file.txt", "mycobacterium tuberculosis", "moxifloxacin", ref_genome_file = "/home/sass/Dev/PhenotypeSeeker/GeneFinder1/RefGenomes/MycoTuber/GCA_000195955.2_ASM19595v2_genomic.fna")

#analyse_genomes.find_kmers_in_genome("TGGCCGCGCGCAA", "/home/sass/Dev/PhenotypeSeeker/MycobacTuberGenes/Moxifloxacin/Genomes/mycobacterium_tuberculosis__moxifloxacin/0_1733.1219.fna")