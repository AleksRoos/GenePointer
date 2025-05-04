from gene_finder import data_collection


SPECIES = "enterococcus faecium"
ANTIBIOTIC = "vancomycin"
DOWNLOAD_DIRECTORY = "./Genomes"

data_collection.make_data_pheno_file(SPECIES, ANTIBIOTIC, DOWNLOAD_DIRECTORY)