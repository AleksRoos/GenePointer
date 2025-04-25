from gene_finder import data_collection, analyse_genomes
from subprocess import call
#Use full paths everywhere

#check which bacteria and antibiotics data are available
data_collection.check_for("" , min_resist=1000, min_suscept=1000, patric_meta_file_path="/home/sass/Dev/PhenotypeSeeker/GeneFinder1/PATRIC_genomes_AMR.txt")