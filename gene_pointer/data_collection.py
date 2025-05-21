"""
patric_tools: A Python package to download data from the PATRIC database
Copyright (C) 2017 Alexandre Drouin
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Example: downloading AMR data sets from the PATRIC database
"""

from joblib import Parallel, delayed
from patric_tools import amr, genomes
import os
from subprocess import call
import random

# if not os.path.exists("/tmp/PATRIC_genomes_AMR.txt"):
#     call(["wget", "ftp://ftp.patricbrc.org/patric2/current_release/RELEASE_NOTES/PATRIC_genomes_AMR.txt", "-O", "/tmp/PATRIC_genomes_AMR.txt"])


MYCOBAC_ISONIA = "/mycobacterium_tuberculosis__isoniazid"
ENTEROC_VANCO = "/enterococcus_faecalis__vancomycin"
PSEUDOM_CIPROFLO = "/pseudomonas_aeruginosa__ciprofloxacin"

species_folders = {"mycobacterium tuberculosis":MYCOBAC_ISONIA, "enterococcus faecalis":ENTEROC_VANCO, "pseudomonas aeruginosa": PSEUDOM_CIPROFLO}
test_species = ["mycobacterium tuberculosis", "enterococcus faecalis", "pseudomonas aeruginosa"]
test_antibiotics = ["Isoniazid", "Vancomycin", "Ciprofloxacin"]
data_pheno_lines = []
data_pheno_lines_susceptible = []
data_pheno_lines_resistant = []
ids_labels = []

THROTTLE = True


#Isoniazid: Mycobacterium tuberculosis
#Vancomycin: Enterococcus faecalis või Enterococcus faecium
#Ciprofloxacin: Pseudomonas aeruginosa (hiljem ka Klebsiella, E.coli)

# Find all antibiotic resistance data sets that match some criteria

#Check if there is data available
#id, species, pheno = amr.get_amr_data_by_species_and_antibiotic(antibiotic=test_antibiotics[0], species=test_species, drop_intermediate=True, amr_metadata_file="/tmp/PATRIC_genomes_AMR.txt")

def choose_random_items(items, x):
    if x > len(items):
        print(f"        Requested {x} items, but only {len(items)} available. Returning all items.")
        return items[:]
    return random.sample(items, x)

def download_genome(g_id, GENOME_FOLDER):
    try:
        genomes.download_genome_contigs(g_id[0], label=str(g_id[1]), outdir=GENOME_FOLDER, throttle=THROTTLE, throttle_time=.5)
    except Exception as e:
        print("         Error downloading genome: ", g_id[0], e)


def make_data_pheno_file(species, antibiotic, DOWNLOAD_DIRECTORY):
    '''
        Makes input file for PhenotypeSeeker in GenePointer out of genomes in specified folder (DOWNLOAD_DIRECTORY). 
        Useful when you already have genomes downloaded. 
        !! genome file names must be prefixed with phenotype value (1 for resistant, 0 for susceptible).
        Example genome file name(resistant genome): 1_genomeID.fna

        Parameters:
            species: str name of species
            antibiotic: str name of antibiotic
            DOWNLOAD_DIRECTORY: path relative path to folder with downloaded genomes (format as specified before) 
        Returns:
            dir_name: str directory with data.pheno file/s
    '''
    
    print(f"        Making data.pheno file for {species} {antibiotic}")
    folder = (species.replace(" ", "_") + "__" + antibiotic).lower()
    

    print("-------Making data.pheno file into ", folder, "-------")
    global data_pheno_lines_resistant, data_pheno_lines_susceptible
    filenames = [f for f in os.listdir(os.path.join(DOWNLOAD_DIRECTORY, folder))]
    
    for filename in filenames:
        address = os.path.abspath(os.path.join(DOWNLOAD_DIRECTORY,folder,filename))
        if filename[0] == "1":
            data_pheno_lines_resistant.append(str(filename[2:-4] + "\t" + address + "\t" + filename[0] + "\n"))
        elif filename[0] == "0":
            data_pheno_lines_susceptible.append(str(filename[2:-4] + "\t" + address + "\t" + filename[0] + "\n"))    
        else:
            continue
    number_of_genomes = len(data_pheno_lines_resistant) + len(data_pheno_lines_susceptible)
    print(f"        Number of genomes: {number_of_genomes}, Num resistant: {len(data_pheno_lines_resistant)}, Num susceptible: {len(data_pheno_lines_susceptible)}")


    species = species.replace(" ", "_")
    antibiotic = antibiotic.replace(" ", "_")
    dir_name = "./GenomPhenoFiles___" + species + "_" + antibiotic
    if not os.path.exists(dir_name):
        call(["mkdir " + dir_name], shell=True)

    
    print("         Writing data.pheno files")
    
    def write_pheno_file(filename, num_resistant, num_susceptible):
        with open(os.path.join(dir_name,filename), "w") as file:
            lines = choose_random_items(data_pheno_lines_resistant, num_resistant)
            if len(lines) < num_resistant:
                num_susceptible += num_resistant - len(lines)
            elif len(lines) > num_resistant:
                num_susceptible -= len(lines) - num_resistant
            lines.extend(choose_random_items(data_pheno_lines_susceptible, num_susceptible))
            file.write(f"SampleID\tAddress\t{antibiotic.capitalize()}\n")
            for line in lines:
                file.write(line)
    
    if number_of_genomes > 2999:
        write_pheno_file("data_3000.pheno", 1500, 1500)
        write_pheno_file("data_2000.pheno", 1000, 1000)
        write_pheno_file("data_1600.pheno", 800, 800)
    if number_of_genomes > 1000:
        write_pheno_file("data_1000.pheno", 500, 500)
        write_pheno_file("data_350.pheno", 175, 175)
        write_pheno_file("data_200.pheno", 100, 100)
        write_pheno_file("data_500.pheno", 250, 250)
    if number_of_genomes > 100:
        write_pheno_file("data_100.pheno", 50, 50)
        write_pheno_file("data_30.pheno", 15, 15)
        write_pheno_file("data_10.pheno", 5, 5)

    with open(os.path.join(dir_name, "data_all.pheno"), "w") as file:

        print("         Resistant: ", len(data_pheno_lines_resistant))
        print("         Susceptible: ", len(data_pheno_lines_susceptible))
        file.write(f"SampleID\tAddress\t{antibiotic.capitalize()}\n")
        for line in data_pheno_lines_resistant:
            file.write(line)
        for line in data_pheno_lines_susceptible:
            file.write(line)
        
    #os.chdir("/home/sass/Dev/MyEnv/patric_tools-master/examples")
        
    data_pheno_lines_resistant = []
    data_pheno_lines_susceptible = []
    file.close()
    return dir_name

def check_for(_species, min_resist=20, min_suscept=20, patric_meta_file_path="/tmp/PATRIC_genomes_AMR.txt"):
    print("------Checking for available data sets for species: ", _species, "------")

    amr_datasets = amr.list_amr_datasets(min_resistant=min_resist,
                                     min_susceptible=min_suscept,
                                     single_species=True, amr_metadata_file=patric_meta_file_path)
    print ("        Available antibiotic resistance data sets: for species ", _species)

    for species, antibiotic in amr_datasets:
        if species[0].lower() == _species.lower():
            print("         Species:", species[0].title(), "  Antibiotic:", antibiotic)
        elif _species == "":
            print("         Species:", species[0].title(), "  Antibiotic:", antibiotic)

    print("\n" * 2)

#take one species and antibiotic, download data and make data.pheno file for it
def download(_species, _antibiotic, patric_meta_file_path="/tmp/PATRIC_genomes_AMR.txt" , max_genomes=100, intermediate_phenos=False, num_threads = 4, output_dir="./Genomes/"):
    '''
        Downloads genome files from patric, even number of resistant and susceptible, prefixes them with phenotype values and makes data.pheno file for PhenotypeSeeker
        Only for binary phenotypes currently
        Parameters:
            _species: str species name
            _antibiotic: str antibiotic name
            patric_meta_file_path: path path to PATRIC_genomes_AMR.txt 
                (downloadable from ftp://ftp.patricbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt)
            max_genomes: int max number of genomes to download (includes both phenotypes)
            intermediate_phenos: bool if intermediate phenotype values should be downloaded (requires additional processing)
            num_threads: int threads to use when downloading in parallel
            output_dir: path folder where genome files will be downloaded
        Returns:
            0 if no genomes could be downloaded
            data.pheno folder path if successful
    '''

    _species = _species.lower()
    _antibiotic = _antibiotic.lower()
    
    print("------Downloading genomes for species: ", _species, " and antibiotic: ", _antibiotic, "------")
    DOWNLOAD_DIRECTORY = output_dir
    GENOME_FOLDER = os.path.join(DOWNLOAD_DIRECTORY, _species.replace(" ", "_") + "__" + _antibiotic.lower())
    _, genome_ids, labels = \
        amr.get_amr_data_by_species_and_antibiotic(species=[_species.lower()],
                                                antibiotic=_antibiotic.lower(),
                                                drop_intermediate=(not intermediate_phenos), amr_metadata_file=patric_meta_file_path)

    resistant = []
    susceptible = []
    #take even amount of resistant and susceptible genomes
    for id, phenotype in zip(genome_ids, labels):
        ids_labels.append((id, phenotype))
        if phenotype == 1 and len(resistant) < int(max_genomes / 2):
            resistant.append((id, phenotype))
        elif phenotype == 0 and len(susceptible) < int(max_genomes / 2):
            susceptible.append((id, phenotype))
    even_ids_labels = resistant + susceptible

    print("     Available Genomes: {}".format(len(labels)))
    print("         Resistant: {}".format((labels == 1).sum()))
    print("         Sensible: {}".format((labels == 0).sum()))

    if len(labels) == 0:
        print("     NO GENOMES FOUND")
        return 0
    
    files = []
    for item in ids_labels:
        files.append(str(item[1])+"_"+str(item[0]) + ".fna")

    #(all(x in files for x in os.listdir(GENOME_FOLDER)) and  TEMP REMOVAL FROM BELOW CONDITIONAL
    if os.path.exists(GENOME_FOLDER) and (len(os.listdir(GENOME_FOLDER)) == max_genomes):
        inp = input("     Directory already exists and includes genomes. Re-Download? -- y/n")
        if inp == "y":
            print("     Outputting genomes to: ", GENOME_FOLDER)
            print("     DOWNLOADING....")
            Parallel(n_jobs=num_threads)(delayed(download_genome)(g_id, GENOME_FOLDER) for g_id in even_ids_labels)
            print("     Done downloading genomes")
    else:
        print("     Fetching genome sequences:")
        if not os.path.exists(GENOME_FOLDER):
            print("     Creating directory: ", GENOME_FOLDER)
            call(["mkdir -p " + GENOME_FOLDER], shell=True)
        print("     Outputting genomes to: ", GENOME_FOLDER)
        print("     DOWNLOADING....")
        Parallel(n_jobs=num_threads)(delayed(download_genome)(g_id, GENOME_FOLDER) for g_id in even_ids_labels)
        print("     Done downloading genomes")

    
    return make_data_pheno_file(_species, _antibiotic, DOWNLOAD_DIRECTORY)