from gene_pointer import data_collection
import os
import posixpath
from urllib.parse import urlsplit, unquote
from subprocess import call



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
"""

PATRIC_FTP_GENOMES_METADATA_URL = "ftp://ftp.patricbrc.org/RELEASE_NOTES/PATRIC_genomes_AMR.txt"
def download_file_from_url(url, outdir):
    """
    Download a file and save it to some output directory

    Parameters:
    -----------
    url: str
        The URL of the file to download
    outdir: str
        The path to the output directory

    Returns:
    --------
    exception: str
        Empty string if no error, url + exception otherwise

    Notes:
    ------
    * Will automatically skip files that have already been downloaded.

    """
    url = url.strip()
    try:
        from os import system
        print("wget --quiet -o /dev/null -O {0!s} --continue --timeout 20 {1!s}".format(os.path.join(outdir, url_extract_file_name(url)), url))
        system("wget --quiet -o /dev/null -O {0!s} --continue --timeout 20 {1!s}".format(os.path.join(outdir, url_extract_file_name(url)), url))
        return ""
    except Exception as e:
        print(e)
        return url + str(e)

def url_extract_file_name(url):
    """
    Return basename corresponding to url

    """
    urlpath = urlsplit(url).path
    basename = posixpath.basename(unquote(urlpath))
    if os.path.basename(basename) != basename:
        raise ValueError(url)
    return basename

def get_latest_metadata(outdir):
    """
    Downloads the latest genome metadata (not to be confused with AMR metadata)

    Parameters:
    -----------
    outdir: str
        The path to the output directory

    """
    exception = download_file_from_url(PATRIC_FTP_GENOMES_METADATA_URL, outdir)
    if exception != '':
        raise RuntimeError("Failed to download the latest AMR metadata: {0!s}".format(exception))

    return os.path.join(outdir, url_extract_file_name(PATRIC_FTP_GENOMES_METADATA_URL))

def main():

    SPECIES = "mycobacterium tuberculosis" #Species to download genomes for
    ANTIBIOTIC = "ethambutol" #Antibiotic to download genomes for

    inp = input("Are you in the desired directory? Press Enter to continue... n to exit.")
    if inp == "n":
        exit(0)

    inp = input("Do you want to re-download the metafile? Press Enter to continue... y to download.")
    if inp == "y":
        print("Downloading the latest metadata into current folder...")
        get_latest_metadata("./")

    #download genomes and make appropriate data.pheno file in current folder
    dataPhenofolder = data_collection.download(SPECIES, ANTIBIOTIC, 
                    max_genomes=150, #Number of genomes to download
                    patric_meta_file_path="/home/sass/Dev/PhenotypeSeeker/GenePointer/PATRIC_genomes_AMR.txt", #Currently required
                    intermediate_phenos=False,
                    num_threads=8,
                    output_dir="./Genomes") #Whether intermediate phenotype genomes are downloaded

    SPECIES = SPECIES.replace(" ", "_") #Replace spaces with underscores for folder names

    print("DataPheno folder: ", dataPhenofolder)
    print(dataPhenofolder.split("_"))
    DATAPHENO_PATH = dataPhenofolder + "/data_all.pheno" #Path to the data.pheno file. if downloaded genomes is larger than 99 the program makes different data.pheno files for different genome sizes. This is the one that is used in the analysis by default containing all available. example: data_all.pheno, data_100.pheno, data_200.pheno, data_30.pheno contain the respective number of genomes in a equal as possible distribution of R/S phenotypes.
    
    print("DataPheno path: ", DATAPHENO_PATH)

if __name__ == "__main__":
    main()
