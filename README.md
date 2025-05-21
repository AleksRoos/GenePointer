### GenePointer (WIP)

#### Pipeline for automated genetic region pinpointing from statistically significant kmers in bacterial genomes

GenePointer is an extension program for phenotypeseeker, written by Erki Aun. GenePointer makes use of the output of phenotypeseeker to identify genes associated with the resistance mechanism analysed and provides the genes in a list for further application in research. The program, at its default settings, simple to use and provides results in a matter of minutes depending on the input data size. 

### Requirements
- python3.8 (https://www.python.org/downloads/release/python-380/)
- patric-tools (API for downloading from the PATRIC database) (https://github.com/aldro61/patric_tools)
- PhenotypeSeeker installed into the same environment (https://github.com/bioinfo-ut/PhenotypeSeeker)
    - GenomeTester4 which comes included with the phenotypeseeker install (https://github.com/bioinfo-ut/GenomeTester4)
- reference genome for the bacteria under investigation (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2,2157&reference_only=true)
- .gff file for the reference genome that has all the genes with annotations (downloadable along with the reference sequences)
- bowtie2

#### Install
- download folder into desired location
- cd GenePointer
- pip install .

The module is now in the environment or on your computer globally depending on whther you installed it in an environment or outside.

#### Usage
The module has functions that help in downloading data from PATRIC and functions to analyse the input genomes.

* Data
Genomes should be in FASTA/FASTQ format and labelled with phenotypes if not downloaded using the included script ("download_genomes.py").
PhenotypeSeeker needs a file that contains, line by line, each input genome with lines in tab separated format: SampleID(GenePointer takes this as the id of the genome downloaded form PATRIC but it can be anything), Address (preferably the complete path to the genome), Antibiotic_Name (the binary(1/0)) value representing resistance or susceptibility.

all scripts in the GenePointer folder should be run in the folder that analysis results end up in. Its best to use species and antibiotic names for folders. Example: run python scripts in /SPECIES/ANTIBIOTIC/{results end up here}

download_genomes.py creates the folder Genomes for downloaded genomes and GenomPhenoFiles___{species_name}_{antibiotic} for all the data_*.pheno files.
Sometimes genome files downloaded from PATRIC can be empty and the program will run into errors with these included in the data.pheno files. A bash script that takes the Genome directory and data.pheno path as input can fix this by checking if the file is empty and then removing it from the corresponding line in data.pheno.

find_genes.py and find_gene2.py are analysis scripts and should be run in the folder where analysis results should go. They can take time and so it's best to run them using no hangup.

If you already have the genomes with their corresponding prefixes ("1_" for positive phenotype, "0_" for negative phenotype) just run the make_genome_input_file.py script with the name of the species and the name of the antibiotic to make the data_pheno file that is used as input to PhenotypeSeeker. Intermediate values are not currently processable by GenePointer.

summarise.py should be run after the find_genes*.py scripts have finished. It requires the output files: genes.csv, intergenic.csv and chi2_results_* from find_genes.py and alignments.csv, no_alignments.csv and chi2_results from find_genes2.py.  

Once you have the genomes downloaded and their corresponding ID,Address,Phenotype table file. Pass the table file's absolute path to the find_genes.py script where the DATAPHENO_PATH is defined and specify the species you are analysing to the SPECIES variable. Upon running the script everything should happen without further input.

The results of analysis will be written to the same folder the script is run in.

#### File name descriptions (in order of creation) ""WIP NAMES""
- ##### kmers_and_coefficients_in_(classifier)_model_(antibiotic).txt
    - Created by phenotypeseeker modelling.
    - Contains all k-mers the model considered statistically associated to the phentype.
    - Columns:
        - k-mers
        - coefficient in regression model
        - No. of samples with the k-mer
        - names of samples containing k-mer

- ##### ref_genome_kmer_locs.txt
    - Created by GenomeTester4 when indexing the reference genome.
    - Contains all k-mers (and k-mers with 1 mismatch) in the reference genome with their respective locations.

- ##### pheno_kmers_ref_genome.txt
    - Created when matching significant k-mers (filtered_kmers_and_coeffs.txt) to reference genome k-mers.
    - Contains all k-mers that the reference genome has in common with the significant k-mer list.

- ##### genes.csv, intergenic.csv
    - Created after locating all kmers on the reference genome and identifying the genes or intergenic regions they belong to. put into unidentified.txt if not anywhere on reference genome.
    
    - genes.csv Columns: k-mer, element class, element start, k-mer location, element end, element annotation
    - intergenic.csv Columns: k-mer, element class, element start, k-mer location, element end

- ##### alignments.csv, no_alignments,csv
    - Created by find_genes2.py and contains the results of aligning all the extended k-mers to the reference genome.

    - alignments.csv columns: kmer_id, coefficient of k-mer in output ML model, genome id, reference genome id, extended k-mer middle base position, sequence, strand, gene annotations
    - no_alignments.csv columns: kmer_id, coefficient of k-mer in output ML model, genome id, reference genome id, extended k-mer middle base position, sequence, strand(+,-)

- ##### gene_summary_align.csv, unaligned_summary.csv
    - Created after analysing k-mers using alignment to the reference genome in find_genes2.py

    - gene_summary_align.csv columns: Gene (total number of genes found), Minimum of k-mer p-values, Signif. kmer prevalence (Tot/Res/Sus), Num of Unique k-mers (total k-mers aligned), Gene prevalence in genomes (Tot/Res/Sus), Unique k-mers in gene, Signif. k-mer location, Genomes with gene.
    - unaligned_summary.csv columns: K-mer, k-mer P-value, k-mer Prevalence, Extended k-mer Sequence

- ##### gene_summary_noalign.csv, intergenic_summary.csv
    - Created after placing k-mers using direct mapping in the find_genes.py script.

    - gene_summary_noalign.csv columns: Genes (total number of genes found), Signif. k-mer p-value, Signif. k-mer prevalence, num of unique k-mers, unique k-mers, Locations
    - intergenic_summary.csv columns: intergenic region id(inter_Start...End), Signif. k-mer p-value, Signif. k-mer prevalence, K-mer, Location