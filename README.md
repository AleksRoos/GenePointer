### GENEFINDER1

##### Pipeline for automated gene identification from statistically significant kmers in bacterial genomes

GeneFinder1 is an extension program for phenotypeseeker written by Erki Aun. GeneFinder1 makes use of the output of phenotypeseeker to identify genes associated with the resistance mechanism analysed and provides the genes in a list for further analysis in research. The program is at its default settings simple to use and provides results in a matter of minutes depending on the input data amounts. 

#### Requirements
- python3.8 (https://www.python.org/downloads/release/python-380/)
- patric-tools (API for downloading from the PATRIC database) (https://github.com/aldro61/patric_tools)
- phenotypeseeker installed into the same environment (https://github.com/bioinfo-ut/PhenotypeSeeker)
    - GenomeTester4 which comes included with the phenotypeseeker install (https://github.com/bioinfo-ut/GenomeTester4)
- reference genome for the bacteria under investigation (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=2,2157&reference_only=true)
- .gff file for the reference genome that has all the genes with annotations (downloadable along with the reference sequences)

#### Install
- download folder into desired location
- cd GeneFinder1
- pip install .

The module is now in the environment or on your computer globally depending on your setup.

#### Usage
The module has functions that help in downloading data from PATRIC and functions that help to a analyse the input genomes.

Genomes should be in FASTA/FASTQ format if not downloaded using the included scripts.
phenotypeseeker needs a

#### File name descriptions (in order of creation) ""WIP NAMES""
- ##### kmers_and_coefficients_in_(classifier)_model_(antibiotic).txt
    - Created by phenotypeseeker modelling.
    - Contains all k-mers the model considered statistically associated to the phentype.
    - Columns:
        - k-mers
        - coefficient in regression model
        - No. of samples with the k-mer
        - names of samples containing k-mer

- ##### filtered_kmers_and_coeffs.txt
    - Contains only k-mers with coefficients more than 0.
    - Columns:
        - kmer
        - coefficient
- ##### ref_genome_kmer_locs.txt
    - Created by GenomeTester4 when indexing the reference genome.
    - Contains all k-mers (and k-mers with 1 mismatch) in the reference genome with their respective locations.
- ##### pheno_kmers_ref_genome.txt
    - Created when matching significant k-mers (filtered_kmers_and_coeffs.txt) to reference genome k-mers.
    - Contains all k-mers that the reference genome has in common with the significant k-mer list.
- ##### filtered_gff.txt
    - Created by reducing the number of columns and rows in the complete reference genome .GFF file (gene annotation file)
    - Contains genetic elements, their locations and descriptions in filtered form.(only immediately valuable information parsed)
    - Columns:
        - genetic element
        - start base
        - end base
        - information
- ##### genes.txt, intergenic.txt, unidentified.txt
    - Created after locating all kmers on the reference genome and identifying the genes or intergenic regions they belong to. put into unidentified.txt if not anywhere on reference genome.
    - Format:
    ´´´
        "------------------------------------------"
        1
        KMER: AAAAAGGCAGGAC
        ELEMENT: gene,START: 2600731, LOCUS: 2601574, END: 2601879
        ID: ID=gene-Rv2328;Name=PE23;gene=PE23;gene_biotype=protein_coding;gene_biotype=protein_coding 
        "------------------------------------------"
    ´´´