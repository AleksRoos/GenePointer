import pandas as pd
from collections import Counter, defaultdict
import re

def summaries():
    
    locations_dataframe = pd.read_csv("locations_matrix.csv", sep=",")
    colheads = locations_dataframe.columns[1:]
    rowheads = locations_dataframe.index
    number_of_genomes = len(rowheads)
    number_of_resistant = sum(1 for x in rowheads if x[0] == "1")
    unique_locations = []
    num_unique_locations = []
    for i in range(len(colheads)):
        locations = []
        for item in locations_dataframe.iloc[:, i]:
            try:
                locations.append(item.strip().split("|"))
            except:
                locations.append(item.strip())
        
        locations = sorted(set([item.strip() for sublist in locations for item in sublist if item not in ["NoGenes", ""]]), key=lambda x: int(x))
        num_unique_locations.append(len(locations))
        unique_locations.append(locations)

    genes_dataframe = pd.read_csv("gene_matrix.csv", sep=",")
    unique_genes = []
    num_unique_genes = []
    for i in range(len(colheads)):
        genes = []
        for item in genes_dataframe.iloc[:, i]:
            try:
                genes.append(item.strip().split("|"))
            except:
                genes.append(item.strip())
        genes = set([item.strip() for sublist in genes for item in sublist if item not in ["(NoGenes)", ")"]])
        num_unique_genes.append(len(genes))
        unique_genes.append(list(genes))
    
    alignments_dataframe = pd.read_csv("alignment_matrix.csv", sep=",")
    unique_alignments = []
    num_unique_alignments = []

    for i in range(len(colheads)):
        alignments = []
        for item in alignments_dataframe.iloc[:,i]:
            try:
                alignments.append(item.strip().split("|"))
            except:
                alignments.append(item.strip())
        alignments = set([item.strip() for sublist in alignments for item in sublist if item not in ["0;-----NotPresent-----"]])
        unique_alignments.append(alignments)
        num_unique_alignments.append(len(alignments))

    with  open("gene_finder_summary.csv", "w") as summary_file:
        summary_file.write(f"Total genomes analysed {number_of_genomes} of  which resistant {number_of_resistant}\n")
        summary_file.write(f"K-mer and reverse complement,Coefficient in ML model,Prevalence in input genomes,Num of unique locations,Unique locations,Num of unique genes,Unique genes or locations in genes,Num of unique alignments,Unique alignments and their locations\n")
        for i in range(len(colheads)):
            summary_file.write(colheads[i].strip().split(" ")[0] + ",")
            summary_file.write(colheads[i].strip().split(" ")[4] + ",")
            summary_file.write(colheads[i].strip().split(" ")[1] + colheads[i].strip().split(" ")[2] + ",")
            summary_file.write(f"{num_unique_locations[i]},")
            summary_file.write(f"{' '.join(unique_locations[i])},")
            summary_file.write(f"{num_unique_genes[i]},")
            summary_file.write(f"{'--------'.join(unique_genes[i])},")
            summary_file.write(f"{num_unique_alignments[i]},")
            summary_file.write(f"{'    '.join(unique_alignments[i])},")
            summary_file.write("\n")
    ###

    return 1

def summarise():

    kmers2, pvalues = readPvalue("chi2_results_Isoniazid.tsv")

    print(len(pvalues))

    descriptions = []
    kmers = []

    with open("genes.csv", "r") as genes_file:
        lines = genes_file.readlines()

        for line in lines:
            line = line.strip().split(",")
            kmers.append(line[0].strip())
            descriptions.append(line[5].strip())

    genes_kmers_dict = defaultdict(list)
    genes_pvalues_dict = defaultdict(float)
    for i in range(len(descriptions)):
        match = re.search(r'Name=([\w.-]+)', descriptions[i])
        if match:
            genes_kmers_dict[match.group(1)].append(kmers[i])
            genes_pvalues_dict[match.group(1)] += float(pvalues[kmers2.index(kmers[i])])

    genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))

def summarise2():

    data_frame = pd.read_csv("alignments.csv", sep=",")
    col8 = data_frame.iloc[:, 8]
    col1 = data_frame.iloc[:, 1]
    col2 = data_frame.iloc[:, 2]

    unique_kmers = set(data_frame.iloc[:,1])

    genes = []
    for i in range(len(col8)):
        genes.append(re.findall(r"'Name'\s*:\s*'([^']+)'", col8[i]))

    genes = [item for sublist in genes for item in sublist if item != []]
    unique_genes = list(set(genes))
    

    kmers_coeffs = defaultdict(list)
    total_significance = sum(set(col2))
    print(total_significance)
    for i in range(len(col2)):
        kmers_coeffs[col1[i]] = float(col2[i])


    genes_kmers_dict = defaultdict(list)
    for i in range(len(col8)):
        genes = re.findall(r"'Name'\s*:\s*'([^']+)'", col8[i])
        for gene in genes:
                genes_kmers_dict[gene].append(col1[i])

    genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))

    genes_Ukmers_dict = defaultdict(list)
    genes_significance = defaultdict(list)
    for gene in genes_kmers_dict:
        genes_Ukmers_dict[gene] = len(set(genes_kmers_dict[gene]))
        for kmer in set(genes_kmers_dict[gene]):
            genes_significance[gene] += kmers_coeffs[kmer]

    with open("gene_summary.csv", "w") as csv_summary:
        csv_summary.write(F"Gene({len(unique_genes)}),#Unique k-mers({len(unique_kmers)}), #Total matching alignments\n")
        for key in genes_kmers_dict:
            csv_summary.write(f"{key}, {genes_Ukmers_dict[key]}, {len(genes_kmers_dict[key])}\n")

    return 0

def readPvalue(chi_results):
    kmers = []
    pvalues = []

    with open(chi_results, "r") as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                kmer = parts[0]
                pvalue = parts[2]
                kmers.append(kmer)
                pvalues.append(pvalue)

    return kmers, pvalues

def make_manhattan():
    pass

def main():

    summarise()



if __name__ == "__main__":
    main()



# import matplotlib.pyplot as plt
# plt.imshow(sequence_matrix, cmap='viridis')
# plt.colorbar()
# plt.show()cd 