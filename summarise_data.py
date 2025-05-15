import pandas as pd
from collections import Counter, defaultdict
import re

import matplotlib.pyplot as plt
import numpy as np



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

def summarise(chi_results):

    kmers2, pvalues = readPvalue(chi_results)

    descriptions = []
    locations = []
    kmers = []

    with open("genes.csv", "r") as genes_file:
        lines = genes_file.readlines()

        for line in lines:
            line = line.strip().split(",")
            kmers.append(line[0].strip())
            locations.append(int(line[3].strip()))
            descriptions.append(line[5].strip())
            #genomes.append(line[4].strip())

    genes_kmers_dict = defaultdict(list)
    genes_pvalues_dict = defaultdict(float)
    genes_locations = defaultdict(list)
    genes_genomes_dict = defaultdict(list)
    for i in range(len(descriptions)):
        match = re.search(r'Name=([\w.-]+)', descriptions[i])
        if match:
            genes_kmers_dict[match.group(1)].append(kmers[i])
            genes_pvalues_dict[match.group(1)] += float(pvalues[kmers2.index(kmers[i])])
            genes_locations[match.group(1)].append(str(locations[i]))
            #genes_genomes_dict[match.group(1)].append(str(genomes[i]))

    # genes_Ugenomes_dict = defaultdict(list)
    # for gene, genomes in genes_genomes_dict.items():
    #     res = 0
    #     sus = 0
    #     total = len(genomes)
    #     for genome in genomes:
    #         if genome[0] == "1":
    #             res += 1
    #         else:
    #             sus += 1
    #     genes_Ugenomes_dict[gene].append([total,res,sus])
    


    genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))

    with open("gene_summary_noalign.csv", "w") as summary:
        summary.write(f"Genes {len(genes_kmers_dict.keys())}, sum of p-values, #Unique k-mers, Unique k-mers + locations + p-values\n")
        for gene in genes_kmers_dict:
            kmers = genes_kmers_dict[gene]
            summary.write(f"{gene}, {genes_pvalues_dict[gene]},{len(genes_kmers_dict[gene])}, ({' '.join(genes_kmers_dict[gene])})\n")
            summary.write(f"-----------,-------------,--------------,({' '.join(genes_locations[gene])})\n")
            summary.write(f"-----------,-------------,--------------,({' '.join([pvalues[kmers2.index(kmer)] for kmer in kmers if kmer in kmers2])})\n")

def summarise2(chi_results):

    data_frame = pd.read_csv("alignments.csv", sep=",")
    col8 = data_frame.iloc[:, 8]
    col1 = data_frame.iloc[:, 1]
    col3 = data_frame.iloc[:,3]
    col5 = data_frame.iloc[:, 5]
    #col2 = data_frame.iloc[:, 2]

    unique_kmers = set(data_frame.iloc[:,1])

    genes = []
    for i in range(len(col8)):
        genes.append(re.findall(r"'Name'\s*:\s*'([^']+)'", col8[i]))

    genes = [item for sublist in genes for item in sublist if item != []]
    unique_genes = list(set(genes))
    chi_kmers, chi_pvalues = readPvalue(chi_results)


    # kmers_coeffs = defaultdict(list)
    # total_significance = sum(set(col2))
    # print(total_significance)
    # for i in range(len(col2)):
    #     kmers_coeffs[col1[i]] = float(col2[i])


    genes_kmers_dict = defaultdict(list)
    kmers_locations_dict = defaultdict(list)
    genes_genomes_dict = defaultdict(list)
    for i in range(len(col8)):
        genes = re.findall(r"'Name'\s*:\s*'([^']+)'", col8[i])
        for gene in genes:
                genes_kmers_dict[gene].append(col1[i])
                kmers_locations_dict[col1[i]].append(str(col5[i]))
                genes_genomes_dict[gene].append(col3[i])

    genes_locations_dict = defaultdict(list)

    for gene, kmers in genes_kmers_dict.items():
        kmers = set(kmers)
        genes_kmers_dict[gene] = kmers
        genomes = set(genes_genomes_dict[gene])
        genes_genomes_dict[gene] = genomes
        for kmer in kmers:
            genes_locations_dict[gene].append(kmers_locations_dict[kmer][0] if len(set(kmers_locations_dict[kmer])) < 2 else kmers_locations_dict[kmer][0] + "*")

    genes_Ugenomes_dict = defaultdict(list)
    for gene, genomes in genes_genomes_dict.items():
        res = 0
        sus = 0
        total = len(genomes)
        for genome in genomes:
            if genome[0] == "1":
                res += 1
            else:
                sus += 1
        genes_Ugenomes_dict[gene].append([total,res,sus])

    genes_Ukmers_dict = defaultdict(list)
    for gene in genes_kmers_dict:
        genes_Ukmers_dict[gene] = len(genes_kmers_dict[gene])

    with open("gene_summary_aligned.csv", "w") as csv_summary:
        csv_summary.write(F"Gene({len(unique_genes)}), Min of k-mer p-values, #Unique k-mers({len(unique_kmers)}), #Unique genomes (R/S) , Unique k-mers + locations + pvalues, Unique genomes\n")
        
        genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))

        for key in genes_kmers_dict:
            kmers = genes_kmers_dict[key]
            sum_of_pvalues = min(
                float(chi_pvalues[chi_kmers.index(kmer)]) for kmer in kmers if kmer in chi_kmers
            )
            csv_summary.write(f"{key}, {sum_of_pvalues}, {genes_Ukmers_dict[key]}, {genes_Ugenomes_dict[key][0][0]}({genes_Ugenomes_dict[key][0][1]}/{genes_Ugenomes_dict[key][0][2]}),({' '.join(genes_kmers_dict[key])}), {' '.join(genes_genomes_dict[key])}\n")
            csv_summary.write(f"----------,---------------,-------------,-------------,({' '.join(genes_locations_dict[key])}), ----------------------------------------------------------------------------------------------------\n")
            csv_summary.write(f"----------,---------------,-------------,-------------,({' '.join([chi_pvalues[chi_kmers.index(kmer)] for kmer in kmers if kmer in chi_kmers])})\n")
    
    with open("positions_pvalues.csv", "w") as manhattan_data:
        manhattan_data.write("kmer,position,pval\n")
        for i in range(1, len(chi_kmers)):
            if list(set(kmers_locations_dict[chi_kmers[i]])) == []:
                manhattan_data.write(f"{chi_kmers[i]},-1,{chi_pvalues[i]}\n")
            else:
                manhattan_data.write(f"{chi_kmers[i]},{list(set(kmers_locations_dict[chi_kmers[i]]))[0]},{chi_pvalues[i]}\n")
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

def make_manhattan(pos_pval):
    # Step 1: Load the data
    df = pd.read_csv(pos_pval)
    # Step 2: Filter invalid positions
    df = df[df["position"] != -1]

    # Step 3: Compute -log10(pval)
    df["neg_log10_p"] = -np.log10(df["pval"])

    # Step 4: Sort by position
    df = df.sort_values("position")

    # Step 5: Plot
    plt.figure(figsize=(10, 5))
    plt.scatter(df["position"], df["neg_log10_p"], c="blue", s=10)
    plt.title("Manhattan Plot")
    plt.xlabel("Genomic Position")
    plt.ylabel("-log10(p-value)")
    plt.grid(True)
    plt.xlim(0, 3149213)
    plt.savefig("manhattan_plot.png")
    plt.show()

def main():

    summarise("chi2_results_Isoniazid_top1000.tsv")
    summarise2("chi2_results_Isoniazid_top1000.tsv")

    make_manhattan("positions_pvalues.csv")

if __name__ == "__main__":
    main()



# import matplotlib.pyplot as plt
# plt.imshow(sequence_matrix, cmap='viridis')
# plt.colorbar()
# plt.show()cd 