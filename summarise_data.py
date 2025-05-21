import pandas as pd
from collections import Counter, defaultdict
import re

import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq



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

def summarise4(chi_results):
    kmers2, pvalues = readPvalue(chi_results)

    locations = []
    kmers = []
    starts =[]
    ends = []
    with open("intergenic.csv", "r") as intergenic_file:
        lines = intergenic_file.readlines()
        for line in lines:
            line = line.strip().split(",")
            kmers.append(line[0].strip())
            locations.append(str(line[3].strip()))
            starts.append(str(line[2].strip()))
            ends.append(str(line[4].strip()))
        
    intergenicID_dict = defaultdict(list)
    for i in range(len(locations)):
        intergenic_region_name = "inter_" + starts[i] + "..." + ends[i]
        intergenic_kmer_pvla = pvalues[kmers2.index(kmers[i])]
        intergenic_kmer_location = locations[i]
        intergenic_kmer = kmers[i]
        intergenicID_dict[intergenic_region_name].append([intergenic_kmer, intergenic_kmer_pvla, intergenic_kmer_location])
    
    kmer_prevalence_dict = get_kmer_prevalence(kmers, "kmers_genomes_sequences_table.csv")

            #genomes.append(line[4].strip())
    with open("intergenic_summary.csv", "w") as summary:
        summary.write(f"inter_Start...End, Signif. k-mer p-value, Signif. k-mer prevalence, K-mer, Location\n")
        for key,value in intergenicID_dict.items():
            summary.write(f"{key},{value[0][1]},{"/".join(kmer_prevalence_dict[value[0][0]])},{value[0][1]},{value[0][2]}\n")

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
        if descriptions[i].strip().split(",")[0] == "intergenic":
            match = str("intergenic_" + descriptions[i].strip().split(",")[2].strip())
        else:
            match = re.search(r'Name=([\w.-]+)', descriptions[i])
        if match:
            genes_kmers_dict[match.group(1)].append(kmers[i])
            genes_pvalues_dict[match.group(1)] += float(pvalues[kmers2.index(kmers[i])])
            genes_locations[match.group(1)].append(str(locations[i]))
    
    kmer_prevalence_dict = get_kmer_prevalence(kmers, "kmers_genomes_sequences_table.csv")


    genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))

    with open("gene_summary_noalign.csv", "w") as summary:
        summary.write(f"Genes {len(genes_kmers_dict.keys())}, Signif. k-mer p-value, Signif. k-mer prevalence, #Unique k-mers, Unique k-mers, Locations\n")
        for gene in genes_kmers_dict:
            kmers = genes_kmers_dict[gene]
            p_values = []
            for kmer in kmers:
                p_values.append(pvalues[kmers2.index(kmer)])
            min_pval_kmer = kmers[p_values.index(min(p_values))]
            min_pval = min(p_values)
            summary.write(f"{gene},{min_pval},{"/".join(kmer_prevalence_dict[min_pval_kmer])},{len(genes_kmers_dict[gene])},{' '.join(genes_kmers_dict[gene])},{' '.join(genes_locations[gene])}\n")

def summarise2(chi_results):

    data_frame = pd.read_csv("alignments.csv", sep=",")
    descriptions = data_frame.iloc[:, 8]
    kmer_names = data_frame.iloc[:, 1]
    genomes = data_frame.iloc[:,3]
    locations = data_frame.iloc[:, 5]
    chi_kmers, chi_pvalues = readPvalue(chi_results)
    kmer_prevalence_dict = get_kmer_prevalence(kmer_names, "kmers_genomes_sequences_table.csv")

    genes_kmers_dict = defaultdict(list)
    kmers_locations_dict = defaultdict(dict)
    genes_genomes_dict = defaultdict(list)
    for i in range(len(descriptions)):
        if descriptions[i].strip().split(",")[0] == "intergenic":
            genes = [str("inter_" + descriptions[i].strip().split(",")[2].strip() + "..." + descriptions[i].strip().split(",")[4].strip()[:-1])]
        else:
            genes = re.findall(r"'Name'\s*:\s*'([^']+)'", descriptions[i])
        for gene in genes:
            genes_kmers_dict[gene].append(kmer_names[i])
            genes_genomes_dict[gene].append(genomes[i])
            kmers_locations_dict[gene][kmer_names[i]] = locations[i]

    for gene, kmers in genes_kmers_dict.items():
        kmers = set(kmers)
        genes_kmers_dict[gene] = kmers
        genomes = set(genes_genomes_dict[gene])
        genes_genomes_dict[gene] = genomes

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
        genes_Ugenomes_dict[gene].append([str(total),str(res),str(sus)])

    genes_Ukmers_dict = defaultdict(list)
    for gene in genes_kmers_dict:
        genes_Ukmers_dict[gene] = len(genes_kmers_dict[gene])

    with open("gene_summary_aligned.csv", "w") as csv_summary:
        csv_summary.write(F"Gene({len(genes_kmers_dict.keys())}), Min of k-mer p-values, Signif. kmer prevalence (T/R/S), #Unique k-mers({len(kmer_prevalence_dict.keys())}), Gene prevalence in genomes (T/R/S) , Unique k-mers, Signif. k-mer location, Unique genomes\n")
        
        genes_kmers_dict = dict(sorted(genes_kmers_dict.items(), reverse=True))
        for key in genes_kmers_dict:
            kmers = list(genes_kmers_dict[key])
            p_values = []
            for kmer in kmers:
                p_values.append(chi_pvalues[chi_kmers.index(kmer)])
            min_of_pvalues = min(p_values)
            min_pval_kmer = kmers[p_values.index(min_of_pvalues)]
            csv_summary.write(f"{key}, {min_of_pvalues}, {'/'.join(kmer_prevalence_dict[min_pval_kmer])}, {genes_Ukmers_dict[key]}, {'/'.join(genes_Ugenomes_dict[key][0])},{' '.join(genes_kmers_dict[key])},{kmers_locations_dict[key][min_pval_kmer]},{' '.join(genes_genomes_dict[key])}\n")     
    
    # with open("positions_pvalues.csv", "w") as manhattan_data:
    #     manhattan_data.write("kmer,position,pval\n")
    #     for i in range(1, len(chi_kmers)):
    #         if list(set(kmers_locations_dict[chi_kmers[i]])) == []:
    #             manhattan_data.write(f"{chi_kmers[i]},-1,{chi_pvalues[i]}\n")
    #         else:
    #             manhattan_data.write(f"{chi_kmers[i]},{list(set(kmers_locations_dict[chi_kmers[i]]))[0]},{chi_pvalues[i]}\n")
    return 0

def get_kmer_prevalence(kmers, kmers_genomes_sequences):
    with open(kmers_genomes_sequences, "r") as file:
        lines = file.readlines()
        kmer_prev_line = lines[0].strip().split(",")
        kmer_prevalence_dict = defaultdict(list)
        for kmer in kmers:
            kmer = Seq(kmer)
            reverse_complement = str(kmer.reverse_complement())
            kmer = str(kmer)
            for item in kmer_prev_line:
                if kmer in item or reverse_complement in item:
                    match = re.search(r'T:(\d+); R:(\d+)/S:(\d+)', item)
                    if match:
                        t, r, s = match.groups()
            kmer_prevalence_dict[kmer] = [str(t),str(r),str(s)]
    return kmer_prevalence_dict
 
def summarise3(chi_results):
    data_frame = pd.read_csv("alignments.csv", sep=",")
    kmers = data_frame.iloc[:, 1]
    sequences = data_frame.iloc[:, 6]	
    chi_kmers, chi_pvalues = readPvalue(chi_results)
    kmers_prevalence = get_kmer_prevalence(kmers, "kmers_genomes_sequences_table.csv")

    kmer_pvals_dict = defaultdict(float)
    for i in range(len(kmers)):
        kmer_pvals_dict[kmers[i]] = chi_pvalues[chi_kmers.index(kmers[i])]


    with open("unaligned_summary.csv", "w") as summary:
        summary.write(f"K-mer, P-value, Prevalence, Extended Sequence\n")
        lines = []
        for i in range(len(kmers)):
            kmer = kmers[i]
            pval = kmer_pvals_dict[kmer]
            prevalence = kmers_prevalence[kmer]
            sequence = sequences[i]
            lines.append(f"{kmer}, {pval}, {"/".join(prevalence)}, {sequence}\n")
        lines = set(lines)
        for line in lines:
            summary.write(line)


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

    antibiotic = "Ciprofloxacin"

    summarise(f"chi2_results_{antibiotic}_top1000.tsv")
    summarise2(f"chi2_results_{antibiotic}_top1000.tsv")
    summarise3(f"chi2_results_{antibiotic}_top1000.tsv")
    summarise4(f"chi2_results_{antibiotic}_top1000.tsv")

    try:
        make_manhattan("positions_pvalues.csv")
    except:
        print("skipping Manhattan plot generation.")
if __name__ == "__main__":
    main()



# import matplotlib.pyplot as plt
# plt.imshow(sequence_matrix, cmap='viridis')
# plt.colorbar()
# plt.show()cd 