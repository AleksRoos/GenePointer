import pandas as pd
from collections import defaultdict
import re
import os

import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq

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

def get_kmer_prevalence(kmers, kmers_genomes_sequences):
    with open(kmers_genomes_sequences, "r") as file:
        lines = file.readlines()
        kmer_prev_line = lines[0].strip().split(",")
        kmer_prevalence_dict = defaultdict(list)
        for kmer in kmers:
            kmer = Seq(kmer)
            reverse_complement = str(kmer.reverse_complement())
            kmer = str(kmer)
            t,r,s = 0,0,0
            for item in kmer_prev_line:
                if kmer in item or reverse_complement in item:
                    match = re.search(r'T:(\d+); R:(\d+)/S:(\d+)', item)
                    if match:
                        t, r, s = match.groups()
            kmer_prevalence_dict[kmer] = [str(t),str(r),str(s)]
    return kmer_prevalence_dict

def summarise_intergenic_of_mapped_kmers(chi_kmers, chi_pvals, kmer_prevalence_dict):
    kmers2, pvalues = chi_kmers, chi_pvals

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

            #genomes.append(line[4].strip())
    with open("intergenic_summary.csv", "w") as summary:
        summary.write(f"inter_Start...End, Signif. k-mer p-value, Signif. k-mer prevalence, K-mer, Location\n")
        for key,value in intergenicID_dict.items():
            summary.write(f"{key},{value[0][1]},{'/'.join(kmer_prevalence_dict[value[0][0]])},{value[0][1]},{value[0][2]}\n")

def summarise_genes_of_mapped_kmers(chi_kmers, chi_pvals, kmer_prevalence_dict):

    kmers2, pvalues = chi_kmers, chi_pvals

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
            summary.write(f"{gene},{min_pval},{'/'.join(kmer_prevalence_dict[min_pval_kmer])},{len(genes_kmers_dict[gene])},{' '.join(genes_kmers_dict[gene])},{' '.join(genes_locations[gene])}\n")

def summarise_aligned_kmers(chi_kmers, chi_pvals, kmer_prevalence_dict):

    data_frame = pd.read_csv("alignments.csv", sep=",")
    descriptions = data_frame.iloc[:, 8]
    kmer_names = data_frame.iloc[:, 1]
    genomes = data_frame.iloc[:,3]
    locations = data_frame.iloc[:, 5]
    chi_kmers, chi_pvalues = chi_kmers, chi_pvals

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

def summarise_unaligned_kmers(chi_kmers, chi_pvals, kmer_prevalence_dict):
    data_frame = pd.read_csv("no_alignments.csv", sep=",")
    kmers = data_frame.iloc[:, 1]
    sequences = data_frame.iloc[:, 6]	
    chi_kmers, chi_pvalues = chi_kmers, chi_pvals

    kmer_pvals_dict = defaultdict(float)
    for i in range(len(kmers)):
        kmer_pvals_dict[kmers[i]] = chi_pvalues[chi_kmers.index(kmers[i])]


    with open("unaligned_summary.csv", "w") as summary:
        summary.write(f"K-mer, P-value, Prevalence, Extended Sequence\n")
        lines = []
        for i in range(len(kmers)):
            kmer = kmers[i]
            pval = kmer_pvals_dict[kmer]
            prevalence = kmer_prevalence_dict[kmer]
            sequence = sequences[i]
            lines.append(f"{kmer}, {pval}, {'/'.join(prevalence)}, {sequence}\n")
        lines = set(lines)
        for line in lines:
            summary.write(line)

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

def extract_antibiotics_from_folder(folder='.'):
    pattern = re.compile(r'k-mers_and_coefficients_in_(.+?)_model_(.+?)\.txt')
    antibiotics = []

    for fname in os.listdir(folder):
        match = pattern.match(fname)
        if match:
            classifier, antibiotic = match.groups()
            antibiotics.append(antibiotic)
            print(f"{fname} → Antibiotic: {antibiotic}")
    
    return antibiotics[0]

def summarise_all(antibiotic = None):
    print("---------- Summarising results ----------")

    try:
        antibiotic = extract_antibiotics_from_folder(".").lower().capitalize()
        #print(f"Extracted antibiotic: {antibiotic}")
    except Exception as e:
        # If extraction fails, default to a specific antibiotic
        antibiotic = input("Enter the antibiotic name (e.g., ethambutol): ").strip().lower().capitalize()
    print(f"    Using antibiotic: {antibiotic}")

    chi_kmers, chi_pvals = readPvalue(f"chi2_results_{antibiotic}_top1000.tsv")
    kmer_prevalence_dict = get_kmer_prevalence(chi_kmers, "kmers_genomes_sequences_table.csv")

    #summarise_genes_of_mapped_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)
    #summarise_intergenic_of_mapped_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)
    
    summarise_aligned_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)
    summarise_unaligned_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)

    try:
        make_manhattan("positions_pvalues.csv")
    except:
        print("    Skipping Manhattan plot generation.")