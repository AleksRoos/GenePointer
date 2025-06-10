import pandas as pd
from collections import defaultdict
import re
import os
from gene_pointer import analysis
import pandas as pd
from glob import glob

import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
import ast


#USED MINREQ
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
#USED MINREQ
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



#USED MINREQ
def summarise_aligned_kmers(alignments_file, chi_kmers, chi_pvals, kmer_prevalence_dict):

    data_frame = pd.read_csv(alignments_file, sep=",")
    descriptions = data_frame.iloc[:, 8]
    kmer_names = data_frame.iloc[:, 1]
    genomes = data_frame.iloc[:,3]
    locations = data_frame.iloc[:, 5]
    chi_kmers, chi_pvalues = chi_kmers, chi_pvals
    seqID = alignments_file.split("_")[-2]
    print(seqID)

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

    with open(f"Summary_aligned_kmers{seqID}.csv", "w") as csv_summary:
        csv_summary.write(F"Genes (Total: {len(genes_kmers_dict.keys())}), Minimum k-mer p-value, Minimum p-value kmer prevalence (Tot/Res/Sus), #Unique k-mers({len(kmer_prevalence_dict.keys())}), Gene prevalence in genomes (T/R/S) , Unique k-mers, Minimum p-value k-mer location, Unique genomes\n")
        
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
#USED MINREQ
def summarise_unaligned_kmers(unaligned_file, chi_kmers, chi_pvals, kmer_prevalence_dict):
    data_frame = pd.read_csv(unaligned_file, sep=",")
    kmers = data_frame.iloc[:, 1]
    sequences = data_frame.iloc[:, 6]	
    chi_kmers, chi_pvalues = chi_kmers, chi_pvals

    kmer_pvals_dict = defaultdict(float)
    for i in range(len(kmers)):
        kmer_pvals_dict[kmers[i]] = chi_pvalues[chi_kmers.index(kmers[i])]


    with open("Summary_unaligned_kmers.csv", "w") as summary:
        summary.write(f"K-mer (Total: {len(kmers)}), P-value, Prevalence, Extended Sequence\n")
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


def extract_info(entry):
    # 1. FeatureType: first element before the first comma
    feature_type = entry.split(",")[0].strip()

    # 2. All 'Name' values using regex
    names = re.findall(r"'Name':\s*'([^']+)'", entry)
    name_joined = "-".join(names) if names else None

    # 3. First ID inside 'Note'
    id_match = re.search(r"ID:([^%~;\s]+)", entry)
    gene_id = id_match.group(1) if id_match else None

    return pd.Series([feature_type, name_joined, gene_id])

def summarise_multi_refseq_alignments(alignments_folder, chi_kmers, chi_pvals, kmer_prevalence_dict):
    
    
    
    # Load all matching files for aligned sequences
    files = glob(f"{alignments_folder}/Aligned*.csv")
    df_list = [pd.read_csv(file) for file in files]
    alignments_data = pd.concat(df_list, ignore_index=True)
    

    # Load all matching files for unaligned sequences
    files = glob(f"{alignments_folder}/Un_aligned*.csv")
    df_list = [pd.read_csv(file) for file in files]
    unaligned_data = pd.concat(df_list, ignore_index=True)

    # Make k-mer p-values dict from chi_kmers and chi_pvals
    kmer_pval_dict = dict(zip(chi_kmers, chi_pvals))

    #Add k-mer p-values to dataframes
    alignments_data["K-mer p-value"] = alignments_data["K-mer"].map(kmer_pval_dict)
    unaligned_data["K-mer p-value"] = unaligned_data["K-mer"].map(kmer_pval_dict)

    #Add k-mer prevalence to dataframes
    alignments_data["K-mer prevalence"] = alignments_data["K-mer"].map(kmer_prevalence_dict)
    unaligned_data["K-mer prevalence"] = unaligned_data["K-mer"].map(kmer_prevalence_dict)

    # Make values numeric 
    alignments_data["K-mer p-value"] = pd.to_numeric(alignments_data["K-mer p-value"], errors="coerce")
    unaligned_data["K-mer p-value"] = pd.to_numeric(unaligned_data["K-mer p-value"], errors="coerce")
    alignments_data["AlignmentLeftSide_Pos"] = pd.to_numeric(alignments_data["AlignmentLeftSide_Pos"], errors="coerce")
    
    # Move useful annotation values from "Genes from base sequence" into separate columns
    alignments_data[['FeatureType', 'Name', 'ID']] = alignments_data['Genes from base sequence'].apply(extract_info)

    # Move p-value column next to k-mer column
    cols = alignments_data.columns.tolist()
    kmer_index = cols.index("K-mer")
    cols.insert(kmer_index + 1, cols.pop(cols.index("K-mer p-value")))
    alignments_data = alignments_data[cols]

    # Move p-value column next to k-mer column in unaligned data
    cols = unaligned_data.columns.tolist()
    kmer_index = cols.index("K-mer")
    cols.insert(kmer_index + 1, cols.pop(cols.index("K-mer p-value")))
    unaligned_data = unaligned_data[cols]




    min_pos = alignments_data["AlignmentLeftSide_Pos"].min()
    max_pos = alignments_data["AlignmentLeftSide_Pos"].max()
    bin_width = 4000
    # Create bins based on the min and max positions
    bin_edges = list(range(int(min_pos) // bin_width * bin_width,
                        int(max_pos) // bin_width * bin_width + bin_width,
                        bin_width))

    # Use pd.cut with custom integer bins
    alignments_data["PositionBin"] = pd.cut(
        alignments_data["AlignmentLeftSide_Pos"],
        bins=bin_edges,
        right=False  # Optional: choose whether the right edge is inclusive
    )


    #alignments_data["PositionBin"] = pd.cut(alignments_data["AlignmentLeftSide_Pos"], bins=1000)
    position_summary = alignments_data.groupby("PositionBin")["K-mer"].count()
    position_summary.to_csv("position_kmer_count_summary.csv", header=True)
    
    # Step 3: Group by bin and compute average p-value
    epsilon = 1
    alignments_data["K-mer log p-value"] = -np.log10(alignments_data["K-mer p-value"].fillna(epsilon))
    position_pval_avg = alignments_data.groupby("PositionBin")["K-mer log p-value"].sum()
    # Step 4: Export to CSV
    position_pval_avg.to_csv("position_kmer_pval_summary.csv", header=True)


    def get_top_genes(series):
        num_genes = 5
        counts = series.dropna().value_counts()
        if counts.empty:
            return "Unknown"
        return ", ".join(counts.head(num_genes).index)

    # Get 3 most common genes per position bin
    top_genes = (
        alignments_data
        .groupby("PositionBin")["Name"]
        .agg(get_top_genes)
        .reset_index()
        .rename(columns={"Name": "TopGenes"})
    )






    binned_df = pd.read_csv("position_kmer_count_summary.csv")
    binned_df["PositionBin"] = binned_df["PositionBin"].astype(str)
    
    binned_df = position_pval_avg.reset_index()
    binned_df.columns = ["PositionBin", "K-mer"]
    binned_df = binned_df.merge(top_genes, on="PositionBin", how="left")
    binned_df.to_csv("position_kmer_count_summary.csv", index=False)

    plt.figure(figsize=(20, 5))
    plt.bar(range(len(binned_df)), binned_df["K-mer"], color="skyblue")
    
    N = 1
    for i, (gene, value) in enumerate(zip(binned_df["TopGenes"], binned_df["K-mer"])):
        if i % N == 0:  # show label every N bars
            plt.text(i, value + 0.1, gene, ha='center', va='bottom', fontsize=1, rotation=90)

    plt.xticks(
        ticks=range(0, len(binned_df), N),
        labels=binned_df["PositionBin"][::N],
        rotation=90,
        fontsize=1
    )
    plt.xlabel("Position Bins across genome")
    plt.ylabel("K-mer Count")
    plt.title("K-mer Distribution by Binned Positions with Top Genes")
    plt.tight_layout()
    plt.savefig("kmer_count_position_binned_histogram_with_genes.png", dpi=1000)
    plt.close()

    binned_df = pd.read_csv("position_kmer_pval_summary.csv")
    binned_df["PositionBin"] = binned_df["PositionBin"].astype(str)

    # Step 2: Merge with binned log p-value dataframe
    binned_df = position_pval_avg.reset_index()
    binned_df.columns = ["PositionBin", "K-mer log p-value"]
    binned_df = binned_df.merge(top_genes, on="PositionBin", how="left")

    binned_df.to_csv("position_kmer_pval_summary.csv", index=False)

    plt.figure(figsize=(20, 5))
    plt.bar(range(len(binned_df)), binned_df["K-mer log p-value"], color="skyblue")
    
    N = 1
    for i, (gene, value) in enumerate(zip(binned_df["TopGenes"], binned_df["K-mer log p-value"])):
        if i % N == 0:  # show label every N bars
            plt.text(i, value + 0.1, gene, ha='center', va='bottom', fontsize=1, rotation=90)

    
    # Custom ticks
    plt.xticks(
        ticks=range(0, len(binned_df), N),
        labels=binned_df["PositionBin"].astype(str)[::N],
        rotation=90,
        fontsize=1
    )

    plt.xlabel("Position Bins across genome")
    plt.ylabel("K-mer log p-val sum")
    plt.title("K-mer Distribution by Binned Alignment Position with Top Genes")
    plt.tight_layout()
    plt.savefig("kmer_sum_pval_position_binned_histogram_with_genes.png", dpi=1000)
    plt.close()

    alignments_data.to_csv("Aligned_kmer_results.csv", index=False)
    unaligned_data.to_csv("Un_aligned_kmer_results.csv", index=False)


#USED
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
#USED
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
#USED MINREQ
def summarise_all(antibiotic = None):
    print("---------- Summarising results ---------- OUTPUT FILES: Summary_aligned_kmers.csv, Summary_unaligned_kmers.csv")

    
    try:
        antibiotic = extract_antibiotics_from_folder(".").lower().capitalize()
        #print(f"Extracted antibiotic: {antibiotic}")
    except Exception as e:
        # If extraction fails, default to a specific antibiotic
        antibiotic = input("Enter the antibiotic name (e.g., ethambutol): ").strip().lower().capitalize()
    print(f"    Using antibiotic: {antibiotic}")

    chi_kmers, chi_pvals = readPvalue(f"chi2_results_{antibiotic}_top1000.tsv")
    kmer_prevalence_dict = get_kmer_prevalence(chi_kmers, "kmers_genomes_sequences_table.csv")

    summarise_multi_refseq_alignments("Alignments", chi_kmers, chi_pvals, kmer_prevalence_dict)


    # #summarise_genes_of_mapped_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)
    # #summarise_intergenic_of_mapped_kmers(chi_kmers, chi_pvals,kmer_prevalence_dict)
    # alignment_dir = "Alignments"
    # alignments_files_dict = analysis.group_files_by_suffix(alignment_dir)
    # alignment_file = ""
    # unaligned_file = ""
    # for suffix, pair in alignments_files_dict.items():
    #     print(f"\n    Summarising pair for seqID '{suffix}': {pair}")
    #     for filename in pair:
    #         # You can set different separators depending on the file name
    #         if filename[-4:] == ".gff":
    #             GFF = os.path.join(alignment_dir, filename)
    #         elif filename[-4:] == ".fna":
    #             REFSEQ = os.path.join(alignment_dir, filename)

    # summarise_aligned_kmers(alignment_file, chi_kmers, chi_pvals,kmer_prevalence_dict)
    # summarise_unaligned_kmers(unaligned_file, chi_kmers, chi_pvals,kmer_prevalence_dict)

    # try:
    #     make_manhattan("positions_pvalues.csv")
    # except:
    #     print("    Skipping Manhattan plot generation.")