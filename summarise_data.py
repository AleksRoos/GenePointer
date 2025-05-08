import pandas as pd

def main():


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
        




if __name__ == "__main__":
    main()



# import matplotlib.pyplot as plt
# plt.imshow(sequence_matrix, cmap='viridis')
# plt.colorbar()
# plt.show()cd 