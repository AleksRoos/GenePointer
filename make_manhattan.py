import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Step 1: Load the data
df = pd.read_csv("positions_pvalues.csv")
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
plt.show()
