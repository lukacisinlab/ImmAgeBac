#%%
import matplotlib.pyplot as plt
import pandas as pd
import scipy.stats
import seaborn as sns
import gseapy as gp

data = pd.read_excel("graph_excel.xlsx")
data.sort_values(by="Real val", inplace=True)

left_values = [min(1, x) for x in data["Real val"]]

width_values = abs(data["Real val"] - 1)

fig, ax = plt.subplots(figsize=(5, 5))
sns.despine(fig=fig)
ax.barh(data["Name"], width_values, color='black', left=left_values, log=True)
plt.savefig("./imm_age_logscale.svg")
# %%
expression_profiles = pd.read_csv("expresion_profiles.csv")

interesting_genes = ["HMGB1", "HMGB2", "TNIP2", "MYD88", "CASP1"] # favorite HMGB2 TNIP2 MYD88


all_genes = list(expression_profiles.columns)
all_genes.remove("Unnamed: 0")
all_genes = set(all_genes)

genes_lpps = set(gp.read_gmt("/home/hruzko/Desktop/GitHub/MAD-BAC/prerank_output_correlation_df_spearman_sorted_GO_Biological_Process_2023/gene_sets.gmt")["Cellular Response To Lipopolysaccharide (GO:0071222)"])
all_genes = all_genes.difference(genes_lpps)
print(all_genes.intersection(genes_lpps))

for gene in interesting_genes: #all_genes
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.despine(fig=fig)
    plt.title(f"Mean gene expression for gene: {gene}")
    
    ax.barh(
        expression_profiles["Unnamed: 0"], 
        expression_profiles[gene], 
        color='black',
        log=True
    )
    
    
    ax.invert_yaxis()
    
    plt.savefig(f"./interesting_genes_good_log/{gene}_barchart.svg") #all_genes_log
    
# %%
expression_profiles = pd.read_csv("expresion_profiles.csv", index_col=0)
expression_profiles["Sum"] = expression_profiles.sum(axis=1, numeric_only=True)
fig, ax = plt.subplots(figsize=(5, 5))
ax.barh(expression_profiles.index, expression_profiles["Sum"])
ax.invert_yaxis()
# %%
import scipy
print(scipy.stats.spearmanr(expression_profiles["Sum"], list(range(1, 21))))

# %%
