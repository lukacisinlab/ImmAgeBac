#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns



df_res_s_GBP = pd.read_csv("/home/hruzko/Desktop/GitHub/MAD-BAC/figures/fig_one/gseapy_results_viz/reactome_gsea.csv")
gsea_results = df_res_s_GBP[['Term', 'NES', 'NOM p-val', 'FDR q-val', 'Lead_genes', 'Tag %', 'Gene %']].sort_values('NES', ascending=False)

plt.figure(figsize=(5, 5))  
barplot = sns.barplot(x='NES', y='Term', data=pd.concat([gsea_results.head(5), gsea_results.tail(5)]), color="black")



plt.xlabel('Normalized Enrichment Score (NES)')
plt.ylabel('Pathway')
plt.title('Top Enriched Pathways')


sns.despine()
plt.grid(axis='x', linestyle='--', alpha=0.7)
plt.tight_layout()

plt.savefig("./reactome_enrichement.svg")
plt.show()
# %%
