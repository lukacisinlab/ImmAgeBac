import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd


NUMBER_OF_PATHWAYS = 7

gsea_results = pd.read_csv("gsea_results_pearson_gbp_v_0.2.csv")
gsea_results = gsea_results.sort_values('NES')
print(gsea_results.head(7)[['NES', 'FDR q-val']])
print(gsea_results.tail(7)[['NES', 'FDR q-val']])


plt.figure(figsize=(14, 10)) 
subset = pd.concat([gsea_results.head(NUMBER_OF_PATHWAYS), gsea_results.tail(NUMBER_OF_PATHWAYS)]) 
print(subset[['NES', 'FDR q-val']])
ax = barplot = sns.barplot(x='NES', y='Term', data = subset, color="black")

for i, (value, fdr) in enumerate(zip(subset['NES'], subset['FDR q-val'])):
    ax.text(
        x=value + 0.1 if value > 0 else value - 0.1,
        y=i,
        s=f"FDR={fdr:.5f}",
        va='center',
        ha='left' if value > 0 else 'right',
        color='black',
        fontsize=8,
    )

plt.xlabel('Normalized Enrichment Score (NES)')
plt.ylabel('Pathway')
plt.title('GSEA Results for Pearson Correlation with GO Biological Process (GBP)', fontsize=16)

sns.despine(left=True, bottom=True)  
plt.grid(axis='x', linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('pearson_gbp_bar_chart.pdf', format='pdf')
