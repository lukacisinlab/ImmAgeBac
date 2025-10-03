#%%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
#%%

df = pd.read_excel("cell_counts.xlsx", index_col="Unnamed: 0")
df_no_bias = df.loc[~df["Donor"].isin(["P1", "P2", "P4"])]
df_no_bias
#%%
SIZE_DEFAULT = 14
SIZE_LARGE = 16
plt.rc("font", size=SIZE_DEFAULT)
plt.rc("axes", titlesize=SIZE_LARGE)
plt.rc("axes", labelsize=SIZE_LARGE)  
plt.rc("xtick", labelsize=SIZE_DEFAULT - 5)  
plt.rc("ytick", labelsize=SIZE_DEFAULT)  

cell_types = df_no_bias['Cell_Type'].unique()
mapping = {cell: i for i, cell in enumerate(cell_types)}
df_no_bias['x_numeric'] = df_no_bias['Cell_Type'].map(mapping)

np.random.seed(42)  
df_no_bias['x_jitter'] = df_no_bias['x_numeric'] + np.random.uniform(-0.15, 0.15, size=len(df_no_bias))

fig, ax = plt.subplots(figsize=(12, 8))
sns.scatterplot(
    data=df_no_bias, 
    x='x_jitter', 
    y='Count', 
    style='Donor', 
    markers=["^", "o", "s", "d", "X"],
    color='black',  
    legend='full'
)

ax.set_xticks(list(mapping.values()))
ax.set_xticklabels(list(mapping.keys()), rotation=45, ha='right')

sns.despine(top=True, right=True, left=True)
plt.title('Swarm Plot of Immune Cell Counts per Donor', fontsize='large')
plt.ylabel('Count frequency')
plt.xlabel('Cell type')
plt.tight_layout()
plt.savefig(fname = 'swarmplot_updated.png',format = 'png', orientation = 'portrait')
plt.show()

# %%
