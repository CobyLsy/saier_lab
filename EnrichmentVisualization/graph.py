import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import cm
import seaborn as sns
import textwrap
import pandas as pd

def fraction_to_float(fraction_str):
    numerator, denominator = fraction_str.split('/')
    return float(numerator) / float(denominator)

df = pd.read_csv('results.tsv', delimiter='\t')

df['ratio_in_study'] = df['ratio_in_study'].apply(fraction_to_float)
df['ratio_in_pop'] = df['ratio_in_pop'].apply(fraction_to_float)

#Calculate enrichment fold by dividing ratio in study by ratio in population
df['fold'] = df['ratio_in_study'] / df['ratio_in_pop']

#plot the first 30 for demonstration
df = df[0:30]

fig, ax = plt.subplots(figsize = (0.5,2.75))

cmap = mpl.cm.bwr_r
norm = mpl.colors.Normalize(vmin = df.p_fdr_bh.min(), vmax = df.p_fdr_bh.max())

mapper = cm.ScalarMappable(norm = norm, cmap = cm.bwr_r)

cbl = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, orientation = 'vertical')

plt.figure(figsize = (2,4))

ax = sns.barplot(data = df, x = 'fold', y = 'name', palette = mapper.to_rgba(df.p_fdr_bh.values))
ax.set_yticklabels([textwrap.fill(e,70) for e in df['name']])

plt.show()
