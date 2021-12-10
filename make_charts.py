import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('perf/data.csv')
df['Seconds'] = ((df['Day']-1)*86400+df['Hour']) / 86400

# df.plot(kind='scatter', x='Seconds', y='Queries')

sns.boxplot(x='Day', y='Queries', data=df, meanline=True,
            showcaps=True, showfliers=False, color='blue')
sns.despine(offset=10, trim=True)
plt.show()
