import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import r2_score


'''
#ETE3 Vs GCf
x=[226,437,1952,1952,437,637,1358,554,554,1358,894,1527,625,1527,585,585,637,708,120,715,127]
y= [96.4,96.6,26.9,51.3,90.5,92.9,76.9,57.3,90.3,71.7,58.9,68.7,40.2,35.4,34.7,52,84.2,56.8,85.5,54.9,74.5]
'''


'''
#reconcILS NNI all paralogy vs NNI SCO

x=[0,2509,1121,5955,9758,4018,3442,6307,2181,6295,2214,3135,3971,1011,4591,7389,9768,5881,6906,13187,4965,7071]
y=[0,210,77,859,1368,474,440,885,301,872,301,440,615,123,654,1121,1376,803,899,1922,699,965]

 9: (3,3)
 8: (3,2)
 7 :(default)
 6: (2,2)
 5:(5,5)

'''

# reconciLS NNI vs gcf smith et al .
#x=[210,77,1368,859,474,440,885,301,872,301,440,615,123,654,1121,1376,803,899,1922,699,965]
#y= [96.4,96.6,39.1,39.1,93.5,92.9,76.9,57.3,90.3,71.7,58.9,68.7,40.2,35.4,34.7,52,84.2,56.8,85.5,54.9,74.5]
y= [96.4,96.6,26.9,51.3,90.5,84.2,56.8,85.5,54.9,74.5,92.9,76.9,90.3,57.3,71.7,58.9,68.7,40.2,35.4,52,34.7]

#9
x=[219,79,1370,860,485,457,891,307,880,303,466,639,130,662,1153,1393,818,907,1945,704,976]
#8
#x=[223,85,1371,863,500,462,897,312,881,304,488,670,137,673,1184,1409,828,910,1957,705,982]

#old param
x=[223,86,1374,865,499,465,896,312,881,304,489,676,144,693,1178,1408,855,921,1941,746,1040]

#5
x=[223,85,1731,863,500,462,897,312,881,304,488,670,134,673,1184,1409,828,910,1957,705,982]

'''
#Concordance correlation between dupcoal and simphy
x= [99.99,86.3,69.8,49,50.1,26.1,29.2,43.1,86.8,81.1,40,73.2,38.6,55.6,54.7,98.4,99.6,49.3,66.6,97.3,60.8,56.1,56.5,64,54.9,54.4,88.6,82,72.8,39.1,32.5,45.7,76.4,66.5,49.3,54,50.6,59.7]
y=[99.6,86.6,72,49.8,48.3,29.1,29.5,45.2,85.3,82.9,41.2,73.4,37.8,56.1,54,97.5,99.4,48.9,65.3,96.9,59.2,58.3,59.6,61.4,55.3,53.5,89.4,80.7,72.1,39.6,31.8,44.8,75.8,64.2,48.4,49.5,56.7,54.5]
'''

print(len(x))
print(len(y))

for i in range(len(y)):
    print(i, y[i],x[i])

x =[abs(i) for i in x]
y =[abs(i) for i in y]

print(len(x)==len(y))
#x = [(1-(i / 1810))*100 for i in x]

df = pd.DataFrame({'X': x, 'Y': y})

sns.set(style="whitegrid")
plt.figure(figsize=(10, 6))
sns.regplot(x='X', y='Y', data=df, scatter_kws={'s': 50}, line_kws={'color': 'red'})
#plt.title('Correlation between X and Y')
#plt.xlabel('X')
#plt.ylabel('Y')

plt.title('Correlation between NNI in reconcILS and gCF in IQ-Tree Dup_cost=5 , Loss_cost=5')
plt.xlabel('NNI in reconcILS')
plt.ylabel('gCF in IQ-Tree')


spearman_correlation, p = spearmanr(x, y)

r2 = r2_score(x, y)

print('spearman_correlation',spearman_correlation,p)
print('R2',r2)

correlation_coefficient = df.corr().iloc[0, 1]
print(correlation_coefficient)

plt.text(100, 100, f'Correlation: {correlation_coefficient:.2f} spearmanr: {spearman_correlation:.2f} r2_score: {r2:.2f}', fontsize=12, color='blue', bbox=dict(facecolor='white', alpha=0.5))

plt.show()


