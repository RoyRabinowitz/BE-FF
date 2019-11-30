import seaborn as sns; sns.set()
import numpy as np
import csv
import matplotlib
import matplotlib.pyplot as plt


file1="c:/Python/BE AA heatmap/all_transversions.csv" #all SNPs
file2="c:/Python/BE AA heatmap/transversions_results.csv" #results file

AA=['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','W','Y','V','*']

ALL=np.zeros(441).reshape(21,21)
results=np.zeros(441).reshape(21,21)
ratio=np.zeros(441).reshape(21,21)
rows_in_ALL=0
rows_in_results=0


with open(file1, 'r') as f1:
    for line in f1:  #f1 is the All SNPs file
        rows_in_ALL+=1
        line=line.split(',')
        aa=line[8]
        
        aa=aa.split('>')
        try:
            Ref_aa=aa[0]
            var_aa=aa[1]
        #print(rows_in_ALL)
            for x in range(len(AA)):
                if AA[x]==Ref_aa:
                    for y in range(len(AA)):
                        if AA[y]==var_aa:
                            ALL[x][y]+=1
        except:
            pass

with open(file2, 'r') as f2:
    for line in f2:  #f2 is the results file
        rows_in_results+=1
        line=line.split(',')
        aa=line[8]
        aa=aa.split('>')
        Ref_aa=aa[0]
        var_aa=aa[1]
        for x in range(len(AA)):
            if AA[x]==Ref_aa:
                for y in range(len(AA)):
                    if AA[y]==var_aa:
                        results[x][y]+=1

for x in range(len(AA)):
    for y in range(len(AA)):
        ratio[x][y]=results[x][y]/ALL[x][y]

cmap = matplotlib.colors.ListedColormap(["#FFFFFF"] + sns.color_palette('Reds', 100)[1:])
ax = sns.heatmap(ratio,vmin=0, vmax=1, cmap=cmap, xticklabels=AA,yticklabels=AA, square=True, linewidths=1, linecolor='whitesmoke',)
ax.invert_yaxis()
plt.xlabel('Variant AA')
plt.ylabel('Reference AA')
plt.savefig('transversions.png', dpi=300)




