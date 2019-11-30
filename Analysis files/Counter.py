
file="c:/Python/BE counter/allSNPs.csv" #all SNPs
results="c:/Python/BE counter/Results_CSV.csv"

#Counters
onlyperfect=0
onlysyn=0
both=0
lines=0

#transition counters
transitionperfect=0
transitionsyn=0
transitionboth=0

#transversion counters
transversionperfect=0
transversionsyn=0
transversionboth=0


lst=['A>G','A>C','A>T','C>T','C>G','C>A','G>A','G>C','G>T','T>A','T>C','T>G']
transitions=['A>G','C>T','G>A','T>C']
transversions=['A>C','A>T','C>G','C>A','G>C','G>T','T>A','T>G']
counter=[0,0,0,0,0,0,0,0,0,0,0,0]
matches=[0,0,0,0,0,0,0,0,0,0,0,0]



with open(file, 'r') as f1:
    for line in f1:
        line=line.split(',')
        for i in range(len(lst)):
            if line[7]==lst[i]:
                counter[i]+=1

with open(results, 'r') as f2:
    for line in f2:
        lines+=1
        matchlist=[0,0] #onlyperfect, onlysyn
        line=line.split(',')
        
        for l in range (len(line)):
            if line[l]=='Perfect':
                matchlist[0]+=1
            elif line[l]=='Synonymous':
                matchlist[1]+=1
        
        if matchlist[0]!=0 and matchlist[1]==0:
            onlyperfect+=1
            for a in transitions:
                if line[7]==a:
                    transitionperfect+=1
            for b in transversions:
                if line[7]==b:
                    transversionperfect+=1
        elif matchlist[0]==0 and matchlist[1]!=0:
            onlysyn+=1
            for a in transitions:
                if line[7]==a:
                    transitionsyn+=1
            for b in transversions:
                if line[7]==b:
                    transversionsyn+=1
        elif matchlist[0]!=0 and matchlist[1]!=0:
            both+=1
            for a in transitions:
                if line[7]==a:
                    transitionboth+=1
            for b in transversions:
                if line[7]==b:
                    transversionboth+=1
        
        for i in range(len(lst)):
            if line[7]==lst[i]:
                matches[i]+=1

        
total=0
totalresults=0

for i in range(len(lst)):
    total+=counter[i] 
    totalresults+=matches[i]

totalratio=totalresults/total

print ('Total SNPs checked: ',total)
print ('Total editable SNPs: ',totalresults)
print('This is ',totalratio*100, 'percent of the SNPs')

for i in range(len(lst)):
    print('For',lst[i],' SNPs, there were total of',counter[i],'SNP. ', matches[i],' are editable and the ratio is ',matches[i]/counter[i])





print('Total reads: ',lines)
print('In general:')
print('SNPs corrected by perfect editing: ',onlyperfect)
print('SNPs corrected by synonyms editing: ',onlysyn)
print('SNPs corrected by both perfect and synonyms editing: ',both)  
print('Transitions:')
print('Transition perfect corrections: ',transitionperfect)
print('Transition synonymous corrections: ',transitionsyn)
print('Transition both corrections: ',transitionboth)
print('Transversions:')
print('Transversions perfect corrections: ',transversionperfect)
print('Transversions synonymous corrections: ',transversionsyn)
print('Transversions both corrections: ',transversionboth)
