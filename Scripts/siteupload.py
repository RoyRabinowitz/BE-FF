#AUTHOR: Shiri Almog , shirialmog1@gmail.com
from Scripts.BEseqType import *
from Bio.Seq import Seq
import csv
from Scripts.baseEditorsTable import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter
from Scripts.transverse import SpecialCleanMatch,origPro

def importSNPS(SNPfile):
    rslts=[]
    with open(SNPfile,"r") as csv_file:
        try:
            csv_read = csv.reader(csv_file, delimiter=",")
        except MemoryError:
            print("MemoryError, try splitting your file into smaller files")
            exit()

        for row in csv_read:
            if row[0]=="Sample ID/ name":
                continue
            snpID=row[0]
            if len(row[2])>1:
                continue

            mutation=row[2]
            wt=row[1][25]
            sequence5 = row[1][0:25]
            sequence3 = row[1][26:]
            readingFrame=row[3]
            aaPosition=0
            geneName=0
            geneID=0
            newSNP=BEseqType(snpID,sequence5, sequence3, wt, mutation, readingFrame, aaPosition, geneName, geneID)
            rslts.append(newSNP)
        rsltsDic={}
        num=0
        for snp in rslts:
            rsltsDic[snp.snpID]=snp.geneName
            num=num+1

        return rslts, rsltsDic

def matchBE(snp, BElist):
    matches_list = []
    matches = {}
    for BE in BElist:
        if BElist[BE][5] == "U":
            sequence=snp.seq3
        else:
            sequence=snp.seq5
            sequence=sequence[::-1] # get reverse

        PAM = BElist[BE][0]
        # find if PAM matches in correct place
        if len(sequence) > BElist[BE][2]:
            match = False
            start = BElist[BE][1] -1
            end = BElist[BE][2]-1
            window = (end-start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<=lenPAM):
                    if sequence[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        match = True
                        break
            if match is True and BE!= "eA3A-BE3":
                matches_list.append(BE)
            elif match is True and snp.seq5[-1]=='T':
                matches_list.append(BE)
        matches[snp.snpID] = matches_list
    return matches

def cleanMatch(snp, Matches, BElist,rev):
    cleanMatches={}
    clean_list=[]
    quietDic={}
    quiet_list=[]
    locations_dic={}
    for match in Matches:
        for BE in Matches[match]:
            if BElist[BE][5] == "U":
                seq3 = snp.seq3
                seq5=snp.seq5
            else:
                seq3 = snp.seq5[::-1]  # get reverse
                seq5 = snp.seq3[::-1]  # get reverse

            totalSeq1 = seq5 + snp.mutation + seq3
            totalSeq2 = seq5 + snp.wt + seq3
            printSeq5 = seq5  # for results
            printSeq3 = seq3  # for results
            diff = 0

            if rev == False:
                protein_seq = Seq(totalSeq2).translate()
            else:
                protein_seq = Seq(totalSeq2).reverse_complement().translate()
            origProtein = protein_seq

            PAM=BElist[BE][0]
            start=BElist[BE][1]-1
            end=BElist[BE][2]-1
            locations=[]
            #find locations of PAM. for each location, check if clean
            window = (end - start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                j=0
                while (j<=lenPAM):
                    if seq3[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        locations.append(i + start) #maybe j is not neccesary?
                        break
            locations_dic[BE] = locations

            for loc in locations_dic[BE]:
                protein_match = False
                locFromEnd = len(printSeq3) - loc
                loc = loc + len(printSeq5)
                activation_window = totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - end:len(
                    printSeq3) - locFromEnd + len(printSeq5) - start + 1]
                num = 0  # number of times the variant appears within the activation window
                max_num = 0
                new_AW = []
                for i in range(len(activation_window)):
                    if activation_window[i]==snp.mutation:
                        num=num+1
                        new_AW.append(snp.wt)
                    else:
                        new_AW.append(activation_window[i])
                new_AW=''.join(new_AW)
                if num>max_num:
                    max_num=num

                if rev==False:
                    finalSeq = Seq(totalSeq1[0:loc - end] + str(new_AW) + totalSeq1[loc - start + 1:])
                else:
                    beginningP = Seq(
                        totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end]).reverse_complement()
                    new_AW = Seq(new_AW).reverse_complement()
                    endP = Seq(
                        totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:]).reverse_complement()
                    finalSeq = endP + new_AW + beginningP
                protein_seq_new = finalSeq.translate()
                if max_num==1:
                    clean_list.append(BE)
                if protein_seq==protein_seq_new:
                    protein_match=True

                if protein_match==True:
                    quiet_list.append(BE)
        cleanMatches[snp.snpID]=clean_list
        quietDic[snp.snpID]=quiet_list
    return cleanMatches, quietDic,origProtein

def beginningCut(snp):
    if snp.readingFrame == "1":
        if len(snp.seq5) % 3 == 1:
            snp.seq5 = snp.seq5[1:]
        elif len(snp.seq5) % 3 == 2:
            snp.seq5 = snp.seq5[2:]
    elif snp.readingFrame == "2":
        if len(snp.seq5) % 3 == 0:
            snp.seq5 = snp.seq5[2:]
        elif len(snp.seq5) % 3 == 2:
            snp.seq5 = snp.seq5[1:]
    elif snp.readingFrame == "3":
        if len(snp.seq5) % 3 == 0:
            snp.seq5 = snp.seq5[1:]
        elif len(snp.seq5) % 3 == 1:
            snp.seq5 = snp.seq5[2:]

    if (len(snp.seq5)+len(snp.seq3)+1) % 3 == 1:
        snp.seq3 = snp.seq3[:-1]
    if (len(snp.seq5)+len(snp.seq3)+1) % 3 == 2:
        snp.seq3 = snp.seq3[:-2]
    return snp

def cutFromRF(snp,totalSeq1,totalSeq2,printSeq5,printSeq3):
    if snp.readingFrame == "1":
        if len(snp.seq5) % 3 == 1:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq2[1:]
            printSeq5 = printSeq5[1:]
        elif len(snp.seq5) % 3 == 2:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
    elif snp.readingFrame == "2":
        if len(snp.seq5) % 3 == 0:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
        elif len(snp.seq5) % 3 == 2:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq1[1:]
            printSeq5 = printSeq5[1:]
    elif snp.readingFrame == "3":
        if len(snp.seq5) % 3 == 0:
            totalSeq1 = totalSeq1[1:]
            totalSeq2 = totalSeq2[1:]
            printSeq5 = printSeq5[1:]
        elif len(snp.seq5) % 3 == 1:
            totalSeq1 = totalSeq1[2:]
            totalSeq2 = totalSeq2[2:]
            printSeq5 = printSeq5[2:]
    if len(totalSeq1) % 3 == 1:
        totalSeq1 = totalSeq1[:-1]
        totalSeq2 = totalSeq2[:-1]
        printSeq3 = printSeq3[:-1]
    if len(totalSeq1) % 3 == 2:
        totalSeq1 = totalSeq1[:-2]
        totalSeq2 = totalSeq2[:-2]
        printSeq3 = printSeq3[:-2]
    return totalSeq1,totalSeq2,printSeq5,printSeq3



def getRevComp(snp):
    len3=len(snp.seq3)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len3])
    seq3=str(rc_seq[len3+1:])
    new_snp=BEseqType(snp.snpID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame, snp.aaPosition,snp.geneName,snp.geneID)
    return new_snp

def checkRF(snp):
    #this function will check the other 2 bases in the reading frame to see whether fixing them may result in a quiet result
    if snp.readingFrame=="1":
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.snpID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 1, 1, 1, 1)
        first_seq5=snp.seq5+snp.mutation
        first_mutation=snp.seq3[0]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq3[1:]
        snp1=BEseqType(snp.snpID,first_seq5,first_seq3,first_wt,first_mutation,2,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation+snp.seq3[0]
        second_mutation = snp.seq3[1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[2:]
        snp2 = BEseqType(snp.snpID, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0, snp1,snp2

    elif snp.readingFrame=="2":
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.snpID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 2, 1, 1, 1)
        first_seq5=snp.seq5[:-1]
        first_mutation=snp.seq5[-1]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.mutation+snp.seq3
        snp1=BEseqType(snp.snpID,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation
        second_mutation = snp.seq3[0]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[1:]
        snp2 = BEseqType(snp.snpID, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0,snp1,snp2

    # elif snp.readingFrame==3:
    else:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(snp.snpID, zero_seq5, zero_seq3, zero_wt, zero_mutation, 3, 1, 1, 1)
        first_seq5=snp.seq5[:-2]
        first_mutation=snp.seq5[-2]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq5[-1]+snp.mutation+snp.seq3
        snp1=BEseqType(snp.snpID,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5[:-1]
        second_mutation = snp.seq5[-1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.mutation+snp.seq3
        snp2 = BEseqType(snp.snpID, second_seq5, second_seq3,second_wt, second_mutation,2,1,1,1)
        return snp0, snp1,snp2

def find_cor(base):
    if base=="C":
        return "T"
    if base=="T":
        return "C"
    if base=="A":
        return "G"
    if base=="G":
        return "A"

def Mainsiteupload(DB):
    SNPS,rsltsDic=importSNPS(DB) #parsing cvs file
    matches={}
    cleanMatchdic = {}
    quietMatchdic={}
    # sort DNA. First, determine which bases we wish to replace. 4 cases:
    # 1. C to T: use CBE list       2. A to G: use ABE list
    # 3. T to C: switch to reverse complement and use ABE
    # 4. G to A: switch to RC and use CBE
    for snp in SNPS:

        rev = False
        snp = beginningCut(snp)
        if snp.mutation == "C" and snp.wt == "T":
            BElist = CBElist
            MinorBElist=CBElistMinor
        elif snp.mutation == "A" and snp.wt == "G":
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "T" and snp.wt == "C":
            snp=getRevComp(snp)
            rev = True
            snp.mutation = "A"
            snp.wt="G"
            BElist = ABElist
            MinorBElist = ABElistMinor
        elif snp.mutation == "G" and snp.wt == "A":
            snp = getRevComp(snp)
            rev = True
            snp.mutation = "C"
            snp.wt="T"
            BElist = CBElist
            MinorBElist = CBElistMinor
        else:
            BElist=None
            MinorBElist=None

        try:
          # check for matches in major window
            check_match = matchBE(snp, BElist)
            matches.update(check_match)
        except:
            pass
        try:
            snp0, snp2, snp3 = checkRF(snp)
            rev0,rev2,rev3=False,False,False
            if snp0.mutation == "C":
                BElist0 = CBElist
            elif snp0.mutation == "A":
                BElist0 = ABElist
            elif snp0.mutation == "T":
                snp0 = getRevComp(snp0)
                rev0 = True
                snp0.mutation = "A"
                snp0.wt = "G"
                BElist0 = ABElist
            elif snp0.mutation == "G":
                snp0 = getRevComp(snp0)
                rev0 = True
                snp0.mutation = "C"
                snp0.wt = "T"
                BElist0 = CBElist
            if snp2.mutation == "C":
                BElist2 = CBElist
            elif snp2.mutation == "A":
                BElist2 = ABElist
            elif snp2.mutation == "T":
                snp2 = getRevComp(snp2)
                rev2 = True
                snp2.mutation = "A"
                snp2.wt = "G"
                BElist2 = ABElist
            elif snp2.mutation == "G":
                snp2 = getRevComp(snp2)
                rev2 = True
                snp2.mutation = "C"
                snp2.wt = "T"
                BElist2 = CBElist

            if snp3.mutation == "C":
                BElist3 = CBElist
            elif snp3.mutation == "A":
                BElist3 = ABElist
            elif snp3.mutation == "T":
                snp3 = getRevComp(snp3)
                rev3 = True
                snp3.mutation = "A"
                snp3.wt = "G"
                BElist3 = ABElist
            elif snp3.mutation == "G":
                snp3 = getRevComp(snp3)
                rev3 = True
                snp3.mutation = "C"
                snp3.wt = "T"
                BElist3 = CBElist
        except:
            pass
        try:
            check_match0 = matchBE(snp0, BElist0)
            matches.update(check_match0)
        except:
            pass
        try:
            check_match2 = matchBE(snp2, BElist2)
            matches.update(check_match2)
        except:
            pass
        try:
            check_match3 = matchBE(snp3, BElist3)
            matches.update(check_match3)
        except:
            pass

            #check for matches in minor window
            #check_match_minor = matchBE(snp, MinorBElist)
            #matchesMinor.update(check_match_minor)

        #check for clean match
        try:
            clean, quiet,orig_protein= cleanMatch(snp,check_match,BElist,rev)

        except:
            clean={}
            quiet={}
            orig_protein=origPro(snp,rev)
        try:

            clean0,quiet0=SpecialCleanMatch(snp0,check_match0,BElist0,rev0,orig_protein)

        except:
            quiet0=[]
        try:
            clean2, quiet2 = SpecialCleanMatch(snp2, check_match2, BElist2,rev2,orig_protein)
        except:
            quiet2 = []
        try:
            clean3, quiet3 = SpecialCleanMatch(snp3, check_match3, BElist3,rev3,orig_protein)
        except:
            quiet3 = []

        for key in quiet0:
            if snp.snpID not in quiet:
                quiet[snp.snpID]=[key]
            elif key not in quiet[snp.snpID]:
                quiet[snp.snpID].append(key)
        for key in quiet2:
            if snp.snpID not in quiet:
                quiet[snp.snpID]=[key]
            elif key not in quiet[snp.snpID]:
                quiet[snp.snpID].append(key)
        for key in quiet3:
            if snp.snpID not in quiet:
                quiet[snp.snpID]=[key]
            elif key not in quiet[snp.snpID]:
                quiet[snp.snpID].append(key)
        try:
            quietMatchdic.update(quiet)
            cleanMatchdic.update(clean)

        except:
            pass

    return matches, cleanMatchdic, quietMatchdic, rsltsDic



