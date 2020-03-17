#AUTHOR: Shiri Almog , shirialmog1@gmail.com
from Scripts.BEseqType import *
from Bio.Seq import Seq
from Scripts.baseEditors import CBElist, CBElistMinor, ABElist, ABElistMinor, BEletter
from Scripts.transverse_single import SpecialCleanMatch,origPro

## This function returns all the BE who's PAM'S are found in the given sequence
def matchBE(snp, BElist):
    matches_list = []
    for BE in BElist:
        if BElist[BE][5] == "U": ## U/D for upstream or downstream PAM, relative to mutation
            sequence=snp.seq3
        else:
            sequence=snp.seq5
            sequence=sequence[::-1] # get reverse

        PAM = BElist[BE][0]
        if len(sequence) > BElist[BE][2]:   # find if PAM matches in correct place
            match = False
            start = BElist[BE][1] -1
            end = BElist[BE][2]-1
            window = (end-start+1)
            for i in range(window):
                temp = start+i
                lenPAM= len(PAM)
                if temp+lenPAM>len(sequence):
                    break
                j=0

                while (j<lenPAM ):
                    if sequence[temp] in BEletter[PAM[j]]:
                        j+=1
                        temp+=1
                    else:
                        break
                    if j == len(PAM):
                        match = True
                        break

            if match is True and BE!= "eA3A-BE3": ##unique case, only works if T comes before mutation
                matches_list.append(BE)
            elif match is True and snp.seq5[-1]=='T': #BE== "eA3A-BE3"
                matches_list.append(BE)
    return matches_list

## This function takes as input the found matches_list, and checks for precise corrections
def cleanMatch(snp,Matches, BElist,rev):
    origMutSeq = ''
    clean_dic={}
    quiet_dic={}
    locations_dic = {}
    originalProtein=origPro(snp,rev) ##saves the original AA sequence, for comparison later
    for BE in Matches:
        if BElist[BE][5] == "U":
            seq3 = snp.seq3
            seq5=snp.seq5
        else:
            seq3 = snp.seq5[::-1] # get reverse
            seq5 = snp.seq3[::-1]  # get reverse

        totalSeq1 = seq5 + snp.mutation + seq3
        totalSeq2 = seq5 + snp.wt + seq3
        printSeq5 = seq5  # for results
        printSeq3 = seq3  # for results
        if rev == False:
            protein_seq = Seq(totalSeq2).translate()
        else:
            protein_seq = Seq(totalSeq2).reverse_complement().translate()

        PAM=BElist[BE][0]
        start=BElist[BE][1]-1
        end=BElist[BE][2]-1
        locations_list=[]

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
                    locations_list.append(i + start)
                    break
        locations_dic[BE]=locations_list
        clean_list = []
        quiet_list=[]
        origMutSeq = {}
        printPam = {}
        printRevCorSeq = {}
        finalShowSeq = {}
        off_target={}
        for loc in locations_dic[BE]:
            protein_match = False
            locFromEnd=len(printSeq3)-loc
            loc=loc+len(printSeq5)
            activation_window=totalSeq1[len(printSeq3)-locFromEnd+len(printSeq5)-end:len(printSeq3)-locFromEnd+len(printSeq5)-start+1]
            num = 0  # number of times the variant appears within the activation window
            max_num=0
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
            ##This is for the visual output
            a=printSeq5+snp.mutation+printSeq3
            off_target[loc]=a[loc-19:loc+1]
            if rev==False:
                new = new_AW
                finalSeq = Seq(totalSeq1[0:loc - end] + str(new) + totalSeq1[loc - start+1:])
                beginningP = Seq(
                    totalSeq1[0:len(printSeq3) - locFromEnd + len(printSeq5) - end])
                endP = Seq(totalSeq1[len(printSeq3) - locFromEnd + len(printSeq5) - start + 1:])
                if BElist[BE][5]=="U":
                    origMutSeq[loc] = printSeq5 + "<b>"+snp.mutation+"</b>" + printSeq3 #sequence with bold mutation
                    printPam[loc] = printPamSeq(beginningP, new_AW, endP, locFromEnd, PAM,BElist[BE][5])
                else:
                    origMutSeq[
                        loc] = printSeq3[::-1] + "<b>" + snp.mutation + "</b>" + printSeq5[::-1]  # sequence with bold mutation
                    printPam[loc] = printPamSeq(beginningP, new_AW, endP, locFromEnd, PAM, BElist[BE][5])

            else:
                beginningP=Seq(totalSeq1[0:len(printSeq3)-locFromEnd+len(printSeq5)-end]).reverse_complement()
                old_AW=Seq(activation_window).reverse_complement()
                new_AW=Seq(new_AW).reverse_complement()
                endP=Seq(totalSeq1[len(printSeq3)-locFromEnd+len(printSeq5)-start+1:]).reverse_complement()
                finalSeq = endP + new_AW + beginningP
                if BElist[BE][5] == "U":
                    finalShowSeq[loc] = endP + "<class style='color:blue'><b>" + new_AW + "</b></class>" + beginningP
                    origMutSeq[loc] = Seq(printSeq3).reverse_complement() + "<b>" + Seq(snp.mutation).reverse_complement() + "</b>" + Seq(printSeq5).reverse_complement()
                    printRevCorSeq[loc]= printPamSeq(beginningP.reverse_complement(), new_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM,BElist[BE][5])
                    printPam[loc]= printPamSeq(beginningP.reverse_complement(), old_AW.reverse_complement(), endP.reverse_complement(), locFromEnd, PAM,BElist[BE][5])
                else:
                    finalShowSeq[loc] = beginningP[::-1] + "<class style='color:blue'><b>" + new_AW[::-1] + "</b></class>" + endP[::-1]
                    origMutSeq[loc] = Seq(printSeq5).reverse_complement()[::-1] + "<b>" + Seq(
                        snp.mutation).reverse_complement() + "</b>" + Seq(printSeq3).reverse_complement()[::-1]
                    printRevCorSeq[loc] = printPamSeq(beginningP.reverse_complement(), new_AW.reverse_complement(),
                                                      endP.reverse_complement(), locFromEnd, PAM, BElist[BE][5])
                    printPam[loc] = printPamSeq(beginningP.reverse_complement(), old_AW.reverse_complement(),
                                                endP.reverse_complement(), locFromEnd, PAM, BElist[BE][5])

            protein_seq_new=finalSeq.translate()
            if protein_seq==protein_seq_new:
                protein_match=True
            if max_num==1:
                clean_list.append(loc)
            elif protein_match == True:
                quiet_list.append(loc)

        clean_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM,clean_list,rev,off_target]
        quiet_dic[BE]=[origMutSeq,printPam,printRevCorSeq,finalShowSeq,PAM, quiet_list,rev,off_target]

    return clean_dic,quiet_dic,origMutSeq,locations_dic,originalProtein
##For visual output
def printPamSeq(seq5,activationWindow,seq3, locfromEnd,PAM,direction):
    if direction=='U':
        len3=len(seq3)
        printPAM=seq5+"<class style='color:blue'><b>"+activationWindow+"</b></class>"+seq3[0:len3-locfromEnd]+"<b><class style='color:#DD96F0'>"+seq3[len3-locfromEnd:len3-locfromEnd+len(PAM)]+"</b></class>"+seq3[len3-locfromEnd+len(PAM):]
        return printPAM
    elif direction=='D':
        len3 = len(seq3)
        printPAM=seq3[len3 - locfromEnd + len(PAM):][::-1]+ "<b><class style='color:#DD96F0'>" + seq3[len3 - locfromEnd:len3 - locfromEnd + len( PAM)] [::-1]+"</b></class>"+seq3[0:len3 - locfromEnd][::-1]+ "<class style='color:blue'><b>" + activationWindow[::-1] + "</b></class>"+seq5[::-1]
        return printPAM

def getRevComp(snp):
    len5=len(snp.seq5)
    len3 = len(snp.seq3)
    totalSeq= snp.seq5+snp.mutation+snp.seq3
    total_seq=Seq(totalSeq)
    rc_seq=total_seq.reverse_complement()
    seq5=str(rc_seq[0:len3])
    start=len(rc_seq)-len5
    seq3=str(rc_seq[start:])
    new_snp=BEseqType(snp.snpID,seq5 ,seq3,snp.mutation,snp.wt, snp.readingFrame, snp.aaPosition,snp.geneName,snp.geneID)
    return new_snp

#this function will check the other 2 bases in the reading frame to see whether fixing them may result in a synonymous correction
def checkRF(snp):
    if snp.readingFrame=='1':
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 1, 1, 1, 1)
        first_seq5=snp.seq5+snp.mutation
        first_mutation=snp.seq3[0]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq3[1:]
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,2,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation+snp.seq3[0]
        second_mutation = snp.seq3[1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[2:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0, snp1,snp2

    elif snp.readingFrame=='2':
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 2, 1, 1, 1)
        first_seq5=snp.seq5[:-1]
        first_mutation=snp.seq5[-1]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5 + snp.mutation
        second_mutation = snp.seq3[0]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.seq3[1:]
        snp2 = BEseqType(12, second_seq5, second_seq3, second_wt,second_mutation, 3,1,1,1)
        return snp0,snp1,snp2

    # elif snp.readingFrame==3:
    else:
        zero_seq5 = snp.seq5
        zero_mutation = snp.mutation
        zero_wt = find_cor(zero_mutation)
        zero_seq3 = snp.seq3
        snp0 = BEseqType(12, zero_seq5, zero_seq3, zero_wt, zero_mutation, 3, 1, 1, 1)
        first_seq5=snp.seq5[:-2]
        first_mutation=snp.seq5[-2]
        first_wt=find_cor(first_mutation)
        first_seq3=snp.seq5[-1]+snp.mutation+snp.seq3
        snp1=BEseqType(12,first_seq5,first_seq3,first_wt,first_mutation,1,1,1,1)
        second_seq5 = snp.seq5[:-1]
        second_mutation = snp.seq5[-1]
        second_wt = find_cor(second_mutation)
        second_seq3 = snp.mutation+snp.seq3
        snp2 = BEseqType(12, second_seq5, second_seq3,second_wt, second_mutation,2,1,1,1)
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

def delete_chars(field):
    field_new=field.replace(" ","")
    for char in field_new:
        if char!='A' and char!='G' and char!='T' and char!='C':
            field_new=field_new.replace(char,"")
    return field_new


def MainBE(upSeq, downSeq, mutation, wt, readingFrame,personalPAM,start,end,fromto,stream):
    upSeq_=delete_chars(upSeq.upper())
    downSeq_=delete_chars(downSeq.upper())
    wt_=delete_chars(wt.upper())
    mutation_=delete_chars(mutation.upper())
    snp=BEseqType(1234,upSeq_, downSeq_, wt_, mutation_, readingFrame, 1, 12, 12)
    refSeq=snp.seq5+"<b>"+snp.wt+"</b>"+snp.seq3
    mutSeq = snp.seq5 + "<b>" + snp.mutation + "</b>" + snp.seq3
    ## initiating vars
    minor_clean_list=[]
    minor_quiet_list=[]
    origMutSeq = {}
    locations_dic = {}
    rev,rev0,rev2,rev3=False,False,False,False
    snp=beginningCut(snp)
    if snp.mutation == "C" and snp.wt == "T":
        BElist = CBElist
        MinorBElist=CBElistMinor
        if fromto == '1' and personalPAM!="" and start!="" and end!="":
            BElist.update({'User customized BE': [personalPAM, start, end, 'C', 'T', stream]})
    elif snp.mutation == "A" and snp.wt == "G":
        BElist = ABElist
        MinorBElist = ABElistMinor
        if fromto == '2' and personalPAM!="" and start!="" and end!="":
            BElist.update({'User customized BE': [personalPAM, start, end, 'A', 'G', stream]})
    elif snp.mutation == "T" and snp.wt == "C":
        snp=getRevComp(snp)
        rev=True
        snp.mutation = "A"
        snp.wt="G"
        BElist = ABElist
        MinorBElist = ABElistMinor
    elif snp.mutation == "G" and snp.wt == "A":
        snp = getRevComp(snp)
        rev=True
        snp.mutation = "C"
        snp.wt="T"
        BElist = CBElist
        MinorBElist = CBElistMinor
    else:
        BElist={}
        MinorBElist = None

    snp0,snp2,snp3=checkRF(snp)
    rev0, rev2, rev3 = False, False, False
    BElist0,BElist2,BElist3={},{},{}
    if snp0.mutation == "C":
        BElist0 = CBElist
        if fromto == '1' and personalPAM!="" and start!="" and end!="":
            BElist0.update({'User customized BE': [personalPAM, start, end, 'C', 'T', stream]})
    elif snp0.mutation == "A":
        BElist0 = ABElist
        if fromto == '2' and personalPAM!="" and start!="" and end!="":
            BElist0.update({'User customized BE': [personalPAM, start, end, 'A', 'G', stream]})
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
        if fromto == '1' and personalPAM!="" and start!="" and end!="":
            BElist2.update({'User customized BE': [personalPAM, start, end, 'C', 'T', stream]})
    elif snp2.mutation == "A":
        BElist2 = ABElist
        if fromto == '2' and personalPAM!="" and start!="" and end!="":
            BElist2.update({'User customized BE': [personalPAM, start, end, 'A', 'G', stream]})
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
        if fromto == '1' and personalPAM!="" and start!="" and end!="":
            BElist3.update({'User customized BE': [personalPAM, start, end, 'C', 'T', stream]})
    elif snp3.mutation == "A":
        BElist3 = ABElist
        if fromto == '2' and personalPAM!="" and start!="" and end!="":
            BElist3.update({'User customized BE': [personalPAM, start, end, 'A', 'G', stream]})
    elif snp3.mutation =="T":
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

    try:
        check_match = matchBE(snp, BElist)
    except:
        pass
    try:
        check_match_minor = matchBE(snp, MinorBElist)
    except:
        pass
    try:
        check_match0 = matchBE(snp0, BElist0)
    except:
        pass
    try:
        check_match2=matchBE(snp2,BElist2)
    except:
        pass
    try:
        check_match3=matchBE(snp3,BElist3)
    except:
        pass
    #check for clean match
    try:
        clean_dic,quiet_dic,origMutSeq, locations_dic,orig_protein= cleanMatch(snp, check_match, BElist,rev)
    except:
        clean_dic = {}
        quiet_dic = {}
        orig_protein = origPro(snp, rev)

    try:
        clean_dic_m, quiet_dic_m, origMutSeq_m, locations_dic_m, orig_protein_m = cleanMatch(snp, check_match_minor, MinorBElist, rev)
        for BE in clean_dic_m:
            if BE in clean_dic:
                minor_clean_list.append(BE)
        for BE in clean_dic_m:
            if BE in clean_dic:
                minor_quiet_list.append(BE)
    except:
        pass
    try:
        clean_dic0,quiet_dic0,origMutSeq0,locations_dic0 = SpecialCleanMatch(snp0, check_match0, BElist0, rev0, orig_protein)
    except:
        quiet_dic0={}
    try:
        clean_dic2,quiet_dic2,origMutSeq2,locations_dic2 = SpecialCleanMatch(snp2, check_match2, BElist2, rev2, orig_protein)
    except:
        quiet_dic2={}
    try:
        clean_dic3,quiet_dic3,origMutSeq3,locations_dic3 = SpecialCleanMatch(snp3, check_match3, BElist3, rev3, orig_protein)
    except:
        quiet_dic3={}
    ## FOR VISUAL OUTPUT
    syn_quiet={}
    for key in quiet_dic0:
        if key in quiet_dic:
            for loc in quiet_dic0[key][5]:
                if loc not in quiet_dic[key][5]:
                    syn_quiet.update({key:quiet_dic0[key]})
        else:
            syn_quiet.update({key: quiet_dic0[key]})
    for key in quiet_dic2:
        if key in quiet_dic:
            for loc in quiet_dic2[key][5]:
                if loc not in quiet_dic[key][5]:
                    syn_quiet.update({key: quiet_dic2[key]})
        elif key not in syn_quiet:
            syn_quiet.update({key: quiet_dic2[key]})
    for key in quiet_dic3:
        if key in quiet_dic:
            for loc in quiet_dic3[key][5]:
                if loc not in quiet_dic[key][5]:
                    syn_quiet.update({key:quiet_dic3[key]})
        elif key not in syn_quiet:
            syn_quiet.update({key: quiet_dic3[key]})
    try:
        return clean_dic, quiet_dic,refSeq,mutSeq, origMutSeq,locations_dic,syn_quiet,minor_clean_list,minor_quiet_list
    except:
        return [],[],[],[],[],[],[],[]


