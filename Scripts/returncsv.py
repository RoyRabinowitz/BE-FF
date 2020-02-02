#AUTHOR: Shiri Almog , shirialmog1@gmail.com
import csv

def returncsv(matches, cleanMatchdic, quietMatchdic,rsltsDic,filepath):
    with open(filepath, mode='w',newline='') as f:
        file_writer=csv.writer(f, delimiter=',',quoting=csv.QUOTE_MINIMAL)
        file_writer.writerow(["snpID","BE1", "BE2", "BE3", "HF-BE3", "BE4(max)", "BE4-Gam","YE1-BE3","YEE-BE3", "VQR-BE3","VRER-BE3","SaBE3", "SaBE4", "SaBE4-Gam", "Sa(KKH)-BE3","Cas12a-BE","Target-AID","Target-AID-NG","xBE3","eA3A-BE3","BE-PLUS","CP-CBEmax variants","evoAPOBEC1-BE4max", "evoFERNY-BE4max","evoCDA1-BE4max", "ABE 7.9","ABE 7.10","ABE 7.10*","xABE","NG-ABEmax" ,"ABESa","VQR-ABE","VRER-ABE","Sa(KKH)-ABE","CP-ABEmax variants"])
        keyList = matches.keys()
        beList=["BE1", "BE2", "BE3", "HF-BE3", "BE4(max)", "BE4-Gam","YE1-BE3","YEE-BE3", "VQR-BE3","VRER-BE3","SaBE3", "SaBE4", "SaBE4-Gam", "Sa(KKH)-BE3","Cas12a-BE","Target-AID","Target-AID-NG","xBE3","eA3A-BE3","BE-PLUS","CP-CBEmax variants","evoAPOBEC1-BE4max", "evoFERNY-BE4max","evoCDA1-BE4max", "ABE 7.9","ABE 7.10","ABE 7.10*","xABE","NG-ABEmax" ,"ABESa","VQR-ABE","VRER-ABE","Sa(KKH)-ABE","CP-ABEmax variants"]
        for key in keyList:
            mlist=[key]
            for BE in beList:
                temp=''
                if key in quietMatchdic:
                    if BE in quietMatchdic[key]:
                        temp='Quiet'
                if key in cleanMatchdic:
                    if BE in cleanMatchdic[key]:
                        temp='Clean'
                mlist.append(temp)

            file_writer.writerow(mlist)

