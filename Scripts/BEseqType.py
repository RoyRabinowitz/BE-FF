#AUTHOR: Shiri Almog , shirialmog1@gmail.com
class BEseqType:
    def __init__(self,snpID,seq5,seq3, wt,mutation,readingFrame, aaPosition, geneName,geneID):
        self.snpID = snpID
        self.seq5=seq5
        self.seq3=seq3
        self.wt = wt
        self.mutation=mutation
        self.aaPosition = aaPosition
        self.readingFrame = readingFrame
        self.geneName=geneName
        self.geneID=geneID



    def __repr__(self):
        return "Id's:"+self.snpID+" , "\
               " \nsequence:"+self.seq5+ self.mutation + self.seq3+\
               " \nwt:"+self.wt

