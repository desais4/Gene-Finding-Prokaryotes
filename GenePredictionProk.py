# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 11:20:21 2020

@author: Shivani
"""

#1.PRE-PROCESSING
#open file containinng genome
MT=open("D:\MEng\SPRING\INTRO BIOINFOR\MTG.txt", "r")
contents=MT.read()      #store contents of file in an object
forw=[]
for i in range(0, len(contents)):
    if  contents[i]=='A': 
         forw.append('A')
    elif contents[i]=='T':
         forw.append('T')
    elif contents[i]=='G': 
          forw.append('G')
    elif contents[i]=='C':
          forw.append('C')

forward= ''.join([str(elem) for elem in forw])
reverse=contents[::-1]
comp=[]
#get complementary strand
for i in range(0, len(reverse)):
    if  reverse[i]=='A': 
        comp.append('T')
    elif reverse[i]=='T':
          comp.append('A')
    elif reverse[i]=='G': 
          comp.append('C')
    elif reverse[i]=='C':
          comp.append('G')
         
rev = ''.join([str(elem) for elem in comp]) 
#-----------------------------------------------------------------------------
#initialize objects
#arrays to store ORFs that store positions of start and stop codons 
ORF1start=[]
ORF1stop=[]
ORF2start=[]
ORF2stop=[]
ORF3start=[]
ORF3stop=[]
ORF1Cstart=[]
ORF1Cstop=[]

#arrays to organize start and stop codons in their respective reading frames
RF1start=[]
RF2start=[]
RF3start=[]
RF1stop=[]
RF2stop=[]
RF3stop=[]

RF1Cstart=[]
RF2Cstart=[]
RF3Cstart=[]
RF1Cstop=[]
RF2Cstop=[]
RF3Cstop=[]
#----------------------------------------------------------------------------
#2. FIND START AND STOP CODONS
#finding positions of start codons
startFor = [i for i in range(len(forward)) if forward.startswith("ATG", i) | forward.startswith("ATA", i)] 
startRev= [i for i in range(len(rev)) if rev.startswith("ATG", i) | rev.startswith("ATA", i)]  

#finding positions of stop codons
stopFor = [i for i in range(len(forward)) if forward.startswith("TAA", i) | forward.startswith("TAG", i)| forward.startswith("AGA", i)| forward.startswith("AGG", i)] 
stopRev=  [i for i in range(len(rev)) if rev.startswith("TAA", i) | rev.startswith("TAG", i)| rev.startswith("AGA", i)| rev.startswith("AGG", i)] 
#--------------------------------------------------------------------------------
#3. ASSIGN START AND STOP CODONS TO READING FRAMES:
#assigning startand stop codons to their respective reading frame
def findRF(startFor,stopFor, startRev, stopRev):
    for i in startFor:
        if i%3==0:
            RF1start.append(i)
        elif i%3==1:
            RF2start.append(i)
        elif i%3==2:
            RF3start.append(i)
            
    #assigning stop codons to their respective reading frame   
    for i in stopFor:
        if i%3==0:
            RF1stop.append(i)
        elif i%3==1:
            RF2stop.append(i)
        elif i%3==2:
            RF3stop.append(i)
            
    for i in startRev:
        if i%3==0:
            RF1Cstart.append(i)
        elif i%3==1:
            RF2Cstart.append(i)
        elif i%3==2:
            RF3Cstart.append(i)
            
    #assigning stop codons to their respective reading frame   
    for i in stopRev:
        if i%3==0:
            RF1Cstop.append(i)
        elif i%3==1:
            RF2Cstop.append(i)
        elif i%3==2:
            RF3Cstop.append(i)
#-----------------------------------------------------------------------------
#4. FIND ORFS IN ALL READING FRAMES
#call functions for all reading frames
def getORF(startframe,stopframe, ORFstart, ORFstop):
    i=0
    j=0
    while(i<len(startframe) and j<len(stopframe)): 
        if startframe[i]>stopframe[j] :    
            j+=1
        elif startframe[i]<stopframe[j] :
            ORFstart.append(startframe[i])
            ORFstop.append(stopframe[j])
            i+=1
#-----------------------------------------------------------------------------
#5. GET NON-OVERLAPPING ORFs
def getnonoverlapORFs(ORFstop, ORFstart):
    ORFstopn = list(set(ORFstop))
    ORFstopn.sort() 
    j=0
    for i in ORFstopn:
        latest = ORFstart[j]
        while(ORFstart[j] <= i):
            ORFstart[j] = latest
            if(j<len(ORFstart)-1):
                j +=1
            else:
                break
    ORFstartn = list(set(ORFstart))
    ORFstartn.sort()
    return (ORFstartn,ORFstopn)
#-----------------------------------------------------------------------------
#6. FILTERING ORFS BASED ON LENGTH AND GC CONTENT
def filterORFs(strand, ORFstartn,ORFstopn):
    for i, j in zip(ORFstartn,ORFstopn):
        if (j+3-i)>200:
            GC=strand[i:j+3]
            Ccont=GC.count('C')
            Gcont=GC.count('G')
            GCcontent=((Ccont+Gcont)/len(GC))*100
            if(GCcontent>30):
                print("PUTATIVE GENE START (bp):"+ str(i))
                print("PUTATIVE GENE END (bp):"+ str(j+3))
                print("Length="+ str(j+3-i))
                print("% GC Content="+ str(GCcontent))
                print(strand[i:j+3])
                print("---------------xxxxxx------------------")
#-----------------------------------------------------------------------------              
#7. FUNCTION CALLS FOR READING FRAMES 
# a.FORWARD STRAND
findRF(startFor, stopFor,startRev, stopRev) 
print("READING FRAME +1")           
getORF(RF1start, RF1stop, ORF1start,ORF1stop) 
ORF1startn,ORF1stopn=getnonoverlapORFs(ORF1stop,ORF1start)
filterORFs(forward, ORF1startn, ORF1stopn)

print("READING FRAME +2")           
getORF(RF2start, RF2stop, ORF2start,ORF2stop) 
ORF2startn,ORF2stopn=getnonoverlapORFs(ORF2stop,ORF2start)
filterORFs(forward, ORF2startn, ORF2stopn)

print("READING FRAME +3")           
getORF(RF3start, RF3stop, ORF3start,ORF3stop) 
ORF3startn,ORF3stopn=getnonoverlapORFs(ORF3stop,ORF3start)
filterORFs(forward, ORF3startn, ORF3stopn)

#b.COMPLEMENTARY STRAND
findRF(startFor, stopFor,startRev, stopRev) 
print("READING FRAME -1")           
getORF(RF1Cstart, RF1Cstop,  ORF1Cstart,ORF1Cstop) 
ORFstartn,ORFstopn=getnonoverlapORFs(ORF1Cstop,ORF1Cstart)
filterORFs(rev, ORFstartn, ORFstopn)
#-----------------------------------------------------------------------------
#8.ACCURACY OF PREDICTION

total=0

#Actual values computed as a sum of length + Percent GC content for each of the 13 genes
Actual=[1003.7,1084.99,1588.24,730.2,246.61,725.2,830.56,386.17,340.1,1422.27,1856.92,567.67,1187.28]

#Predicted/Epected values computed as a sum of length + Percent GC content for each of the 13 genes
Predicted=[1004.64,981.35, 1588.23,730.19 ,246.61 , 725.19, 896.83,747.83 ,343.46 , 1468.78, 1865.75,570.42 ,1206.99 ]

#assign score based on difference between actual and predicted value, more the difference, less the score
for i, j in zip(Actual,Predicted):
    if abs(i-j)<5:
        score=4
    elif abs(i-j)<10 and abs(i-j)>=5:
        score=3
    elif abs(i-j)>15 and abs(i-j)<=20:
        score=2
    else:
        score=1
    total+=score

#compute accuracy
print("Score earned (out of "+ str(4*len(Actual))+ ")="+ str(total))
print("Accuracy(in %)=Score earned/Total score=")
print(total/(4*len(Actual))*100)
#------------------------------------------------------------------------------        
#close file
MT.close()
