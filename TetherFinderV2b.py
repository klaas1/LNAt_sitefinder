# -*- coding: utf-8 -*-
"""
@author: Nicolaas Hermans
nhermans@leiden.physics.nl
"""

#finds specefic (restriction) sequence in Genome, outputs n basepairs sequence 
#around the restriction site. Note, 5' and 3' overhangs need to be defined by 
#the user. For now, only palindromic sites can be used.

folder = 'Human' #folder with chromosome sequence files (note, do not put other files in this folder)
probe1 = 'GCACATATACACCATGG' #define site (palindrome, if not palindrome, do sequences independently).
probe2 = 'TCCATTTGCTTGGTAGATC' #t+ccat+tt+gct+tgg+ta+gatc BglII, note: check strands for tethering!
restrict1 = 'CCATGG' #define restriction site (palindrome, if not palindrome, do sequences independently).
restrict2 = 'AGATCT'

Min=1000 #min size of fragment
Max=8000 #max size of fragments
Width=1 #width of the outputbins in bp

import re #regular expressions lib
import os #filenames
import time # just a timer
import matplotlib.pyplot as plt #graphs
import numpy as np

def comp_string(base):
    if(base == "T"):
        return "A" 
    if(base == "A"):
        return "T" 
    if(base == "G"):
        return "C" 
    if(base == "C"):
        return "G"  

rprobe1=''
rprobe2=''        
for x in probe1:
    rprobe1=comp_string(x)+rprobe1
for x in probe2:
    rprobe2=comp_string(x)+rprobe2
print(rprobe1, rprobe2)

start = time.time()
filenames = os.listdir(folder)
os.chdir( folder )

totbp=0
count=0
count2=0
fragment=0

fragments = []  
FragSeq='' 
for filename in filenames:
    #Import genome as a single string
    my_file = open(filename)    
    file_contents = my_file.read()    
    #sequence = file_contents.strip('\n') #this doesn't work
    sequence = file_contents.replace("\n", "")
    #sequence = re.sub(r'\d+|\W+', '', sequence) #removes numbers etc. makes stuff slower for large sequence files
    #sequence = sequence.upper()      #changes everything to upper case, idem 
    
    cutslist= re.finditer(rprobe1, sequence) 
    indices = [m.start(0) for m in cutslist] 
    ProbeListSense=sorted(indices)
    cutslist= re.finditer(probe2, sequence) #use tether for oposing strand
    indices = [m.start(0) for m in cutslist]
    ProbeListAntisense=sorted(indices)
    rv=len(rprobe2)
    for p in ProbeListSense:
        oldtest=Max
        for c in ProbeListAntisense:
            test=c-p
            if test>0 and test>Min and test<oldtest:  #cycles through all cutsites untill the nearest neighbour
                oldtest=test #stop the forloop next instance
                x=p+1 #remove restriction sites in the probe from fragment
                y=c+rv-1
                a = re.search(restrict1, sequence[x:y]) 
                b = re.search(restrict2, sequence[x:y])
                if a==None and b==None:
                    FragmentLength=c-p+rv                    
                    fragments.append(FragmentLength) #add length of the 3' probe to the length of the fragment
                    FragSeq=FragSeq+sequence[x:y]+'\n'                    
                    count2=count2+1
                count=count+1
            
    cutslist= re.finditer(rprobe2, sequence) 
    indices = [m.start(0) for m in cutslist]
    ProbeListSense=sorted(indices)
    cutslist= re.finditer(probe1, sequence) 
    indices = [m.start(0) for m in cutslist] 
    ProbeListAntisense=sorted(indices)
    rv=len(probe1)
    for p in ProbeListSense:
        oldtest=Max
        for c in ProbeListAntisense:
            test=c-p
            if test>0 and test>Min and test<oldtest:  #cycles through all cutsites untill the nearest neighbour
                oldtest=test 
                x=p+1 
                y=c+rv-1
                a = re.search(restrict1, sequence[x:y]) 
                b = re.search(restrict2, sequence[x:y])
                if a==None and b==None:
                    FragmentLength=c-p+rv  
                    fragments.append(FragmentLength)
                    FragSeq=FragSeq+sequence[x:y]+'\n'                     
                    count2=count2+1
                count=count+1
   
    bp = len(sequence)
    fragments = sorted(fragments)
    #print 'sequence length =',bp ,'bp'
    print('total fragments analyzed:', count,'in file  ',filename,'(',bp,'bp)')
    totbp=totbp+bp
    my_file.close()
    re.purge()

end = time.time()
print('time =', end - start)
print('Genome was', totbp, 'bp, found', count,'possible fragments, resulting in' ,count2,'tethers')

Bins=(Max-Min)/Width
hist,bin_edges=np.histogram(fragments,bins=Bins, range=(Min,Max))
f = open('Tether_Fragments.txt', 'w')
Header=probe1+' '+rprobe2+' '+rprobe1+' '+probe2+'\n'
f.write(Header)
for k in hist:
    f.write(str(Min))    
    f.write(" %s\n" % str(k)) 
    Min=Min+Width
f.close()

f = open('Tether_Sequences.txt', 'w')
f.write(FragSeq)
f.close()

plt.hist(fragments, bins=Bins)
plt.title("fragments sizes")
plt.xlabel("Fragments sizes (bp)")
plt.ylabel("Frequency")

fig = plt.gcf()