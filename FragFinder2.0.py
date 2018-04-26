# -*- coding: utf-8 -*-
"""
@author: Nicolaas Hermans
nhermans@leiden.physics.nl
"""

#finds specefic (restriction) sequence in Genome, outputs n basepairs sequence 
#around the restriction site. Note, 5' and 3' overhangs need to be defined by 
#the user. For now, only palindromic sites can be used.

restrict1 = 'CCATGG' #define restriction site (palindrome, if not palindrome, do sequences independently).
restrict2 = 'AGATCT'
#restrict3 = 'CTCGAG'
RestrictList = [restrict1,restrict2]
folder = 'Human' #folder with chromosome sequence files (note, do not put other files in this folder)
#probe = 'GCACATATACACCATGG' #define site (palindrome, if not palindrome, do sequences independently).
probe = 'TCCATTGGCTTGGTAGATC' #BglII probe.
Max=8000 #max size of fragments
Min=1000 #min size of fragment
Width=1 #width of the outputbins for the historgram in bp

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

rprobe=''        
for x in probe:
    rprobe=comp_string(x)+rprobe
print probe, rprobe
start = time.time()
filenames = os.listdir(folder)
os.chdir( folder )

totbp=0
count=0
print RestrictList

fragments = []    
for filename in filenames:
    #Import genome as a single string
    my_file = open(filename)    
    file_contents = my_file.read()    
    #sequence = file_contents.strip('\n') #this doesn't work
    sequence = file_contents.replace("\n", "")
    #sequence = re.sub(r'\d+|\W+', '', sequence) #removes numbers etc. makes stuff slower for large sequence files
    #sequence = sequence.upper()      #changes everything to upper case, idem 
    allindices = []    
    for s in RestrictList:
        cutslist= re.finditer(s, sequence) 
        indices = [m.start(0) for m in cutslist] 
        allindices = allindices+indices
    cuts=sorted(allindices)
    
    probelist= re.finditer(probe, sequence)
    ProbeInd=  [m.start(0) for m in probelist] 
    
    rv=len(probe)
    for p in ProbeInd:
        oldtest=Max
        for c in cuts:
            test=p-c
            if test>0 and test>Min and test<oldtest:  #cycles through all cutsites untill the nearest neighbour
                oldtest=test #stop the forloop next instance
                FragmentLength=test+rv               
                fragments.append(FragmentLength) #add length of the 3' probe to the length of the fragment
                count=count+1    
                   
    probelist= re.finditer(rprobe, sequence)
    ProbeInd=  [m.start(0) for m in probelist] 

    for p in ProbeInd:
        oldtest=Max
        for c in cuts:
            test=c-p
            if test>0 and test>Min and test<oldtest:  #cycles through all cutsites untill the nearest neighbour
                oldtest=test #stop the forloop next instance
                fragments.append((test+6)) #add length of the 3' probe to the length of the fragment
                count=count+1    
    bp = len(sequence)
    fragments = sorted(fragments)
    #print 'sequence length =',bp ,'bp'
    print 'total fragments analyzed:', count,'in file  ',filename,'(', bp,'bp)'
    totbp=totbp+bp
    my_file.close()
    re.purge()

end = time.time()
print 'time =', end - start
print 'total', totbp, 'bp, and', count,'fragments'  

Bins=(Max-Min)/Width
hist,bin_edges=np.histogram(fragments,bins=Bins, range=(Min,Max))
f = open('Pulldown_Fragments.txt', 'w')
Header=restrict1+' '+restrict2+'\n'
f.write(Header)
for k in hist:
    f.write(str(Min))    
    f.write(" %s\n" % str(k))
    Min=Min+Width
f.close()

plt.hist(fragments, bins=Bins)
plt.title("fragments sizes")
plt.xlabel("Fragments sizes (bp)")
plt.ylabel("Frequency")

fig = plt.gcf()