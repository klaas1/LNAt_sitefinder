# -*- coding: utf-8 -*-
"""
Created on Mon May 11 10:47:17 2015
@author: Nicolaas Hermans
nhermans@leiden.physics.nl
"""

#finds specefic (restriction) sequence in Genome, outputs n basepairs sequence 
#around the restriction site. Note, 5' and 3' overhangs need to be defined by 
#the user. For now, only palindromic sites can be used.

restrict = 'GGTAACC' #define restriction site (palindrome, if not palindrome, do sequences independently).
folder = 'Human' #folder with chromosome sequence files (note, do not put other files in this folder)
probe = 'TGGGCTTTCCTTTGAGGGTAACC' #without restriction site
stickyend = 4
#probe = 'TTAGGGTTAGGGTTAGGGTT'
#Depending on 5' or 3' overhang, one of the output sequences needs to be converted to the complementary strand  
import re #regular expressions lib
import os
import time # just a timer

start = time.time()
filenames = os.listdir(folder)
os.chdir( folder )
threeprime = []
fiveprime = []
positions = []
tot=0
totbp=0
restrict = restrict.upper() 
probe = probe.upper()
for filename in filenames:
    #Import genome as a single string
    my_file = open(filename)    
    file_contents = my_file.read()    
    #sequence = file_contents.strip('\n') #this doesn't work
    sequence = file_contents.replace("\n", "")
    #sequence = re.sub(r'\d+|\W+', '', sequence) #removes numbers etc. makes stuff slower for large sequence files
    #sequence = sequence.upper()      #changes everything to upper case, idem 
    count=0
    seqlist = re.finditer(restrict, sequence)
    for m in seqlist:
            count = count+1
            #print restrict,' @ ', m.start() 
            threeprime.append(sequence[(m.start()):(m.start()+len(probe)+len(restrict))]) #outputs all 5' and 3' sequences adiecent to restriction site
            fiveprime.append(sequence[(m.start()-len(probe)):m.start()+len(restrict)]) #empty sequences are near the ends of the DNA
    bp = len(sequence)
    #print 'sequence length =',bp ,'bp'
    print '# restriction sites',restrict,' in ',filename,' =', count   
    tot=tot+count
    totbp=totbp+bp
    my_file.close()
print 'total', restrict,' sites =',tot,'in', totbp, 'bp'    
        
glue = ','  #make single string with all the sequences
fiveprime = glue.join(fiveprime) #idem
threeprime = glue.join(threeprime)

def comp_string(base):  #make antisense
    if(base == "T"):
        return "A" 
    if(base == "A"):
        return "T" 
    if(base == "G"):
        return "C" 
    if(base == "C"):
        return "G"   
   
subprobe=''
subqrode=''
eborp=probe[::-1] #reverse probe
if stickyend == 5: 
    eborp=probe
for x in eborp:
    if stickyend == 5: 
        subprobe = subprobe + x
    else:
        subprobe = x+subprobe  
    y=comp_string(x)    
    if stickyend == 5: 
        subqrode = y + subqrode
    else:
        subqrode = subqrode + y
    findthis = subprobe+restrict    
    findthat = restrict+subqrode      
    if stickyend == 5: 
        findthat = restrict+subprobe 
        findthis = subqrode+restrict    
    count=0    
    #five=  fiveprime.count(findthis) #slightly slower and not compatible with Python 3 
    five=len(re.findall(findthis, fiveprime))
    #three=  fiveprime.count(findthis)     
    three=len(re.findall(findthat, threeprime))
    count = five + three
    print '# sequence compatible with',findthis,'and', findthat,'-->', five,'+',three,"=", count
    #print '# sequence compatible with',findthis,'-->', five
    #print '# sequence compatible with',findthat,'-->',three    
    if count == 0:
        break 

end = time.time()
print 'time =', end - start