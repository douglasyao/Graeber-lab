__author__ = 'NBalanis'

#takes in a .txt file of the form. where the values are the log2 ratio
#Sample	                     Chromosome	  Start	        End	   Num_Probes Segment_Mean
#TCGA-A1-A0SB-01A-11D-A141-01	1	     3218610	247813706	129072	     0.0034
#TCGA-A2-A0D2-01A-21D-A036-01	7	44642801	45813188	640	0.1997
#TCGA-A2-A0D2-01A-21D-A036-01	7	45813795	45862003	7	-0.4057
#species is either 'h' or 'm'
import sys,os,math
from collections import defaultdict
from collections import OrderedDict
from operator import itemgetter

#argument 1 is the input .seg file 
#argument 2 is either h or m
# argument 3 is the name you want to prefix the output file

input = sys.argv[1]
species = sys.argv[2]
prefix=sys.argv[3]
chrom_num = 0
d=OrderedDict()
e=OrderedDict()

def any_no_num(s):
    return not all(i.isdigit() for i in s)

if species == 'h':
    chrom_num = 22
elif species == 'm':
	chrom_num=19
else:
    print('Unknown species, %s is not an acceptable species, try h or m' % (argv[2]))

#load CNA score and breakpoints in separate ordered dictionaries linking (e) sample ~ CNA score   ... (f) sample ~ breakpoints
with open(input) as f:
    next (f)
    for line in f:
        line=line.rstrip()
        sample=line.split('\t')[0]
        chromi=line.split('\t')[1]
        chromi=chromi[3:]
        if any_no_num(chromi): 
            continue 
        if line.split('\t')[5] == 'NA':
            continue
        mean= abs(float(line.split('\t')[5]))
        diff= abs(float(line.split('\t')[3]) - float(line.split('\t')[2]))
        CNA_score= float(mean*diff)
        if sample not in d:
            d[sample]=1
            e[sample]=CNA_score
        else:
            d[sample] +=1
            e[sample] += CNA_score
f.close()

#sort the dictionaries
f = OrderedDict(sorted(d.items(),key=itemgetter(1),reverse=True))
g = OrderedDict(sorted(e.items(),key=itemgetter(1),reverse=True))


#write the dictionarys to a file
breakpoints_file= open(prefix + '_bkpts_sorted_by_bkpts_per_chrom.txt', 'w')
breakpoints_file.write('Sample' + '\t' + "bkpt_chrom"+'\t'+'bkpt_samp'+'\t'+'ICNA_score' + '\n')
for key in f:
    bkpt_samp_chrom = round((f[key]-chrom_num)/chrom_num,2)
    total_breakpoints = f[key]
    ICNA_score= round(e[key])
    breakpoints_file.write( str(key) + '\t' + str(bkpt_samp_chrom) + '\t' + str(total_breakpoints) + '\t' + str(ICNA_score) +'\n')
breakpoints_file.close()

#breakpoints_file2= open(prefix + '_bkpts_sorted_by_ICNA_per_sample.txt', 'w')
#breakpoints_file2.write('Sample' + '\t' + "bkpt_chrom"+'\t'+'bkpt_samp'+'\t'+'ICNA_score' + '\n')
#for key in g:
    # total_breakpoints = f[key]
  #  ICNA_score= round(g[key])
 #   breakpoints_file2.write( str(key) + '\t' + str(bkpt_samp_chrom) + '\t' + str(total_breakpoints) + '\t' + str(ICNA_score) +'\n')
#breakpoints_file2.close()





