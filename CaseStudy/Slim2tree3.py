import pyslim
import tskit
import sys
from Bio import Phylo
from io import StringIO
import re
import random
import os
import subprocess
# for i in 30 40 50 60 65 70 75 80 85 90; do python Slim2tree3.py $i 1; done
AGE=sys.argv[1]
REP=sys.argv[2]
InputFile="CHIPage_"+AGE+"_rep"+REP+".trees"
if os.path.exists("CHIPage_"+AGE+"_rep"+REP+"Cellinfo.txt")==True:
   os.remove("CHIPage_"+AGE+"_rep"+REP+"Cellinfo.txt")
os.system('slim -d AGE='+AGE+' -d REP='+REP+' CHIP2EnvN.slim')




SampleN=400
#ts = pyslim.load(InputFile)#"output.trees")
ts = tskit.load(InputFile)      # plain tskit TreeSequence
n_tips=0
for tree in ts.trees():

    n_tips =+ tree.num_samples()  # number of leaf (sample) nodes in this tree
print(f"Tree: {n_tips} tips")
#print (tip_birth_times)
#open('a','r').readlines()
# Get individuals alive at the end (time = 0 in SLiM)
final_gen = 0
alive_inds = pyslim.individuals_alive_at(ts, final_gen)
print(f"Alive individuals at time {final_gen}: {len(alive_inds)}")
print (alive_inds[:5],SampleN)
# These are the individuals SLiM would have used for VCF output
# Assuming you used p1.outputVCFSample(10)
if len(alive_inds)<SampleN:
    sample_inds = alive_inds
else: 
   # sample_inds = alive_inds[:SampleN]
    
    sample_inds = random.sample(list(alive_inds), SampleN)
# Get their first node (haploid) or both (if diploid)
sample_nodes = [ts.individual(i).nodes[0] for i in sample_inds]

# Simplify the tree sequence to only the sampled nodes
ts_simplified = ts.simplify(samples=sample_nodes)
print (len(sample_nodes),len(ts_simplified.samples()),sample_nodes[:3],ts_simplified.samples()[:3])

print (ts_simplified.num_trees,len(sample_nodes))
tree = ts_simplified.first()
#tree.draw_node_labels = True
newicks = [tree.as_newick(root=root) for root in tree.roots]
newick_file = InputFile + ".nwk"
OutF=open(newick_file,'w')
OutF.write('\n'.join(newicks))
OutF.close()

