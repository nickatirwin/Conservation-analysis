# conservation_analysis.py

# usage: python conservation_analysis.py [query protein] [path to database]
#
# note: the query protein should be in fasta format. The database can be eukaryota or viridiplantae
#
# example: python conservation_analysis.py proteinA.fasta /database/viridiplantae
#
# requirements: diamond, mafft, iqtree, ete3

# load modules
import sys
import subprocess
from ete3 import Tree
from collections import Counter
from Bio.Align import substitution_matrices
import numpy as np
from scipy.ndimage import gaussian_filter1d
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns


# load arguments
query = sys.argv[1]
database = str(sys.argv[2]).lower()
if database in ['archaeplastida','eukaryota']:
    database = '/resources/conservation/'+database

# define functions

# alignment trimming function - trim an alignment based on a query sequence
def query_trim(query,aln):
    # load alignment
    alignment = open(aln,'r').read().split('>')[1:]
    # get query
    for seq in alignment:
        if seq.split('\n')[0] == query:
            query_seq = seq.split('\n',1)[1].replace('\n','')
    # get non-gap positions in the query
    positions = {}
    n = 0
    for p in list(query_seq):
        if p != '-':
            positions[n] = True
        n += 1
    # trim the alignment based on query positions
    out = open(aln+'.trim','w')
    for seq in alignment:
        n = 0
        sequence = ''
        for p in list(seq.split('\n',1)[1].replace('\n','')):
            try:
                positions[n]
                sequence = sequence + p
            except:
                pass
            n += 1
        if sequence.count('-')/len(sequence) <= 0.5:
            out.write('>'+seq.split('\n')[0]+'\n'+sequence+'\n')
    out.close()

# identify homologs using Diamond BLAST

# run the blast
# parameters: 8 threads, evalue < 1e-5, query-coverage > 50%, max 2000 hits, sensitive search mode
print('Running BLAST search...')
subprocess.call('diamond blastp --query '+query+' --db '+database+'.dmnd --max-target-seqs 1000 --threads 8 --sensitive --outfmt 6 --evalue 1e-5 --query-cover 50 --out query.blastp',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# extract the hits
# load the sequences from the database
db_seqs = {}
db = open(database+'.fasta','r').read().split('>')[1:]
for seq in db:
    db_seqs[seq.split('\n')[0]] = '>'+seq
q = open(query,'r').read().split('>')[1:]
db_seqs[q[0].split('\n')[0]] = '>'+q[0]

# get the blast hits and output the sequence for each hit
blast = open('query.blastp','r').readlines()
out = open('query.blastp.fasta','w')
for hit in blast:
    out.write(db_seqs[hit.split('\t')[1]])
out.write(db_seqs[q[0].split('\n')[0]])
out.close()

# align the hits and make a phylogeny

# align the sequences using mafft
print('Aligning hits...')
subprocess.call('mafft --thread 4 query.blastp.fasta > query.blastp.fasta.aln',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# trim the alignment based on the query
query_name = q[0].split('\n')[0]
query_trim(query_name,'query.blastp.fasta.aln')

# make a phylogeny
print('Generating initial phylogeny...')
subprocess.call('iqtree -m LG+G4 -fast -alrt 1000 -nt 4 -s query.blastp.fasta.aln.trim -redo',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# analyze tree and extract closely related sequences
print('Extracting related sequences from initial phylogeny...')
t = Tree('query.blastp.fasta.aln.trim.treefile')
# midpoint root the tree
midpoint = t.get_midpoint_outgroup()
t.set_outgroup(midpoint)
# find query and extract related sequences - continue taking sequences to get close to the limit
# assign query sequence and identify the query node
matches = t.search_nodes(name=query_name)
query_node = list(matches)[0]
# go up nodes collecting sequences
goal_sequences = 200
clade = query_node
clade_d = {}
root = False
while root == False:
    m = len(clade.get_leaves())
    clade_d[clade] = m
    clade = clade.up
    if clade.is_root() == True:
        m = len(clade.get_leaves())
        clade_d[clade] = m
        root = True
# get clade with closest number of sequences to the max
clade_d = {c:abs(goal_sequences-clade_d[c]) for c in clade_d}
clade = min(clade_d, key=clade_d.get)
# get the sequences
sequences = [l.name for l in clade.get_leaves()]

# extract sequences, re-align
out = open('query.select.fasta','w')
for hit in sequences:
    out.write(db_seqs[hit])
out.close()

# align sequences
print('Aligning selected sequences...')
subprocess.call('mafft --thread 8 query.select.fasta > query.select.fasta.aln',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
query_trim(query_name,'query.select.fasta.aln')

# make phylogeny and infer site rates
print('Generating secondary phylogeny and inferring rates...')
subprocess.call('iqtree -m LG+G4 -fast -alrt 1000 -nt 4 -s query.select.fasta.aln.trim -wsr',shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

# extract amino acid conservation and plot
#
# return table with: protein\tamino_acid\tresidue\trate\trate_avg\tconservation\tconservation_avg\t%sites\n

# load arguments
query = query_name
rate = 'query.select.fasta.aln.trim.rate'
msa = 'query.select.fasta.aln.trim'

# calculate conservation values
rate = [l for l in open(rate,'r').readlines() if l.startswith('#') == False][1:]
residues = [a.split('\t')[0] for a in rate]
# set window size
window = int(len(residues)/100)
if window < 2:
    window = 3
rates = [float(a.split('\t')[1]) for a in rate]
rolling_rates = list(pd.Series(rates).rolling(window=window).mean())
rolling_rates = [i for i in rolling_rates if str(i) != 'nan']+(['nan']*(window-1))
conservation = []
for i in rates:
    if i != 0:
        conservation.append(1/i)
    else:
        # if the rate is zero (constant site), assign the highest observed conservation
        conservation.append(1/min([v for v in conservation if v != 0]))
rolling_conservation = list(pd.Series(conservation).rolling(window=window).mean())
rolling_conservation = [i for i in rolling_conservation if str(i) != 'nan']+(['nan']*(window-1))
# gaussian smooth the rolling average
roll = []
for v in rolling_conservation:
    if v == 'nan':
        roll.append(np.nan)
    else:
        roll.append(v)
x = np.arange(len(roll))
mask = ~np.isnan(roll)
data_interp = np.interp(x, x[mask], np.array(roll)[mask])
rolling_conservation = gaussian_filter1d(data_interp, sigma=1.0)

# get gap content
msa = open(msa,'r').read().split('>')[1:]
seq_d = {}
for s in msa:
    seq_d[s.split('\n')[0]] = s.split('\n',1)[1].replace('\n','')
gaps = []
n = 0
while n < len(seq_d[msa[0].split('\n')[0]]):
    gap_count = 0
    for seq in seq_d:
        if seq_d[seq][n] == '-':
            gap_count += 1
    gaps.append(gap_count)
    n += 1
sites = [((1-(g/len(seq_d)))*100) for g in gaps]
amino_acids = [a.strip() for a in db_seqs[query].split('\n')[1] if a.strip() != '']

# also calculate shannon entropy and median distance

def shannon_entropy(sites):
    # determine counts and probabilities
    counts = dict(Counter(sites))
    prob = {c:counts[c]/len(sites) for c in counts}
    entropy = -1*np.sum([(prob[c]*np.log2(prob[c])) for c in counts if c != '-'])
    if entropy == -0.0:
        entropy = 0.0
    return entropy

def distance(query_residue,site):
    matrix = substitution_matrices.load("BLOSUM62")
    # skip gaps
    if site == '-':
        d = 'NA'
    else:
        d = matrix[query_residue,site]
    return d

alignment = open('query.select.fasta.aln.trim','r').read().split('>')[1:]
aln_length = len(alignment[0].split('\n',1)[1].replace('\n',''))
query_sequence = db_seqs[query_name].split('\n',1)[1].replace('\n','')

# calculate shannon entropy
p = 0
shannon = []
while p < aln_length:
    site = [s.split('\n',1)[1][p].replace('\n','')[0] for s in alignment]
    entropy = shannon_entropy(site)
    shannon.append(entropy)
    p += 1
rolling_shannon = list(pd.Series(shannon).rolling(window=window).mean())
rolling_shannon = [i for i in rolling_shannon if str(i) != 'nan']+(['nan']*(window-1))
# gaussian smooth the rolling average
roll = []
for v in rolling_shannon:
    if v == 'nan':
        roll.append(np.nan)
    else:
        roll.append(v)
x = np.arange(len(roll))
mask = ~np.isnan(roll)
data_interp = np.interp(x, x[mask], np.array(roll)[mask])
rolling_shannon = gaussian_filter1d(data_interp, sigma=1.0)

# calculate median distance
p = 0
med_dist = []
while p < aln_length:
    distances = []
    query_res = query_sequence[p]
    site = [s.split('\n',1)[1][p].replace('\n','')[0] for s in alignment]
    for s in site:
        d = distance(query_res,s)
        if d != 'NA':
            distances.append(d)
    med_dist.append(np.median(distances))
    p += 1
rolling_dist = list(pd.Series(med_dist).rolling(window=window).mean())
rolling_dist = [i for i in rolling_dist if str(i) != 'nan']+(['nan']*(window-1))
# gaussian smooth the rolling average
roll = []
for v in rolling_dist:
    if v == 'nan':
        roll.append(np.nan)
    else:
        roll.append(v)
x = np.arange(len(roll))
mask = ~np.isnan(roll)
data_interp = np.interp(x, x[mask], np.array(roll)[mask])
rolling_dist = gaussian_filter1d(data_interp, sigma=1.0)

# append everything into a table
protein = [query for i in residues]
df = pd.DataFrame({'protein':protein,'residue':residues,'aa':amino_acids,'rate':rates,'rate_window':rolling_rates,'conservation':conservation,
    'rolling_conservation':rolling_conservation,'entropy':shannon,'rolling_entropy':rolling_shannon,'distance':med_dist,'rolling_distance':rolling_dist,'sites':sites})
df.to_csv(protein[0]+'.rates.tsv', sep='\t', index=False)
df = pd.read_csv(protein[0]+'.rates.tsv',sep='\t')

# plot data
matplotlib.rcParams['pdf.fonttype'] = 42  # Embed fonts as TrueType (Type 42)
matplotlib.rcParams['ps.fonttype'] = 42

fig, ax1 = plt.subplots(figsize=(12, 6))
ax2 = ax1.twinx()
sns.lineplot(df,x='residue',y='rolling_conservation',color='#2166ac',linestyle='-',linewidth=1,ax=ax1)
sns.lineplot(df,x='residue',y='conservation',color='#92c5de',linestyle='-',linewidth=0.5,ax=ax1)
sns.lineplot(df,x='residue',y='sites',color='#bababa',linestyle='-',linewidth=1,ax=ax2)
ax2.set_ylim(0, 102)
ax1.set_xlabel('Residue')
ax1.set_ylabel('Conservation (1/rate)')
ax2.set_ylabel('Sites (%)')
plt.title(protein[0]+' (alignment sequences: ' + str(len(seq_d))+')', loc='left', fontsize=12)
plt.savefig(protein[0]+'.conservation_plot.pdf', format='pdf', bbox_inches='tight')
plt.show()

# plot data
matplotlib.rcParams['pdf.fonttype'] = 42  # Embed fonts as TrueType (Type 42)
matplotlib.rcParams['ps.fonttype'] = 42

# move files and make result directory
subprocess.call('mkdir '+query_name+'_result',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.call('mv *blastp* '+query_name+'_result/',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.call('mv *select* '+query_name+'_result/',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.call('mv *tsv '+query_name+'_result/',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
subprocess.call('mv *pdf '+query_name+'_result/',shell=True,stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


