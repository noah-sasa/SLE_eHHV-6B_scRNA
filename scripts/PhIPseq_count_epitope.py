#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN (shohei.kojima@riken.jp)
Usage: python %prog fq1.gz fq2.gz result.tsv
Version:
    Python 3.7.4
    pandas 1.3.1
'''

import sys,gzip,collections,string
import pandas as pd


fq1 = sys.argv[1]
fq2 = sys.argv[2]
outfilename = sys.argv[3]


######## Parameters ########
trim_read_len = 10
hash_length = 20
min_insert_len = 20
allow_phage_mismatch = 1
allow_read1_mismatch = 1
allow_read2_mismatch = 1
############################


######## Pre-defined sequences ########
HindIII_GGGS = 'GGAGCTGTCGTATTCCAGTCAGGTGTGATGCTCGGGGATCCGAATTCTCCTGCAGGGATATCCCGGGAGCTCGTCGACAAGCTTGGTGGCGGTTCA'
XhoI = 'AACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTAACTAGTTACTCGAG'
spikein = 'ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACC'
HindIII_spikein = 'AACCCCTCAAGACCCGTTTAGAGGCCCCAAGGGGTTAACTAGTTAAAGCTT'
aa_spikein = 'MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTT'
#######################################


len_HindIII_GGGS = len(HindIII_GGGS)
len_XhoI = len(XhoI)


def complement(seq):
    return seq.translate(str.maketrans('ATGCatgc', 'TACGtacg'))[::-1]

def is_same_seq(seq1, seq2, mismatch_threshold):
    mismatch = 0
    for c1, c2 in zip(seq1, seq2):
        if c1 == c2:
            continue
        mismatch += 1
        if mismatch > mismatch_threshold:
            return False
    return True


# load HHV peptide info
id_to_aa_pos = {}
files = ['/path/to/HHV6A_peptides.aa.fa', '/path/to/HHV6B_peptides.aa.fa']
for f in files:
    with open(f) as infile:
        for line in infile:
            ls = line.split()
            if not line[0] == '>':
                continue
            id, prot = line.replace('>', '').split()
            prot, pos = prot.split(':')
            start, end = pos.split('-')
            id_to_aa_pos[id] = (prot, start, end)
# calc protein length
hhv_prot_length = {}
for id in id_to_aa_pos:
    prot, start, end = id_to_aa_pos[id]
    if not prot in hhv_prot_length:
        hhv_prot_length[prot] = []
    hhv_prot_length[prot].append(end)
for prot in hhv_prot_length:
    hhv_prot_length[prot] = max(hhv_prot_length[prot])
# load positive control peptide info
f = '/path/to/positive_control_bacteria.tsv'
ctrl_pep_info = {}
df = pd.read_table(f, index_col = 0)
for i, r, f in zip(df.index, df['uniref'], df['uniref_func']):
    f = f.split(' n=')[0]
    ctrl_pep_info[i] = '%s:%s' % (r, f)


seq_to_id = {}
oligo_ids = []
hhv_ids = set()
# spike-in
seq = spikein
id = 'spikein'
seq_head = seq[:hash_length]
if not seq_head in seq_to_id:
    seq_to_id[seq_head] = []
seq_to_id[seq_head].append((id, seq))
oligo_ids.append((id, seq))
# hhv6 oligos
f = '/path/to/HHV6_oligo_pool.xlsx'
hhv = pd.read_excel(f, index_col = None)
for id, seq in zip(hhv['name'], hhv['sequence']):
    seq = seq[35:-23]
    seq_head = seq[:hash_length]
    if not seq_head in seq_to_id:
        seq_to_id[seq_head] = []
    seq_to_id[seq_head].append((id, seq))
    oligo_ids.append((id, seq))
    hhv_ids.add(id)
# Positive control oligos
f = '/path/to/PC_oligo_pool.xlsx'
hhv = pd.read_excel(f, index_col = None)
for id, seq in zip(hhv['name'], hhv['sequence']):
    seq = seq[35:-23]
    seq_head = seq[:hash_length]
    if not seq_head in seq_to_id:
        seq_to_id[seq_head] = []
    seq_to_id[seq_head].append((id, seq))
    oligo_ids.append((id, seq))

def judge_hhv_id(id):
    if id in hhv_ids:
        return True
    return False



infile1 = gzip.open(fq1, 'rt')
infile2 = gzip.open(fq2, 'rt')

n_processed = 0
n_phage = 0
n_phage_short = 0
n_phage_assigned = 0
n_phage_unassigned = 0
n_not_phage = 0
oligo_found = collections.Counter()
for _, _ in zip(infile1, infile2):
    n_processed += 1
    if n_processed % 500_000 == 0:
        print(n_processed, 'spots processed, n_phage=%d, n_phage_assigned=%d, n_spikein=%d' % (n_phage, n_phage_assigned, oligo_found['spikein']))
    read1 = next(infile1).strip()
    read2 = next(infile2).strip()
    for _ in range(2):
        next(infile1)
        next(infile2)
    is_read1_same = is_same_seq(read1[trim_read_len:len_HindIII_GGGS], HindIII_GGGS[trim_read_len:], allow_phage_mismatch)
    is_read2_same = is_same_seq(read2[trim_read_len:len_XhoI], XhoI[trim_read_len:], allow_phage_mismatch)
    is_read2_spikein = is_same_seq(read2[trim_read_len:len_XhoI], HindIII_spikein[trim_read_len:], allow_phage_mismatch)
    if (is_read1_same and is_read2_same) or (is_read1_same and is_read2_spikein):
        n_phage += 1
        insert1 = read1[len_HindIII_GGGS:]
        insert2 = read2[len_XhoI:]
        insert2 = complement(insert2)
        len_insert1 = len(insert1)
        len_insert2 = len(insert2)
        if len_insert1 < min_insert_len:
            n_phage_short += 1
            continue
        if len_insert2 < min_insert_len:
            n_phage_short += 1
            continue
        if len_insert1 < hash_length:
            n_phage_short += 1
            continue
        seq_head = insert1[:hash_length]
        if not seq_head in seq_to_id:
            n_phage_unassigned += 1
            continue
        match = []
        for id, seq in seq_to_id[seq_head]:
            is_insert1_same = is_same_seq(seq[:len_insert1], insert1, allow_read1_mismatch)
            is_insert2_same = is_same_seq(seq[-len_insert2:], insert2, allow_read2_mismatch)
            if is_insert1_same and is_insert2_same:
                match.append(id)
        if len(match) >= 1:
            n_phage_assigned += 1
            for id in match:
                oligo_found[id] += 1
        else:
            n_phage_unassigned += 1
    else:
        n_not_phage += 1


print('n_processed', n_processed)
print('n_phage', n_phage)
print('n_not_phage', n_not_phage)
print('n_phage_assigned', n_phage_assigned)
print('n_phage_unassigned', n_phage_unassigned)
print('n_phage_short', n_phage_short)
print('n_spikein', oligo_found['spikein'])
print()
stat = '#n_processed=%d;n_phage=%d;n_not_phage=%d;n_phage_assigned=%d;n_phage_unassigned=%d;n_phage_short=%d;n_spikein=%d\n' % (
    n_processed, n_phage, n_not_phage, n_phage_assigned, n_phage_unassigned, n_phage_short, oligo_found['spikein'],
)


out = [stat, 'ID\tis_HHV6\tprotein\tstart\tend\tprotein_length\tcount\tinsert\n']
for id, seq in oligo_ids:
    is_hhv = judge_hhv_id(id)
    if is_hhv:
        short_id = id.split('_', 1)[1]
        prot, start, end = id_to_aa_pos[short_id]
        length = hhv_prot_length[prot]
    else:
        if id == 'spikein':
            prot = aa_spikein
        else:
            short_id = 'oligo_' + id.rsplit('_', 1)[1]
            prot = ctrl_pep_info[short_id]
        start, end = 0, 0
        length = 0
    tmp = '%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\n' % (id, is_hhv, prot, start, end, length, oligo_found[id], seq)
    out.append(tmp)

with open(outfilename, 'w') as outfile:
    outfile.write(''.join(out))
