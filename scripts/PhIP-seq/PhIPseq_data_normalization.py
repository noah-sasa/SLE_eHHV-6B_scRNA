#!/usr/bin/env python

'''
Author: Shohei Kojima @ RIKEN (shohei.kojima@riken.jp)
Usage: python %prog fq1.gz fq2.gz result.tsv
Version:
    Python 3.7.4
    pandas 1.3.1
    numpy 1.21.6
'''

import glob
import pandas as pd
import numpy as np


input_file_base = 'Input'
input_file_bases = ['Input_1', 'Input_2', 'Input_3', 'Input_4']


# data
files = glob.glob('./results/*.tsv')  # e.g., sample_1.oligo_count.tsv
files = sorted(files)

samples = []
for f in files:
    sample = f.split('/')[-1].split('.')[0]
    if 'NSC' in s:
        continue
    if 'Input' in s:
        continue
    samples.append(sample)

replicates = {}
for sample in samples:
    base = sample.split('_')[0]
    if not base in replicates:
        replicates[base] = []
    replicates[base].append(sample)

input_oligo_threshold_factor = 0.05
input_spikein_copy = 2000
outfilename = 'summary_results_by_ranges.hhv6_oligos.tsv'



# load hhv6 protein lengths and split
f = files[0]
df = pd.read_table(f, index_col = 0, comment = '#')
df = df[ df['is_HHV6'] == True ]
# protein index start and end
proteins = {}
for n, protein in enumerate(df['protein']):
    name = protein.split('_prot_', 1)[1].rsplit('_', 1)[0]
    if not name in proteins:
        proteins[name] = []
    proteins[name].append(n)
for name in proteins:
    proteins[name] = [min(proteins[name]), max(proteins[name]) + 1]
# protein length and virus name
protein_lengths = {}
protein_to_virus = {}
for oligoname,row in df.iterrows():
    protein = row['protein'].split('_prot_', 1)[1].rsplit('_', 1)[0]
    if not protein in protein_lengths:
        protein_lengths[protein] = 0
    end = row['end']
    protein_lengths[protein] = max(end, protein_lengths[protein])
    if 'h6a' in oligoname:
        virus = 'h6a'
    elif 'h6b' in oligoname:
        virus = 'h6b'
    protein_to_virus[protein] = virus
# merge
for protein in proteins:
    proteins[protein].append(protein_lengths[protein])
    proteins[protein].append(protein_to_virus[protein])


# load input
inputs = []
for base in input_file_bases:
    f = './results/%s.oligo_count.tsv' % base
    df = pd.read_table(f, index_col = 0, comment = '#')
    df = df[ df['is_HHV6'] == True ]
    df = df[['count']]
    df.columns = [base]
    inputs.append(df)
df = pd.concat(inputs, axis = 1)
df['count'] = df.sum(axis = 1)
input_oligo_count = {}
for oligoname,row in df.iterrows():
    input_oligo_count[oligoname] = row['count']
input_total = sum(list(input_oligo_count.values()))
input_mean = input_total / len(input_oligo_count)
input_oligo_threshold = input_mean * input_oligo_threshold_factor
print('input_mean:', input_mean)
print('input_oligo_threshold:', input_oligo_threshold)
print(df[ df['count'] < input_oligo_threshold ].shape)
print(df[ df['count'] >= input_oligo_threshold ].shape)

input_norm_factors = {}
absent_oligo = set()
for i in input_oligo_count:
    if input_oligo_count[i] < input_oligo_threshold:
        absent_oligo.add(i)
        input_norm_factors[i] = 1
        continue
    input_norm_factors[i] = input_mean / input_oligo_count[i]


# ranges
ranges = []
idf = df.index.tolist()
for protein in proteins:
    s, e, prot_len, virus = proteins[protein]
    for i in range(e - s + 4):
        _s = max(s, s + i - 3)
        _e = min(e, s + i + 1)
        if _s == _e:
            _s -= 1  # last 14-aa
        is_input0 = True
        for j in idf[_s:_e]:
            if not j in absent_oligo:
                is_input0 = False
                break
        r = (virus, protein, i * 14, min((i+1)*14, prot_len))
        ranges.append((r, is_input0))


# load data by sample (replicates)
results = {}
for sample in replicates:
    results[sample] = {}
    dfs = []
    for rep in replicates[sample]:
        f = './results/%s.oligo_count.tsv' % rep
        df = pd.read_table(f, index_col = 0, comment = '#')
        n_spikein = df.iloc[0,5]
        df = df[ df['is_HHV6'] == True ]
        df['count'] = df['count'] * input_spikein_copy / n_spikein  # calc. copy number
        dfs.append(df[['count']])
    df = pd.concat(dfs, axis = 1)
    norm_factors = [ input_norm_factors[i] for i in df.index ]
    df *= np.array(norm_factors).reshape(-1, 1)  # normalize by input
    ndf = df.to_numpy()
    idf = df.index.tolist()
    range_i = 0
    for protein in proteins:
        s, e, prot_len, virus = proteins[protein]
        for i in range(e - s + 4):
            _s = max(s, s + i - 3)
            _e = min(e, s + i + 1)
            if _s == _e:
                _s -= 1  # last 14-aa
            tmp = ndf[_s:_e,:]
            tmp = tmp[[ not j in absent_oligo for j in idf[_s:_e] ]]
            if tmp.shape[0] == 0:
                score = 0
            else:
                score = tmp.mean()
            r, _ = ranges[range_i]
            range_i += 1
            results[sample][r] = score


out = []
samples = []
for sample in results:
    samples.append(sample)
tmp = 'virus\tprotein_id\tstart\tend\tis_input0\t%s\n' % '\t'.join(samples)
out.append(tmp)
for r, is_input0 in ranges:
    scores = []
    for sample in results:
        scores.append('%f' % results[sample][r])
    scores = '\t'.join(scores)
    tmp = '%s\t%s\t%s\t%d\t%d\t%s\n' % (*r, is_input0, scores)
    out.append(tmp)

with open(outfilename, 'w') as outfile:
    outfile.write(''.join(out))
