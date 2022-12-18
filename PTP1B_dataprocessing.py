# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 09:14:13 2022

@author: Andrew Hren, University of Colorado, Boulder
         andrew.hren@colorado.edu
    
"""

from Bio import SeqIO, Align
from Bio.Seq import Seq
import pandas as pd
import csv
import re
import ast
import math
import numpy as np


## Vsearch and usearch are required for multiple steps in this script

## General workflow:
## Hifi reads are separated into separate fastas based on exact barcode matches
## Sequences are translated to proteins then protein seqs are value counted.
## The output (pcounts) is run through usearch to generate the qrowdots string
## Mutations are ID'd, then comparison data is exported to csvs

## PacBio barcodes
## Includes full 16bp barcodes, then 1 and 2bp truncations from the 5' end
## Basically, requires 14bp of the barcode on one end to identify.
## Tons of reads were being tossed despite me going in and being able to tell
## exactly where it came from. Even with 14bp, there's still 5-6 mutations
## required to get misidentified. It's much more likely I non-identify than 
## mis-identify.




def idx_return(barcode,barcode15,barcode14,direction):
    
    barcodes_fwd = ['CACATATCAGAGTGCG','ACACACAGACTGTGAG','ACACATCTCGTGAGAG',
                    'CATATATATCAGCTGT','TCTGTATCTCTATGTG','ACAGTCGAGCGCTGCG',
                    'ACATATCAGAGTGCG','CACACAGACTGTGAG','CACATCTCGTGAGAG',
                    'ATATATATCAGCTGT','CTGTATCTCTATGTG','CAGTCGAGCGCTGCG',
                    'CATATCAGAGTGCG','ACACAGACTGTGAG','ACATCTCGTGAGAG',
                    'TATATATCAGCTGT','TGTATCTCTATGTG','AGTCGAGCGCTGCG']
                    
    barcodes_rev = ['CGCACTCTGATATGTG','CTCACAGTCTGTGTGT','CTCTCACGAGATGTGT',
                    'ACAGCTGATATATATG','CACATAGAGATACAGA','CGCAGCGCTCGACTGT',
                    'GCACTCTGATATGTG','TCACAGTCTGTGTGT','TCTCACGAGATGTGT',
                    'CAGCTGATATATATG','ACATAGAGATACAGA','GCAGCGCTCGACTGT',
                    'CACTCTGATATGTG','CACAGTCTGTGTGT','CTCACGAGATGTGT',
                    'AGCTGATATATATG','CATAGAGATACAGA','CAGCGCTCGACTGT']
    
    if direction == 'fwd':
        if barcode in barcodes_fwd:
            idx = barcodes_fwd.index(barcode)
        if barcode15 in barcodes_fwd:
            idx = barcodes_fwd.index(barcode15)  
        if barcode14 in barcodes_fwd:
            idx = barcodes_fwd.index(barcode14)   
            
    if direction == 'rev':
        if barcode in barcodes_rev:
            idx = barcodes_rev.index(barcode)
        if barcode15 in barcodes_rev:
            idx = barcodes_rev.index(barcode15)
        if barcode14 in barcodes_rev:
            idx = barcodes_rev.index(barcode14)        
        
    if idx >= 12:
        idx = idx - 12
    if idx >= 6:
        idx = idx - 6
    if idx <= 5:
        idx = idx   
    
    return idx
        



def formatter(row):
    seq_id = str(row.name)
    if len(seq_id) != 6:
        while len(seq_id) != 6:
            seq_id = '0' + seq_id   
    final_id = str(seq_id) + '-' + str(row.counts)
    info = [final_id,row.seq]
    return info

def count_extract(row):
    query = str(row.query)
    count = query[7:]
    return count

def translate(row):
    # translate the CDS of ptp1b by finding the six bases upstream of ATG
    prot_seq = str(Seq(str(row.qseq)).translate())
    return prot_seq

def SNP_indexer(row):
    #write out mutations by going through qrowdots output, comparing to original
    base_dna = 'ATGGAGATGGAAAAGGAGTTCGAGCAGATCGACAAGTCCGGGAGCTGGGCGGCCATTTACCAGGATATCCGACATGAAGCCAGTGACTTCCCATGTAGAGTGGCCAAGCTTCCTAAGAACAAAAACCGAAATAGGTACAGAGACGTCAGTCCCTTTGACCATAGTCGGATTAAACTACATCAAGAAGATAATGACTATATCAACGCTAGTTTGATAAAAATGGAAGAAGCCCAAAGGAGTTACATTCTTACCCAGGGCCCTTTGCCTAACACATGCGGTCACTTTTGGGAGATGGTGTGGGAGCAGAAAAGCAGGGGTGTCGTCATGCTCAACAGAGTGATGGAGAAAGGTTCGTTAAAATGCGCACAATACTGGCCACAAAAAGAAGAAAAAGAGATGATCTTTGAAGACACAAATTTGAAATTAACATTGATCTCTGAAGATATCAAGTCATATTATACAGTGCGACAGCTAGAATTGGAAAACCTTACAACCCAAGAAACTCGCGAGATCTTACATTTCCACTATACCACATGGCCTGACTTTGGAGTCCCTGAATCACCAGCCTCATTCTTGAACTTTCTTTTCAAAGTCCGAGAGTCAGGGTCACTCAGCCCGGAGCACGGGCCCGTTGTGGTGCACTGCAGTGCAGGCATCGGCAGGTCTGGAACCTTCTGTCTGGCTGATACCTGCCTCTTGCTGATGGACAAGAGGAAAGACCCTTCTTCCGTTGATATCAAGAAAGTGCTGTTAGAAATGAGGAAGTTTCGGATGGGGCTGATCCAGACAGCCGACCAGCTGCGCTTCTCCTACCTGGCTGTGATCGAAGGTGCCAAATTCATCATGGGGGACTCTTCCGTGCAGGATCAGTGGAAGGAGCTTTCCCACGAGGACCTGGAGCCCCCACCCGAGCATATCCCCCCACCTCCCCGGCCACCCAAACGAATCCTGGAGCCACACAATCTCGAGCATCATCATCATCATCATTGA'
    base_dna = [*base_dna]
    variant = [*row.qrowdots]
    mut_list = []
    char_num = 0
    for char in variant:
        
        if char != '.':
            try:
                mut = base_dna[char_num] + str(char_num + 1) + variant[char_num]
                mut_list.append(mut)
            except:
                mut_list.append('NA')
        char_num += 1
    return mut_list

def extract1(row):
    # extracting data from the fasta titles I created, sample name
    info = row.query
    info_list = info.split(sep='-')
    return info_list[0]
def extract2(row):
    # extracting data from the fasta titles I created, pcounts
    info = row.query
    info_list = info.split(sep='-')
    return info_list[1]
def extract3(row):
    # extracting data from the fasta titles I created, reads
    info = row.query
    info_list = info.split(sep='-')
    return int(info_list[2])

def mutation_indexer(row):
    #write out mutations by going through qrowdots output, comparing to original
    base_protein = 'MEMEKEFEQIDKSGSWAAIYQDIRHEASDFPCRVAKLPKNKNRNRYRDVSPFDHSRIKLHQEDNDYINASLIKMEEAQRSYILTQGPLPNTCGHFWEMVWEQKSRGVVMLNRVMEKGSLKCAQYWPQKEEKEMIFEDTNLKLTLISEDIKSYYTVRQLELENLTTQETREILHFHYTTWPDFGVPESPASFLNFLFKVRESGSLSPEHGPVVVHCSAGIGRSGTFCLADTCLLLMDKRKDPSSVDIKKVLLEMRKFRMGLIQTADQLRFSYLAVIEGAKFIMGDSSVQDQWKELSHEDLEPPPEHIPPPPRPPKRILEPHNLEHHHHHH*AAADPMVRVLEA'
    base_protein = [*base_protein]
    variant = [*row.qrowdots]
    mut_list = []
    char_num = 0
    for char in variant:
        
        if char != '.':
            try:
                mut = base_protein[char_num] + str(char_num + 1) + variant[char_num]
                mut_list.append(mut)
            except:
                mut_list.append('NA')
        char_num += 1
    return str(mut_list)

def log2fc(row):
    lfc = math.log(row.freq_y/row.freq_x,2)
    return lfc

def formatter(index,count,sequence):
    index = str(index)
    
    if len(index) != 6:
        while len(index) != 6:
            index = '0' + index   
    final_id = str(index) + '-' + str(count)
    info = [final_id,sequence]
    return info

def asstring(row):
    litstring = str(row.mutations)
    out = ast.literal_eval(litstring)
    return out

def extract_position(row):
    pos = re.findall(r'\d+', str(row.mutations))   
    pos = pos[0]     
    return pos

def aa_category(row):
    positive = ['H','K','R']
    negative = ['D','E']
    aromatic = ['F','Y','W']
    sulfuric = ['C','M']
    amide = ['N','Q']
    polar = ['S','T']
    nonpolar = ['G','A','V','L','I','P']
    types = [positive,negative,aromatic,sulfuric,amide,polar,nonpolar]
    names = ['positive','negative','aromatic','sulfuric','amide','polar','nonpolar']
    letters = re.findall(r'[a-zA-Z]+',row.mutations[0])
    if len(letters) == 1:
        change = 'deletion'
    else:
        x = 0
        for sublist in types:
            if letters[0] in sublist:
                changeA = names[x]
            if letters[1] in sublist:
                changeB = names[x]
            x += 1
        change = changeA + '-' + changeB
    return change

def mutation_count(row):
    litstring = str(row.mutations)
    out = len(ast.literal_eval(litstring))
    return out


## Step 1: initialize sample names

x = 0
## This setion defines the sample names
samples = ['R1_0','R1_375','R1_75','R2_0','R2_375','R2_75']
output_files = []
for sample_name in samples:
    output_name = './outputs/' + sample_name + '.fasta'
    output_files.append(output_name)

## Step 2: convert hifi reads FASTQ into FASTA using vsearch
## (too big for usearch)
#  vsearch --fastx_filter ../hifi_reads.fastq --fastaout fasta_reads.fasta --fastq_qmax 93   




## Step 3: split into individually barcoded samples
## Reads with 14+ bases of barcoded homology are binned
## in their corresponding files
"""
    barcodes_fwd = ['CACATATCAGAGTGCG','ACACACAGACTGTGAG','ACACATCTCGTGAGAG',
                    'CATATATATCAGCTGT','TCTGTATCTCTATGTG','ACAGTCGAGCGCTGCG',
                    'ACATATCAGAGTGCG','CACACAGACTGTGAG','CACATCTCGTGAGAG',
                    'ATATATATCAGCTGT','CTGTATCTCTATGTG','CAGTCGAGCGCTGCG',
                    'CATATCAGAGTGCG','ACACAGACTGTGAG','ACATCTCGTGAGAG',
                    'TATATATCAGCTGT','TGTATCTCTATGTG','AGTCGAGCGCTGCG']
                    
    barcodes_rev = ['CGCACTCTGATATGTG','CTCACAGTCTGTGTGT','CTCTCACGAGATGTGT',
                    'ACAGCTGATATATATG','CACATAGAGATACAGA','CGCAGCGCTCGACTGT',
                    'GCACTCTGATATGTG','TCACAGTCTGTGTGT','TCTCACGAGATGTGT',
                    'CAGCTGATATATATG','ACATAGAGATACAGA','GCAGCGCTCGACTGT',
                    'CACTCTGATATGTG','CACAGTCTGTGTGT','CTCACGAGATGTGT',
                    'AGCTGATATATATG','CATAGAGATACAGA','CAGCGCTCGACTGT']

with open('fasta_reads.fasta') as f:    
    merge = SeqIO.parse('fasta_reads.fasta', 'fasta')
    for contig in merge:
        barcode = contig.seq[0:16]
        barcode15 = contig.seq[0:15]
        barcode14 = contig.seq[0:14]
        rev_barcode = contig.seq[-16:]
        rev_b15 = contig.seq[-15:]
        rev_b14 = contig.seq[-14:]
        
        if barcode in barcodes_fwd or barcode15 in barcodes_fwd or barcode14 in barcodes_fwd:
            idx = idx_return(barcode,barcode15,barcode14,'fwd')
            sample = samples[idx]
            with open(output_files[idx],'a') as output:
                output.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')
                       
        elif barcode in barcodes_rev or barcode15 in barcodes_rev or barcode14 in barcodes_rev:
            idx = idx_return(barcode,barcode15,barcode14,'rev')   
            sample = samples[idx]
            with open(output_files[idx],'a') as output:
                output.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')
               
        elif rev_barcode in barcodes_fwd or rev_b15 in barcodes_fwd or rev_b14 in barcodes_fwd:
            idx = idx_return(rev_barcode,rev_b15,rev_b14,'fwd')    
            sample = samples[idx]
            with open(output_files[idx],'a') as output:
                output.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')
        
        elif rev_barcode in barcodes_rev or rev_b15 in barcodes_rev or rev_b14 in barcodes_rev:
            idx = idx_return(rev_barcode,rev_b15,rev_b14,'rev')   
            sample = samples[idx]
            with open(output_files[idx],'a') as output:
                output.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')

                
        else:
            with open('bad_output.fasta','a') as output:
                output.write('>'+str(contig.id)+'\n'+str(contig.seq)+'\n')
            x += 1
            if x%1000 == 0:
                print('bad barcode count: ',x)
                print(barcode)



## Step 4:
## Count reads in each fasta file
sample_reads = []
for sample in samples:
    print(sample)
    path1 = './outputs/' + sample + '.fasta'
    records = SeqIO.parse(path1,'fasta')
    sizes = [rec for rec in records]
    print(len(sizes))
    #data = pd.read_csv(path1,names=['qseq','reads'],skiprows=[0])
    #sample_reads.append(data['reads'].sum())

sample_reads = [667568,717434,780752,649875,738381,698341]    

## Step 5:
## Value count unique DNA sequences:
## This code searches for the region directly upstream of 
## the coding sequence. If it's not present, it gets 
## the reverse complement. Then it extracts the coding
## sequence, value counts the unique
## dna sequences, then outputs to csv.



barcode_index = 0
for sample in samples:
    dna_seqs = []
    path = './outputs/' + sample + '.fasta'
    with open(path,'r') as file:
        data = SeqIO.parse(path, 'fasta')
        for contig in data:
            barcode = contig.seq[0:16]
            barcode15 = contig.seq[0:15]
            barcode14 = contig.seq[0:14]
            rev_barcode = contig.seq[-16:]
            rev_b15 = contig.seq[-15:]
            rev_b14 = contig.seq[-14:]
            sequence = contig.seq

            if 'TAAGTGCAGAAAGAGGAGAAATACTAGACCGGT' in sequence:
                backwards = False
            else:
                backwards = True
             
            if backwards == True:  
                sequence = sequence.reverse_complement()
            sequence = str(sequence)
            ptp_start = sequence.find('AGACCGGT') + 8
            ptp_seq = sequence[ptp_start:ptp_start+990]
             
            dna_seqs.append(ptp_seq)
        print(len(dna_seqs))
        dna_seq_df = pd.DataFrame(dna_seqs)
        dna_counts = pd.DataFrame(dna_seq_df.value_counts())
        dna_counts.columns=['counts']
        output_name = './outputs/' + sample + '_dna_counts.csv'
        dna_counts.to_csv(output_name)  
    barcode_index += 1


## Step 6:
## This is the new fasta file creater. 
## Reads the dna counts csv and creates
## a pcounts fasta for the protein alignment usearch run
## Here, unique translations which appear fewer than 5 times are omitted
## for the LFC data sheets. For calculating frequencies, all reads
## are included.
for sample in samples:
    print(sample)
    path1 = './outputs/' + sample + '_dna_counts.csv'
    path2 = './outputs/' + sample + '_pcounts.fasta'
    data = pd.read_csv(path1,names=['qseq','reads'],skiprows=[0])
    print(data['reads'].sum())
    ## Skip the header row and redefine columns
    data['prot_seq'] = data.apply(translate,axis=1)
    print('prot_seqs')
    data = data.astype({'reads': int})
    prot_seqs = data['prot_seq'].value_counts()
    print('value_counted')
    prot_seqs = prot_seqs.reset_index()
    prot_seqs.rename(columns={'index':'prot_seq','prot_seq':'pcount'},inplace=True)
    # Only keep proteins which have 5 or more reads in total
    mixed_data = pd.merge(data,prot_seqs,how='left',on='prot_seq')
    mixed_data = mixed_data.groupby(by=['prot_seq']).sum()
    print('dont')
    
    ## this line is included for LFC calcs. Comment out for freq calcs
    #mixed_data = mixed_data[mixed_data['reads'] >= 5]
    mixed_data.drop('pcount',axis=1,inplace=True)
    mixed_data = pd.merge(mixed_data,prot_seqs,how='left',on='prot_seq')
    
    print('done')
    #drop pcount, add real pcounts
    
    prot_list = list(mixed_data['prot_seq'])
    pcount_list = list(mixed_data['pcount'])
    read_list = list(mixed_data['reads'])
    
    print('writing')
    for x in range(len(prot_list)):

        name = str(x) + '-' + str(pcount_list[x]) + '-' + str(read_list[x])
        seq = str(prot_list[x])
        with open(path2,'a') as output:
            output.write('>'+str(name)+'\n'+str(seq)+'\n')


## Step 7:
## Use USEARCH to get qrowdots output from pcounts fastas
## I used CMD in windows, requires usearch.exe in directory as well as
## the reference protein sequence as a fasta.

usearch --usearch_global R1_0_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R1_0_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq
usearch --usearch_global R1_375_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R1_375_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq
usearch --usearch_global R1_75_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R1_75_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq
usearch --usearch_global R2_0_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R2_0_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq
usearch --usearch_global R2_375_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R2_375_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq
usearch --usearch_global R2_75_pcounts.fasta --db ref_prot.fasta --id 0.0 --userout R2_75_palign.txt --userfields query+ql+id+mism+diffs+caln+qrowdots+qseq




## This interprets qrowdots into mutations, calculates read frequencies
## then creates the comparison files. Read frequences are based on the
## number of reads in each condition's file post-barcode sorting. 

dataframes = []    
for sample in samples:
    path1 = './outputs/' + sample + '_palign.txt'
    data = pd.read_csv(path1,sep='\t',names=['query','ql','id','mism','diffs','caln','qrowdots','qseq'])
    #alignments.to_csv(path2)
    print('reading palign file')
    data['sample'] = data.apply(extract1,axis=1)
    
    col1 = 'pcounts_' + sample
    col2 = 'reads_' + sample
    
    data[col1] = data.apply(extract2,axis=1)
    data[col2] = data.apply(extract3,axis=1)
    data = data.astype({col1:int,col2:int})
    data['mutations'] = data.apply(mutation_indexer,axis=1)
    dataframes.append(data)
    print('mutations applied')
    
    #check = data[data['qseq'] == 'MEMEKEFEQIDKSGSWAAIYQDIRHEASDFPCRVAKLPKNKNRNRYRDVSPFDHSRIKLHQEDNDYINASLIKMEEAQRSYILTQGPLPNTCGHFWEMVWEQKSRGVVMLNRVMEKGSLKCAQYWPQKEDKEMIFEDTNLKLTLISEDIKSYYTVRQLELENLTTQETREILHFHYTTWPDFGIPESPASFLNFLFKVRESGSLSPEHGPVVVHCSAGIGRSGTFCLADTCLLLMDKRKDPSSVDIKKVLLEMRKFRMGLIQTADQLRFSYLAVIEGAKFIMGDSSVQDQWKELSHEDLEPPPEHIPPPPRPPKRILEPHNLEHHHHHHAAADPMVRVLEA']

freq_dfs = []
idx = 0
for df in dataframes:
    col2 = 'reads_' + samples[idx]
    col3 = 'freq_' + samples[idx]
    print('calculating frequency')
    sample_reads = [667568,717434,780752,649875,738381,698341]
    df[col3] = df[col2] / sample_reads[idx]
    df = df.drop(labels=['query','caln','mism','sample'],axis=1)
    freq_dfs.append(df)
    idx += 1


print('creating comparison file')     
## first, combine the individual R1 dataframes with merge
R1_0vs375 = pd.merge(freq_dfs[0],freq_dfs[1],how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations'])
R1_0vs375['log2FC_0vs375'] = np.log2(R1_0vs375['freq_R1_375']/R1_0vs375['freq_R1_0'])

R1_0vs75 = pd.merge(freq_dfs[0],freq_dfs[2],how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations'])
R1_0vs75['log2FC_0vs75'] = np.log2(R1_0vs75['freq_R1_75']/R1_0vs75['freq_R1_0'])

R1_data = pd.merge(R1_0vs375,R1_0vs75,how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations','reads_R1_0','pcounts_R1_0','freq_R1_0'])
R1_data['mutation_count'] = R1_data.apply(mutation_count,axis=1)

## This code searches for duplicate mutation entries (confusing)
## then only keeps the one with the highest pcounts and read counts. 
## A couple of the output rows had the incorrect log2FC calculation, so 
## I recalculated frequencies and log2FC at the end.
dupes = R1_data[R1_data.duplicated(subset=['qseq']) == True]   
baddies = list(set(dupes['qseq'].values))
filler = pd.DataFrame()
for bad in baddies:
    temp = R1_data[R1_data['qseq'] == bad]
    fill = temp.max().to_frame().T
    filler = pd.concat([filler,fill],axis=0)
for index,row in filler.iterrows():
    if pd.isna(row['reads_R1_0']):
        continue
    ## These exceptiosn avoid issues when NAN values show up
    try:
        row['freq_R1_0'] = int(row['reads_R1_0'])/int(sample_reads[0])
    except: 
        pass
    try:
        row['freq_R1_375'] = int(row['reads_R1_375'])/int(sample_reads[1])
    except:
        pass
    try:
        row['freq_R1_75'] = int(row['reads_R1_75'])/int(sample_reads[2])
    except:
        pass
    row['log2FC_0vs375'] = np.log2(row['freq_R1_375']/row['freq_R1_0'])
    row['log2FC_0vs75'] = np.log2(row['freq_R1_75']/row['freq_R1_0'])
    
R1_data.drop_duplicates(subset=['qseq'],inplace=True,keep=False)
R1_data = pd.concat([R1_data,filler],axis=0)

# optional, remove rows where R1_0 has no reads
#R1_data = R1_data[R1_data['reads_R1_0'].notna()]


R1_data.to_csv('R1_data_freq.csv',index=False)

## Now for R2

print('creating comparison file')     

R2_0vs375 = pd.merge(freq_dfs[3],freq_dfs[4],how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations'])
R2_0vs375['log2FC_0vs375'] = np.log2(R2_0vs375['freq_R2_375']/R2_0vs375['freq_R2_0'])

R2_0vs75 = pd.merge(freq_dfs[3],freq_dfs[5],how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations'])
R2_0vs75['log2FC_0vs75'] = np.log2(R2_0vs75['freq_R2_75']/R2_0vs75['freq_R2_0'])

R2_data = pd.merge(R2_0vs375,R2_0vs75,how='outer',on=['qseq','ql','id','diffs','qrowdots','mutations','reads_R2_0','pcounts_R2_0','freq_R2_0'])
R2_data['mutation_count'] = R2_data.apply(mutation_count,axis=1)

## Addressing duplicate entries
dupes = R2_data[R2_data.duplicated(subset=['qseq']) == True]   
baddies = list(set(dupes['qseq'].values))
filler = pd.DataFrame()
for bad in baddies:
    temp = R2_data[R2_data['qseq'] == bad]
    fill = temp.max().to_frame().T
    filler = pd.concat([filler,fill],axis=0)
for index,row in filler.iterrows():
    if pd.isna(row['reads_R2_0']):
        continue
    try:
        row['freq_R2_0'] = int(row['reads_R2_0'])/int(sample_reads[0])
    except:
        pass
    try:
        row['freq_R2_375'] = int(row['reads_R2_375'])/int(sample_reads[1])
    except:
        pass
    try:
        row['freq_R2_75'] = int(row['reads_R2_75'])/int(sample_reads[2])
    except:
        pass
    row['log2FC_0vs375'] = np.log2(row['freq_R2_375']/row['freq_R2_0'])
    row['log2FC_0vs75'] = np.log2(row['freq_R2_75']/row['freq_R2_0'])
R2_data.drop_duplicates(subset=['qseq'],inplace=True,keep=False)
R2_data = pd.concat([R2_data,filler],axis=0)
#R2_data = R2_data[R2_data['reads_R2_0'].notna()]

R2_data.to_csv('R2_data_freq.csv',index=False)


## This part makes the extrapolation datasets. For simplicity, the assumed
## 4.99 reads does not contribute to the total reads in the frequency calculation.
data1 = pd.read_csv('R1_data_freq.csv')
data2 = pd.read_csv('R2_data_freq.csv')
sample_reads = [667568,717434,780752,649875,738381,698341]    

data1['reads_R1_0'].fillna(value=4.99,inplace=True)
data1['freq_R1_0'] = data1['reads_R1_0']/sample_reads[0]
data1['log2FC_0vs375'] = np.log2(data1['freq_R1_375']/data1['freq_R1_0'])
data1['log2FC_0vs75'] = np.log2(data1['freq_R1_75']/data1['freq_R1_0'])
data1.to_csv('R1_extrapolated_publish.csv',index=False)

data2['reads_R2_0'].fillna(value=4.99,inplace=True)
data2['freq_R2_0'] = data2['reads_R2_0']/sample_reads[3]
data2['log2FC_0vs375'] = np.log2(data2['freq_R2_375']/data2['freq_R2_0'])
data2['log2FC_0vs75'] = np.log2(data2['freq_R2_75']/data2['freq_R2_0'])
data2.to_csv('R2_extrapolated_publish.csv',index=False)


"""
## Now, looking at the frequency of mutations in variants with 5 
## or fewer mutations total

R1data = pd.read_csv('R1_data_freq.csv')
R1data = R1data[R1data['mutation_count'] <= 5]
#R1data = R1data[R1data['freq_R1_0'].notna()]
#R1data = R1data.dropna(subset=['freq_R1_0'])
R2data = pd.read_csv('R2_data_freq.csv')
R2data = R2data[R2data['mutation_count'] <= 5]
#R2data = R2data[R2data['freq_R2_0'].notna()]
#R2data = R2data[(R2data['reads_R2_375'].notna()) | (R2data['reads_R2_75'].notna())]
#R2data = R2data.dropna(subset=['freq_R2_0'])

mutations1 = R1data['mutations'].to_list()
mutations2 = R2data['mutations'].to_list()
mutations = mutations1 + mutations2
muts = []
## Runs through mutation list, adds each individual item to list
## Then duplicates are removed by converting to a set
for item in mutations:
    list_item = ast.literal_eval(item)
    for p in list_item:
        muts.append(p)
muts = list(set(muts))
## The asterix of stop codon mutations messes things up, delete from list
muts = [x for x in muts if '*' not in x]
print('unique mutations: ',len(muts))
R10_freqs = []
R1375_freqs = []
R175_freqs = []
R20_freqs = []
R2375_freqs = []
R275_freqs = []
position = []
for mut in muts:
    current1 = R1data[R1data['mutations'].str.contains(mut)]
    current2 = R2data[R2data['mutations'].str.contains(mut)]
    R10_freqs.append(current1['freq_R1_0'].sum())
    R1375_freqs.append(current1['freq_R1_375'].sum())
    R175_freqs.append(current1['freq_R1_75'].sum())
    R20_freqs.append(current2['freq_R2_0'].sum())
    R2375_freqs.append(current2['freq_R2_375'].sum())
    R275_freqs.append(current2['freq_R2_75'].sum())
    position.append(re.findall(r'\d+',mut)[0])
    
# Normalize to total frequency within a condition    
print('normalizing frequencies')
R10_norm = [i/sum(R10_freqs) for i in R10_freqs]
R1375_norm = [i/sum(R1375_freqs) for i in R1375_freqs]
R175_norm = [i/sum(R175_freqs) for i in R175_freqs]
R20_norm = [i/sum(R20_freqs) for i in R20_freqs]
R2375_norm = [i/sum(R2375_freqs) for i in R2375_freqs]
R275_norm = [i/sum(R275_freqs) for i in R275_freqs]

#LFC_R1_375 = np.log2(pd.Series(R1375_norm)/pd.Series(R10_norm)).to_list()
#LFC_R1_75 = np.log2(pd.Series(R175_norm)/pd.Series(R10_norm)).to_list()
#LFC_R2_375 = np.log2(pd.Series(R2375_norm)/pd.Series(R20_norm)).to_list()
#LFC_R2_75 = np.log2(pd.Series(R275_norm)/pd.Series(R20_norm)).to_list()

    
#final = pd.DataFrame(zip(muts,position,R10_norm,R1375_norm,R175_norm,R20_norm,R2375_norm,R275_norm,LFC_R1_375,LFC_R1_75,LFC_R2_375,LFC_R2_75),columns=['mutation','position','freq_R1_0','freq_R1_375','freq_R1_75','freq_R2_0','freq_R2_375','freq_R2_75','LFC_R1_375','LFC_R1_75','LFC_R2_375','LFC_R2_75'])
final2 = pd.DataFrame(zip(muts,position,R10_norm,R1375_norm,R175_norm,R20_norm,R2375_norm,R275_norm),columns=['mutation','position','freq_R1_0','freq_R1_375','freq_R1_75','freq_R2_0','freq_R2_375','freq_R2_75'])

R10_sum = []
R1375_sum = []
R175_sum = []
R20_sum = []
R2375_sum = []
R275_sum = []
pos = []
for x in range(0,329):
    current = final2[final2['position'] == str(x+1)]
    R10_sum.append(current['freq_R1_0'].sum())
    R1375_sum.append(current['freq_R1_375'].sum())
    R175_sum.append(current['freq_R1_75'].sum())
    R20_sum.append(current['freq_R2_0'].sum())
    R2375_sum.append(current['freq_R2_375'].sum())
    R275_sum.append(current['freq_R2_75'].sum())
    pos.append(x+1)

LFC_R1_375 = np.log2(pd.Series(R1375_sum)/pd.Series(R10_sum)).to_list()
LFC_R1_75 = np.log2(pd.Series(R175_sum)/pd.Series(R10_sum)).to_list()
LFC_R2_375 = np.log2(pd.Series(R2375_sum)/pd.Series(R20_sum)).to_list()
LFC_R2_75 = np.log2(pd.Series(R275_sum)/pd.Series(R20_sum)).to_list()

final3 = pd.DataFrame(zip(pos,R10_sum,R1375_sum,R175_sum,R20_sum,R2375_sum,R275_sum,LFC_R1_375,LFC_R1_75,LFC_R2_375,LFC_R2_75),columns=['position','total_R1_0','freq_R1_375','total_R1_75','total_R2_0','total_R2_375','total_R2_75','LFC_R1_375','LFC_R1_75','LFC_R2_375','LFC_R2_75'])
final3.to_csv('position_frequency_all_reads.csv',index=False)
#final.to_csv('mutation_frequency.csv',index=False)







