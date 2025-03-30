import polars as pl
import pandas as pd
import glob
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import pairwise2
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import pyranges as pr
import pyfaidx
import numpy as np
import os
import argparse


parser = argparse.ArgumentParser(
    prog='Map bed-like tracks to the RNA in some structure.',
    description='Problems? Ask Jose.',
    epilog='Hope it works.'
)

parser.add_argument('--inputFile', required=True, help='Path to input file. Must hace 3 columns [Chromosome, Position, Score]')
parser.add_argument('--inputStrand', required=True, help='Strand where the values of the input file are. + or -')
parser.add_argument('--geneOfInterest', required=True, help='Name of the gene to look at, as per the GTF file. E.g., RNA18SN1')
parser.add_argument('--outFile', required=False, help='File to save the output .defattr to.')
parser.add_argument('--gtf', required=True, help="Path to a GTF file. Assumes columns ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] like in the RefSeq GTF.")
parser.add_argument('--fasta', required=True, help='Path to a genome fasta file.')
parser.add_argument('--cifFile', required=True, help='Full path to .cif file to paint.')
parser.add_argument('--skipCheck', required=False, help='Whether to assume that RNA chains have a description containing the keyword "RNA" for speed.', default=False)

in_args = parser.parse_args()

input_bed = in_args.inputFile
input_strand = in_args.inputStrand
cif_file = in_args.cifFile
fasta_path = in_args.fasta
gtf_file = in_args.gtf
geneOfInterest = in_args.geneOfInterest
if not in_args.outFile:
    in_args.outFile = os.path.abspath(os.path.join(in_args.inputFile, os.pardir)) + f"{geneOfInterest}_{input_bed.split('/')[-1].split('.')[0]}_{cif_file.split('/')[-1].split('.')[0]}.defattr"
outFile = in_args.outFile

def extract_sequences_from_cif(cif_file):
    cif_dict = MMCIF2Dict(cif_file)
    seqs = [x for x in cif_dict['_entity_poly.pdbx_seq_one_letter_code']]
    ids = [x for x in cif_dict['_entity_poly.pdbx_strand_id']]
    descs = [x for x in cif_dict['_entity.pdbx_description']]
    sequences = [SeqRecord(Seq(seq.replace('\n', '').upper()), id=i, description=d) for seq, i, d in zip(seqs, ids, descs)]
            
    return sequences

print(f'Loading .cif file at {cif_file}', flush=True)
refSeqs = extract_sequences_from_cif(cif_file)

if in_args.skipCheck:
    refSeqsRNA = refSeqs
else:
    refSeqsRNA = [s for s in refSeqs if 'RNA' in s.description]
    
print(f'Loading input file at {input_bed}', flush=True)
crosslinks_bed = pr.read_bed(input_bed).df
crosslinks_bed['Score'] = crosslinks_bed['End']
crosslinks_bed['Strand'] = input_strand
crosslinks_bed['End'] = crosslinks_bed['Start'] + 1
crosslinks_bed = pr.PyRanges(crosslinks_bed)


print(f'Parsing gtf {gtf_file}', flush=True)
column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
dtypes = {nam:pl.Utf8 for nam in column_names}

gtf_df = pl.read_csv(gtf_file, separator='\t' , comment_prefix='#', has_header=False, new_columns=column_names, schema_overrides=dtypes)
gtf_df = gtf_df.with_columns(
    pl.col('attribute').map_elements(lambda x: x.split('gene_id "')[-1].split('"')[0], return_dtype=str).alias('gene_id'),
    pl.col('attribute').map_elements(lambda x: x.split('gene_name "')[-1].split('"')[0], return_dtype=str).alias('gene_name'),
     pl.col('attribute').map_elements(lambda x: x.split('biotype "')[-1].split('"')[0], return_dtype=str).alias('biotype')
)


class iclip_track():
    
    def __init__(self, gene, crosslinks_bed=crosslinks_bed, fasta_path=fasta_path):
        self.gene_name = gene
        self.fasta_path = fasta_path
        self.crosslinks_bed = crosslinks_bed
        self.gene = gtf_df.filter((pl.col('gene_id') == gene) & pl.col('feature').str.contains('transcript')
                                ).rename(
                                    {'seqname':'Chromosome', 'start':'Start', 'end':'End', 'strand':'Strand'}
                                ).with_columns(
                                    pl.col('Start').cast(pl.Int32) - 1,
                                    pl.col('End').cast(pl.Int32)
                                )
        self.len = self.gene['End'][0] - self.gene['Start'][0]
        self.crosslinks = np.zeros(self.len)
        self.gene_range = pr.PyRanges(self.gene.to_pandas())
        self.score_seq(gtf_df, fasta_path)
        self.best_alignment = None
        self.best_score = -np.Inf
        self.best_seq = None
        
    def score_seq(self, gtf_df=gtf_df, fasta=fasta_path):
        self.sequence = pr.get_sequence(self.gene_range, fasta)[0].upper().replace('T', 'U')
        
        crosslinks_remap = self.crosslinks_bed[(self.crosslinks_bed.Start >= self.gene['Start'][0]) & (self.crosslinks_bed.End < self.gene['End'][0]) & (self.crosslinks_bed.Chromosome == self.gene['Chromosome'][0]) & (self.crosslinks_bed.Strand == self.gene['Strand'][0])]

        for idx, cross in crosslinks_remap.df.iterrows():
            self.crosslinks[cross.Start - self.gene['Start'][0]] = cross.Score
    
    def get_cross_seq(self):
        return self.sequence, self.crosslinks
    
    
    def map_interactor(self, L):
        
        for seq in L:
            alignments = pairwise2.align.globalms(self.sequence, seq.seq, 2, -5, -3, -0.1)
            score = alignments[0][2]

            if score > self.best_score:
                self.best_score = score
                self.best_seq = seq
                self.best_alignment = alignments[0]
            
        return self

    def map_interactor2(self, L):
        best_score = self.best_score
        best_seq = self.best_seq
        best_alignment = self.best_alignment

        def get_best_alignment(seq):
            alignments = pairwise2.align.globalms(self.sequence, seq.seq, 2, -5, -3, -0.1)
            score = alignments[0][2]
            return score, seq, alignments[0]

        with ThreadPoolExecutor() as executor:
            results = executor.map(get_best_alignment, L)

            for score, seq, alignment in results:
                if score > best_score:
                    best_score = score
                    best_seq = seq
                    best_alignment = alignment

        self.best_score = best_score
        self.best_seq = best_seq
        self.best_alignment = best_alignment

        return self
    
    def calc_attr(self, L=None):
        
        if self.best_alignment is None or L is not None: self.map_interactor2(L)
        
        score_new = np.zeros(len(self.best_alignment.seqA))
        score_index = 0
        for i in range(len(score_new)):
            if self.best_alignment.seqA[i] == '-': continue
            score_new[i] = self.crosslinks[score_index]
            score_index += 1
        
        self.cross_attr = np.array([x for x, y in zip(score_new, self.best_alignment.seqB) if y != '-'])
        
        return self 
    
    def write_attr(self, fileName, extraLab=''):
        
        with open(fileName, 'w+') as f:
            f.write(f'attribute: RNAscore_{self.gene_name}_{extraLab}\nmatch mode: any\nrecipient: residues\n')
            for i in range(self.len):
                f.write(f'\t/{self.best_seq.id}:{i+1}\t{self.crosslinks[i]}\n')
            
print(f'Calculating attribute for transcript {geneOfInterest}', flush=True)
t = iclip_track(geneOfInterest).calc_attr(refSeqsRNA)

t.write_attr(f'{outFile}', extraLab=input_bed.split('/')[-1].split('.')[0])
