import numpy as np
import json
import os
import glob
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import pairwise2
from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import pickle as pkl
import argparse

parser = argparse.ArgumentParser(
    prog='Process AF3 inference. Asumes a directory structure of {baseDir}/{runName}_{other}, where {runName}_{other} are the {refProtein}_{possibleInteractor} pairs.',
    description='Problems? Ask Jose.',
    epilog='Hope it works.'
)

parser.add_argument('--runName', required=True, help='Name of run and ref protein to process.')
parser.add_argument('--baseDir', required=False, help='Parent directory holding the AF3 inference output.')
parser.add_argument('--outDir', required=False, help='Directory to save the output.')
parser.add_argument('--cifFile', required=False, help='Full path to .cif file o paint.', default='/home/vzl634/datadir/structures/4ug0.cif')

in_args = parser.parse_args()

if not in_args.baseDir:
    in_args.baseDir = f'/home/vzl634/datadir/af3_runs/{in_args.runName}_run/inference_output'

if not in_args.outDir:
    in_args.outDir = os.path.abspath(os.path.join(in_args.baseDir, os.pardir))
    
runName = in_args.runName
baseDir = in_args.baseDir
cif_file = in_args.cifFile


def extract_sequences_from_cif(cif_file):
    cif_dict = MMCIF2Dict(cif_file)
    seqs = [x for x in cif_dict['_entity_poly.pdbx_seq_one_letter_code']]
    ids = [x for x in cif_dict['_entity_poly.pdbx_strand_id']]
    descs = [x for x in cif_dict['_entity.pdbx_description']]
    sequences = [SeqRecord(Seq(seq.replace('\n', '').upper()), id=i, description=d) for seq, i, d in zip(seqs, ids, descs)]
            
    return sequences

class af3_pair():
    def __init__(self, pred_dir):
        self.pred_dir = pred_dir
        self.structure = MMCIF2Dict(f'{pred_dir}/{pred_dir.split("/")[-1]}_model.cif')
        self.interactor = pred_dir.split('_')[-1]
        self.min_paes = []
        self.iptms = []
        self.ptms = []
        self.chains = []
        self.confidences = None
        self.iCP = None
        self.best_score = -np.inf
        self.best_seq = None
        self.best_alignment = None
        self.summaries = glob.glob(f'{pred_dir}/*/summary_confidences.json')
        self.iCP_attr = None
        
        with concurrent.futures.ThreadPoolExecutor() as executor:
            results = list(executor.map(self.process_summary, self.summaries))
        
        for result in results:
            min_pae, iptm, ptm = result
            self.min_paes.append(min_pae)
            self.iptms.append(iptm)
            self.ptms.append(ptm)
        
        self.mean_iptm = np.mean(self.iptms)
        self.min_pae = np.min(self.min_paes)
    
    def process_summary(self, summary):
        with open(summary, 'r') as f:
            result = json.load(f)
            min_pae = result['chain_pair_pae_min'][0][1]
            iptm = result['iptm']
            ptm = result['ptm']
        return min_pae, iptm, ptm
        
    def load_chains(self):
        cif_dict = self.structure

        three_to_one = {
            'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
            'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
            'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
            'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y',
            'SEC': 'U', 'PYL': 'O', 'ASX': 'B', 'GLX': 'Z', 'XAA': 'X',
            'C':'C', 'U':'U', 'T':'T', 'G':'G', 'A':'A'
        }
                
        chainID = [int(x) for x in cif_dict['_entity_poly_seq.entity_id']]
        aas = [three_to_one[x] for x in cif_dict['_entity_poly_seq.mon_id']]
        
        chainID.append('')
        aas.append('')
        
        i0 = chainID[0]
        currentChain = ''
        for i, aa in zip(chainID, aas):
            if i == i0:
                currentChain += aa
            else:
                self.chains.append(currentChain)
                currentChain = aa
                i0 = i
        
        return self
    
    def calc_iCP(self):
        with open(f'{self.pred_dir}/{self.pred_dir.split("/")[-1]}_confidences.json', 'r') as f:
            self.confidences = json.load(f)
            chain_len = {chain:self.confidences['token_chain_ids'].count(chain) for chain in ['A', 'B']}
            contact_probs = (np.array(self.confidences['contact_probs'])[chain_len['A']:, :chain_len['A']] + np.array(self.confidences['contact_probs'])[:chain_len['A'], chain_len['A']:].T) / 2
            self.iCP = contact_probs.max(1)
        
        return self
        
    def map_interactor(self, L):
        
        for seq in L:
            alignments = pairwise2.align.globalms(self.chains[1], seq.seq, 2, -10, -3, -0.1)
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
            alignments = pairwise2.align.globalms(self.chains[1], seq.seq, 2, -10, -3, -0.1)
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
        if self.iCP is None: self.calc_iCP()
        
        score_new = np.zeros(len(self.best_alignment.seqA))
        score_new[:] = np.nan
        score_index = 0
        for i in range(len(score_new)):
            if self.best_alignment.seqA[i] == '-': continue
            score_new[i] = self.iCP[score_index]
            score_index += 1
        
        self.iCP_attr =  np.array([x for x, y in zip(score_new, self.best_alignment.seqB) if y != '-'])
        
        return self
    
    def get_attr(self):
        return self.iCP_attr
        
    def get_iptm(self):
        return self.mean_iptm
    
    def get_min_pae(self):
        return self.min_pae
    
    def get_interactor(self):
        return self.interactor
    
    def get_match(self):
        return self.best_seq.id, self.best_score
    
    
    
refSeqs = extract_sequences_from_cif(cif_file)
cif_name = cif_file.split('/')[-1].split('.')[0]

results = glob.glob(f'{baseDir}/{runName}_*')
results.sort()

af3_pairs = [af3_pair(pair).load_chains().map_interactor2(refSeqs).calc_attr() for pair in results if os.path.exists(f'{pair}/ranking_scores.csv')]

with open(f'{outDir}/{runName}_mappedPairs.pkl', 'wb') as f:
    pkl.dump(af3_pairs, f)
    
best_dict = {chain:('', -np.inf) for chain in set([x.get_match()[0] for x in af3_pairs])}
for pair in af3_pairs:
    if pair.get_match()[1] > best_dict[pair.get_match()[0]][1]:
        best_dict[pair.get_match()[0]] = (pair.get_interactor(), pair.get_match()[1], pair)
        
assert(len(set(best_dict.keys())) == len(set([x[0] for x in best_dict.values()])))

best_pairs = [x[2] for x in best_dict.values()]

with open(f'{outDir}/{runName}_{cif_name}.defattr', 'w+') as f:
    f.write(f'attribute: iCP_{runName}\nmatch mode: any\nrecipient: residues\n')
    for pair in best_pairs:
        for i in range(len(pair.get_attr())):
            if np.isnan(pair.get_attr()[i]): continue
            f.write(f'\t/{pair.get_match()[0]}:{i+1}\t{pair.get_attr()[i]}\n')