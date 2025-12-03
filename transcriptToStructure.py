import numpy as np
import polars as pl
import pyfaidx

from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import PairwiseAligner

from concurrent.futures import ProcessPoolExecutor
import os
import argparse
from pathlib import Path
import logging



def extract_sequences_from_cif(cif_file):
    cif_dict = MMCIF2Dict(str(cif_file))
    seqs = [x for x in cif_dict['_entity_poly.pdbx_seq_one_letter_code']]
    ids = [x for x in cif_dict['_entity_poly.pdbx_strand_id']]
    descs = [x for x in cif_dict['_entity.pdbx_description']]
    sequences = [SeqRecord(Seq(seq.replace('\n', '').upper()), id=i, description=d) for seq, i, d in zip(seqs, ids, descs)]
            
    return sequences

# Function to extract sequence given coord_idxnates
def get_sequence(chromosome, start, end, strand='+', genome:pyfaidx.Fasta|str|Path =''):
    """
    Extract sequence from genome.
    
    Parameters:
    chromosome : str
        Chromosome name (e.g., 'chr1')
    start : int
        0-based start position
    end : int
        1-based end position (exclusive)
    strand : str
        '+' or '-' for strand
    genome
        genome or path to a genome
    
    Returns:
    str : The requested sequence
    """
    if isinstance(genome, str) or isinstance(genome, Path):
        genome = pyfaidx.Fasta(str(genome))
        
    try:
        seq = genome[chromosome][start:end].seq
        
        # Reverse complement if on minus strand
        if strand == '-':
            seq = pyfaidx.complement(seq)[::-1]
            
        return seq
    except (KeyError, ValueError) as e:
        logging.info(f"Error fetching sequence at {chromosome}:{start}-{end}: {e}")
        return ""


class iclip_track():
    
    def __init__(self, gene, crosslinks_df, gtf_df, fasta_path, mock=False):
        if mock:
            return
        
        self.gene_name = gene
        self.fasta_path = fasta_path
        
        gene_df = gtf_df.filter(
                (pl.col('gene_id') == gene) & 
                pl.col('feature').str.contains('transcript')
            ).with_columns(
                pl.col('start').cast(pl.Int32) - 1, # convert to 0-based coords
                pl.col('end').cast(pl.Int32)
            )
            
        self.chromosome = gene_df['seqname'][0]
        self.start = gene_df['start'][0]
        self.end = gene_df['end'][0]
        self.strand = gene_df['strand'][0]
        self.len = self.end - self.start
        
        self.sequence = get_sequence(
            chromosome=self.chromosome,
            start=self.start,
            end=self.end,
            strand=self.strand, 
            genome=fasta_path
        ).upper().replace('T', 'U')
        
        crosslinks_at_gene = crosslinks_df.filter(
            (pl.col('chromosome') == self.chromosome) &
            (pl.col('crosslink_site') >= self.start) &
            (pl.col('crosslink_site') < self.end) &
            (pl.col('strand') == self.strand)
        )
        self.crosslinks = np.zeros(self.len)
        self.crosslinks[crosslinks_at_gene['crosslink_site'] - self.start] = crosslinks_at_gene['n_crosslinks'] 
        if self.strand == '-':
            self.crosslinks = self.crosslinks[::-1]
        
        assert len(self.crosslinks) == len(self.sequence)
        
        self.best_alignment = None
        self.best_score = -np.Inf
        self.best_seq = None
      
    
    @staticmethod
    def get_best_alignment(chain: str, seq: SeqRecord):
        aligner = trackwiseAligner()
        aligner.mismatch_score = -np.inf # Just put a gap instead
        aligner.match_score = 1
        try:
            alignment = next(aligner.align(chain, seq))
            if alignment.score > (len(chain) * 0.8): # Arbitrarily chosen by none other than myself
                return alignment.score, seq, alignment.aligned
            else:
                return -np.inf, seq, ((), ())
        except:
            return -np.inf, seq, ((), ())
    
    
    @staticmethod
    def get_alignment_wrapper(args):
        chain, seq = args
        return iclip_track.get_best_alignment(chain, seq)
    

    def map_interactor(self, L, n_cpus=1):
        best_score = self.best_score
        best_seq = self.best_seq
        best_alignment = self.best_alignment
                
        args_list = [(self.sequence, seq) for seq in L]

        with ProcessPoolExecutor(max_workers=n_cpus) as executor:
            results = list(executor.map(iclip_track.get_alignment_wrapper, args_list))

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
        
        if self.best_alignment is None and L is not None: self.map_interactor(L)
        if self.best_seq is None: 
            logging.info(f'No mappable attribute for {self.name}')
            return self
        
        score_new = np.zeros(len(self.best_seq.seq))
        score_new[:] = np.nan
        for seqA_chunk, seqB_chunk in zip(self.best_alignment[0], self.best_alignment[1]):
            score_new[seqB_chunk[0]:seqB_chunk[1]] = self.crosslinks[seqA_chunk[0]:seqA_chunk[1]]
        self.crosslinks_attr = score_new
        
        return self
    
    
    def get_attr(self):
        return self.crosslinks_attr
    
    
    def get_match(self):
        if self.best_seq is not None:
            return self.best_seq.id, self.best_score
        return None, self.best_score
    
    
    def get_match_id(self):
        if self.best_seq is not None:
            return self.best_seq.id
        return None
    
    
    def get_match_score(self):
        return self.best_score
    
    
def main():
    
    n_cpus = (
        int(os.environ.get("SLURM_CPUS_ON_NODE", 0)) or
        int(os.environ.get("SLURM_CPUS_PER_TASK", 0)) or
        int(os.environ.get("PBS_NUM_PPN", 0)) or
        os.cpu_count()  # fallback to all available CPUs and hope for the best
    )
    
    logging.info(f"Using {n_cpus} CPUs.")
    
    
    parser = argparse.ArgumentParser(
        prog='Map bed-like tracks to the RNA in some structure.',
        description='Problems? Ask Jose.',
        epilog='Hope it works.'
    )

    parser.add_argument(
        '--inputFiles',
        required=True,
        help='Comma-separated paths to input files. Each must hace 3 columns [Chromosome, Position, Score]'
    )
    parser.add_argument(
        '--inputStrands',
        required=True,
        help='Comma-separated strand where the values of the input files are, in the same order as the input files. + or -'
    )
    parser.add_argument(
        '--genesOfInterest',
        required=True,
        help='Comma-separated names of the gene to look at, as per the GTF file. E.g. RNA18SN1,RNA18SN2'
    )
    parser.add_argument(
        '--outDir',
        required=True, type=Path, 
        help='File to save the output .defattr to.'
    )
    parser.add_argument(
        '--runName',
        required=True, type=str, 
        help='Unique name for this run.'
    )
    parser.add_argument(
        '--gtf', type=Path,
        required=True,
        help="Path to a GTF file. Assumes columns ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] like in the RefSeq GTF."
    )
    parser.add_argument(
        '--fasta', type=Path,
        required=True,
        help='Path to a genome fasta file.'
    )
    parser.add_argument(
        '--cifFile', type=Path,
        required=True,
        help='Full path to .cif file to paint.'
    )
    parser.add_argument(
        '--skipCheck',
        action='store_true',
        help='Whether to assume that RNA chains have a description containing the keyword "RNA" for speed.'
    )

    in_args = parser.parse_args()

    input_beds = in_args.inputFiles.split(',')
    input_strands = in_args.inputStrands.split(',')
    
    genesOfInterest = in_args.genesOfInterest.split(',')
    
    cif_file = in_args.cifFile
    cif_name = cif_file.stem
    fasta_path = pyfaidx.Fasta(in_args.fasta)    
    gtf_file = in_args.gtf
        
    outDir = in_args.outDir
    outDir.mkdir(parents=True, exist_ok=True)
    runName = in_args.runName
    
    logging.info(f'Loading .cif file at {cif_file}')
    refSeqs = extract_sequences_from_cif(cif_file)

    if in_args.skipCheck:
        refSeqsRNA = refSeqs
    else:
        refSeqsRNA = [s for s in refSeqs if 'RNA' in s.description]
    
        
    crosslinks_df = pl.DataFrame({
        'chromosome': pl.Series(pl.utf8),
        'crosslink_site': pl.Series(pl.int64),
        'n_crosslinks': pl.Series(pl.int64),
        'strand': pl.Series(pl.utf8)   
    })
    for input_bed, strand in zip(input_beds, input_strands):
        logging.info(f'Loading input file at {input_bed}, {strand} strand')
        
        temp_bed = pl.read_csv(input_bed, separator='\t', has_header=False) 
        temp_bed = (
            temp_bed
            .rename({old:new for old, new in zip(temp_bed.columns, ['chromosome', 'crosslink_site', 'n_crosslinks'])})
            .with_columns(pl.lit(strand).alias('strand'))
        )
        crosslinks_df = pl.concat([crosslinks_df, temp_bed])
    
    logging.info(f'Parsing gtf {gtf_file}')
    column_names = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    dtypes = {nam:pl.Utf8 for nam in column_names}

    gtf_df = pl.read_csv(gtf_file, separator='\t' , comment_prefix='#', has_header=False, new_columns=column_names, schema_overrides=dtypes)
    gtf_df = gtf_df.with_columns(
        pl.col('attribute').map_elements(lambda x: x.split('gene_id "')[-1].split('"')[0], return_dtype=str).alias('gene_id'),
        pl.col('attribute').map_elements(lambda x: x.split('gene_name "')[-1].split('"')[0], return_dtype=str).alias('gene_name'),
        pl.col('attribute').map_elements(lambda x: x.split('biotype "')[-1].split('"')[0], return_dtype=str).alias('biotype')
    )

    iclip_tracks: list[iclip_track] = []   
    for geneOI in genesOfInterest:
        logging.info(f'Calculating attribute for transcript {geneOI}')
        iclip_tracks.append(
            iclip_track(
                gene=geneOI,
                crosslinks_df=crosslinks_df,
                gtf_df=gtf_df,
                fasta_path=fasta_path
            )
            .map_interactor(refSeqsRNA, n_cpus)
            .calc_attr()
        )
        
    logging.info(f'Done. Got {len(iclip_tracks)} genes of interest, {sum(track.best_seq is not None for track in iclip_tracks)} of which mapped to a chain in {cif_name}. Creating attribute file...')
    
    # Messy chunks. Maybe I'll fix someday
    best_dict = {chain:{'gene_name':'', 'score':0, 'track':None} for chain in set([x.get_match_id() for x in iclip_tracks])}
    for track in iclip_tracks:
        match_chain_id = track.get_match_id()
        current_match_score = track.get_match_score()
        if current_match_score > best_dict[match_chain_id]['score']:
            best_dict[match_chain_id] = {
                'gene_name': track.gene_name,
                'score': current_match_score,
                'track': track
            }
            
    assert(len(set(best_dict.keys())) == len(set([x['gene_name'] for x in best_dict.values()])))

    best_tracks: list[iclip_track|None] = [x['track'] for x in best_dict.values()]
                
    with open(outDir / f'{runName}_{cif_name}.defattr', 'w+') as f:
        f.write(f'attribute: crosslinks_{runName}_{cif_name}\nmatch mode: any\nrecipient: residues\n')
        for track in best_tracks:
            if track is not None:
                crosslinks_vals = track.get_attr()
                match_chain_id = track.get_match_id()
                for i in range(len(crosslinks_vals)):
                    if np.isnan(crosslinks_vals[i]): 
                        continue
                    f.write(f'\t/{match_chain_id}:{i+1}\t{crosslinks_vals[i]}\n')
    
    logging.info(f'Wrote {outDir}/{runName}_{cif_name}.defattr')
    

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(message)s"
    )
    
    logging.info('Mapping crosslinks to structure...')
    
    main()
    
    logging.info('All done!')

