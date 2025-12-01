Scripts to generate .defattr files for chimeraX given a .cif structure

## System requirements
SO: Tested on Red Hat Enterprise Linux 8
CPU: Tested on Intel and AMD CPUs
RAM: At least 4GB

## Dependencies
Note: for AF3 screen to structure, only numpy and biopython are needed

    numpy
    biopython
    polars
    pyfaidx

## Installation
Activate or set up a conda environment with the above packages

    conda create -n write_defattr numpy biopython pandas polars pyranges pyfaidx
    conda activate write_defattr


Navigate to the folder where you can to clone this repo and do

    git clone https://github.com/JoseFcoMH/defattr_painters
    cd defattr_painters


Should be done in a few minutes.    

## Usage

### Usage of AF3screenToStructure.py
    
    usage: Process AF3 inference. [-h] [--runName RUNNAME] --baseDir BASEDIR --outDir OUTDIR --cifFile CIFFILE [--save_processed_pairs]
    
    Problems? Ask Jose.
    
    optional arguments:
      -h, --help            show this help message and exit
      --runName RUNNAME     Name of this run. Can be whatever.
      --baseDir BASEDIR     Directory containing the raw AF3 inference output.
      --outDir OUTDIR       Directory where to save the output.
      --cifFile CIFFILE     Full path to .cif file to paint.
      --save_processed_pairs
                            Whether to save the processed pairs in outDir
    

outDir is the folder where the output, a list of processed af3_pair objects and the .defattr file will be saved to.
The {runName}\_{other} must follow the default local AF3 output file structure.

### Usage of transcriptToStructure.py
    
    usage: Map bed-like tracks to the RNA in some structure. [-h] --inputFiles INPUTFILES --inputStrands INPUTSTRANDS --genesOfInterest GENESOFINTEREST --outDir OUTDIR --runName RUNNAME --gtf GTF --fasta FASTA --cifFile
                                                             CIFFILE [--skipCheck]
    
    Problems? Ask Jose.
    
    options:
      -h, --help            show this help message and exit
      --inputFiles INPUTFILES
                            Comma-separated paths to input files. Each must hace 3 columns [Chromosome, Position, Score]
      --inputStrands INPUTSTRANDS
                            Comma-separated strand where the values of the input files are, in the same order as the input files. + or -
      --genesOfInterest GENESOFINTEREST
                            Comma-separated names of the gene to look at, as per the GTF file. E.g. RNA18SN1,RNA18SN2
      --outDir OUTDIR       Directory to save the output .defattr in.
      --runName RUNNAME     Unique name for this run.
      --gtf GTF             Path to a GTF file. Assumes columns ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] like in the RefSeq GTF.
      --fasta FASTA         Path to a genome fasta file.
      --cifFile CIFFILE     Full path to .cif file to paint.
      --skipCheck           Whether to assume that RNA chains have a description containing the keyword "RNA" for speed.
    


### Example usage
Navigate to the cloned repository, unzip the test_data.zip folder, and run

    python AF3screenToStructure.py --runName zak --baseDir test_data/AF3_output  --outDir . --cifFile test_data/4ug0.cif

Which should finish in a few seconds. 

Or run
    
    python transcriptToStructure.py --inputFiles test_data/RNA_crosslinks/exampleFile.tsv \
                                    --inputStrands + \
                                    --genesOfInterest RNA18SN1 \
                                    --outDir . \
                                    --gtf /path/to/gtf \
                                    --fasta /path/to/fasta \
                                    --runName test_run \
                                    --cifFile test_data/4ug0.cif \
                                    --skipCheck

Which should also finish in a few seconds.

## References

[AlphaFold3_Repo](https://github.com/google-deepmind/alphafold3)
[AlphaFold3_Publication](https://www.nature.com/articles/s41586-024-07487-w)
[ChimeraX](https://onlinelibrary.wiley.com/doi/10.1002/pro.4792)

## Citation
[ZAKα is a sensor of mRNA stasis at the ribosomal exit channel]https://www.biorxiv.org/content/10.1101/2025.11.22.689755v1





