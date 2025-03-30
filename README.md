Scripts to generate .defattr files for chimeraX given a .cif structure

## System requirements
    SO: Tested on Red Hat Enterprise Linux 8
    CPU: Tested on Intel and AMD CPUs
    RAM: At least 4GB

## Dependencies
    numpy
    biopython
    pandas
    polars
    pyranges
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
    
    usage: Process AF3 inference. Asumes a directory structure of {baseDir}/{runName}_{other}, where {runName}_{other} are the {refProtein}_{possibleInteractor} pairs. [-h] --runName RUNNAME [--baseDir BASEDIR] [--outDir OUTDIR] [--cifFile CIFFILE]
    
    Problems? Ask Jose.
    
    optional arguments:
      -h, --help         show this help message and exit
      --runName RUNNAME  Name of run and ref protein to process.
      --baseDir BASEDIR  Parent directory holding the AF3 inference output.
      --outDir OUTDIR    Directory to save the output to.
      --cifFile CIFFILE  Full path to .cif file o paint.
    

outDir is the folder where the output, a list of processed af3_pair objects and the .defattr file will be saved to.
The {runName}\_{other} must follow the default local AF3 output file structure.

### Usage of transcriptToStructure.py
    
    usage: Map bed-like tracks to the RNA in some structure. [-h] --inputFile INPUTFILE --inputStrand INPUTSTRAND --geneOfInterest GENEOFINTEREST [--outFile OUTFILE] [--gtf GTF] [--fasta FASTA] [--cifFile CIFFILE] [--skipCheck SKIPCHECK]
    
    Problems? Ask Jose.
    
    optional arguments:
      -h, --help            show this help message and exit
      --inputFile INPUTFILE
                            Path to input file. Must hace 3 columns [Chromosome, Position, Score]
      --inputStrand INPUTSTRAND
                            Strand where the values of the input file are. + or -
      --geneOfInterest GENEOFINTEREST
                            Name of the gene to look at, as per the GTF file. E.g., RNA18SN1
      --outFile OUTFILE     File to save the output .defattr to.
      --gtf GTF             Path to a GTF file. Assumes columns ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'] like in the RefSeq GTF.
      --fasta FASTA         Path to a genome fasta file.
      --cifFile CIFFILE     Full path to .cif file to paint.
      --skipCheck SKIPCHECK
                            Whether to assume that RNA chains have a description containing the keyword "RNA" for speed.


### Example usage
Navigate to the cloned repository and run

    python AF3screenToStructure.py --runName zak --baseDir test_data/AF3_output  --outDir . --cifFile test_data/4ug0.cif

Which should finish in a few seconds. 

    Or run
    
    python transcriptToStructure.py --inputFile test_data/RNA_crosslinks/exampleFile.tsv \
                                    --inputStrand + \
                                    --geneOfInterest RNA18SN1 \
                                    --outFile RNA_crosslinks.defattr \
                                    --gtf /path/to/gtf \
                                    --fasta /path/to/fasta \
                                    --cifFile test_data/4ug0.cif

Which should also finish in a few seconds.

## References

[AlphaFold3_Repo](https://github.com/google-deepmind/alphafold3)
[AlphaFold3_Publication](https://www.nature.com/articles/s41586-024-07487-w)
[PyRanges_Repo](https://github.com/pyranges/pyranges)
[Pyranges_Publication](https://academic.oup.com/bioinformatics/article/36/3/918/5543103)
[ChimeraX](https://onlinelibrary.wiley.com/doi/10.1002/pro.4792)
