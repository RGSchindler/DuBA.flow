# Readme
Software pipeline for nanopore sequencing

## Usage
This program was tested on Ubuntu 20.04 & Python ver. 3.8.1 with following packages:

Tool | Version | Installation 
-----|---------|-------------
[pyYAML](https://github.com/yaml/pyyaml) | 6.0 | pip install pyYAML
[samtools](https://github.com/samtools/samtools) | 1.9 | conda install -c bioconda samtools=1.9
[minimap2](https://github.com/lh3/minimap2) | 2.24 | conda install -c bioconda minimap2=2.24
[ngmlr](https://github.com/philres/ngmlr) | 0.2.7 | conda install -c bioconda ngmlr=0.2.7
[canu](https://github.com/marbl/canu) | 1.9 | conda install -c bioconda canu=1.9
[longshot](https://github.com/pjedge/longshot) | 0.4.1 | conda install -c bioconda longshot=0.4.1
[Sniffles](https://github.com/fritzsedlazeck/Sniffles) | 1.0.12 | conda install -c bioconda sniffles=1.0.12
[NanoPlot](https://github.com/wdecoster/NanoPlot) | 1.39.0 | pip install NanoPlots
[minibar](https://github.com/calacademy-research/minibar) | 1.0 | GitHub download & install


### Command Line Options
Argument | Required | Type | Default | Example |Help 
---------|----------|------|---------|---------|----
\-conf | True | string | config.yml | -conf config.yml |Config file (yml)
\-all | True | Flag | False | -all | Set this flag to run the complete pipeline
\-s | True | string | None | -s 1234 | Only run parts of the pipeline, i.e.:  <ul> 1. create concat fasta <br /> 2. run minibar <br /> 3. create working directories <br /> 4. run faidx <br /> 5. run nanoplot <br /> 6. run minimap 2 pipeline <br /> 7. run canu pipeline </ul>

Note:   
Either -all or -s must be set in order to run the pipeline (mutually exclusive).

Example:   
> python pipeline.py -conf config.yml -all    
> python pipeline.py -conf config.yml -s 1234

## Known Issues

## To be Implemented
- Multithread minimap/ngmlr to speed up the procedure
- Add on downstream analysis

## References

Wouter De Coster, Svenn D’Hert, Darrin T Schultz, Marc Cruts, Christine Van Broeckhoven, NanoPack: visualizing and processing long-read sequencing data, Bioinformatics, Volume 34, Issue 15, 01 August 2018, Pages 2666–2669, https://doi.org/10.1093/bioinformatics/bty149

Edge, P., Bansal, V. Longshot enables accurate variant calling in diploid genomes from single-molecule long read sequencing. Nat Commun 10, 4660 (2019). https://doi.org/10.1038/s41467-019-12493-y

Sedlazeck, F.J., Rescheneder, P., Smolka, M. et al. Accurate detection of complex structural variations using single-molecule sequencing. Nat Methods 15, 461–468 (2018). https://doi.org/10.1038/s41592-018-0001-7

Twelve years of SAMtools and BCFtools
Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, Heng Li
GigaScience, Volume 10, Issue 2, February 2021, giab008, https://doi.org/10.1093/gigascience/giab008

Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. https://doi.org/10.1093/bioinformatics/bty191

Koren S, Walenz BP, Berlin K, Miller JR, Phillippy AM. Canu: scalable and accurate long-read assembly via adaptive k-mer weighting and repeat separation. Genome Research. (2017).

