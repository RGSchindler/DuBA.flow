import os
import logging
from dataclasses import dataclass, field

import pandas as pd 
import numpy as np 
import process



logger = logging.getLogger(__name__)

@dataclass
class Sample:
    reference: str
    sample: str
    
    reference_path: str = field(init=False)
    sample_fastq: str = field(init=False)
    sample_bam: str = field(init=False)
    workdir: str = field(init=False)

    # Stats for report
    success: bool = field(init=False, default=False)
    total_sequences: int = field(init=False, default=np.NaN)
    reads_mapped: int = field(init=False, default=np.NaN)
    reads_unmapped: int = field(init=False, default=np.NaN)
    error_rate: float = field(init=False, default=np.NaN)
    average_length: float = field(init=False, default=np.NaN)
    average_quality: float = field(init=False, default=np.NaN)


    def __post_init__(self):
        self.workdir = f"/home/Nanopore/input/output/{self.sample}"
        self.reference_path = f"{self.workdir}/{self.reference}"
        self.sample_fastq = f"{self.workdir}/{self.sample}.fastq"
        self.sample_bam = f"{self.workdir}/{self.sample}.bam"
        self.routine()


    def routine(self) -> None:
        if not os.path.exists(self.sample_fastq):
            logger.warning(f"{self.sample}.fastq does not exist")
            self.success = False
        elif not os.path.exists(self.reference_path):
            logger.warning(f"{self.reference} does not exist")
            self.success = False
        else:
            self.process()
            self.create_stats()
            self.set_stats()

    '''
        f"minimap2 -ax map-ont {ref} {barcode} > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam",
        f"samtools view -S -b {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam",
        f"samtools sort {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam -o {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
        f"samtools index {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
    '''


    def process(self) -> None:
        #mapping = process.run(f"/home/Nanopore/bin/micromamba run -n minimap2 minimap2 -ax map-ont {self.reference_path} {self.sample_fastq}")
        #process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools sort -o {self.sample_bam}", stdin=mapping.stdout, chained=False)
        #try:
        #    os.waitpid(os.getpgid(mapping), 0)
        #except:
        #    pass
        process.run(f"/home/Nanopore/bin/micromamba run -n minimap2 minimap2 -ax map-ont {self.reference_path} {self.sample_fastq} -o {self.workdir}/temp.sam", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools view -S -b {self.workdir}/temp.sam -o {self.workdir}/temp.bam", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools sort {self.workdir}/temp.bam -o {self.sample_bam}", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools index {self.sample_bam}", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools consensus -a -A -f pileup --show-ins no --show-del yes {self.sample_bam} -o {self.workdir}/consensus.pileup", chained=False)
        try:
            os.remove(f"{self.workdir}/temp.sam")
            os.remove(f"{self.workdir}/temp.bam")
        except FileNotFoundError:
            pass

    def create_stats(self) -> None:
        stats = process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools stats {self.sample_bam}")
        grep = process.run("grep ^SN", stdin=stats.stdout)
        process.run(f"cut -f 2- > {self.workdir}/stats.txt", stdin=grep.stdout, chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools idxstats {self.sample_bam} > {self.workdir}/idxstats.txt", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n samtools samtools depth -a {self.sample_bam} > {self.workdir}/coverage", chained=False)
        process.run(f"/home/Nanopore/bin/micromamba run -n deeptools plotCoverage -b {self.sample_bam} --plotFile coverage_plot --ignoreDuplicates", cwd=self.workdir, chained=False)
    

    def set_stats(self) -> None:
        ''' Read and parse stats from stats.txt for each sequencing run '''
        try:
            stats_df = pd.read_csv(f"{self.workdir}/stats.txt", sep="\t", header=None) 
            stats_df = stats_df.iloc[:,0:2]
            k = [s.split(":")[0] for s in stats_df.iloc[:,0]]
            stats = dict(zip(k, stats_df.iloc[:,1]))

            self.total_sequences = stats["sequences"]
            self.reads_mapped = stats["reads mapped"]
            self.reads_unmapped = stats["reads unmapped"]
            self.error_rate = stats["error rate"]
            self.average_length = stats["average length"]
            self.average_quality = stats["average quality"]

            if (self.total_sequences>100) and (self.reads_mapped>20):
                self.success = True
            else:
                self.success = False

        except (FileNotFoundError, KeyError) as e: 
            logger.warning(e)
            logger.warning(f"Failed to get stats for {self.sample}")
            self.success = False
