import os
import argparse
from subprocess import Popen, PIPE

import yaml



###########################################################################################
########################          Argument Parser             #############################
###########################################################################################

parser = argparse.ArgumentParser(description="Nanopore Sequencing Pipeline", formatter_class=argparse.RawTextHelpFormatter)

## ADD Arguments
parser.add_argument("-conf", default="config.yml", required=True, type=str, help="Config file (yaml)")
parser.add_argument("-ref", default="reference_map.yml", type=str, help="Map of references (yaml)")


args = parser.parse_args()


###########################################################################################
########################          Auxilary Function           #############################
###########################################################################################

def create_working_dicts(idx):
    os.mkdir(os.path.join(path, "00_analysis"))
    os.mkdir(os.path.join(path, "00_analysis/NanoPlot"))
    os.mkdir(os.path.join(path, "00_analysis/canu"))
    os.mkdir(os.path.join(path, "00_analysis/minimap2"))
    os.mkdir(os.path.join(path, "00_analysis/ngmlr"))

    for i in range(1, idx+1):
        os.mkdir(os.path.join(path, f"00_analysis/canu/barcode{'{0:03}'.format(i)}"))
        os.mkdir(os.path.join(path, f"00_analysis/minimap2/barcode{'{0:03}'.format(i)}"))
        os.mkdir(os.path.join(path, f"00_analysis/ngmlr/barcode{'{0:03}'.format(i)}"))

        
def run_faidx(refs):
    for ref in refs:
        process = Popen(["samtools", "faidx", ref], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        if stderr: print(stderr)


def run_nanoplot():
    cmd = f"NanoPlot --summary sequencing_summary.txt --loglength -o {path}/00_analysis/NanoPlot"

    process = Popen(args=cmd.split(" "),
                    stdout=PIPE, stderr=PIPE)

    stdout, stderr = process.communicate()

    if stderr: print(stderr)


def run_minimap2(ref, barcode, idx, print_out=False):
    ## Minimap2
    cmds = [
        f"minimap2 -ax map-ont {ref} {barcode} > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam",
        f"samtools view -S -b {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam",
        f"samtools sort {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam -o {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
        f"samtools index {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
        f"samtools stats {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam | grep '^SN' | cut -f 2- > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_map.txt",
        f"samtools idxstats {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}-reads-chr.txt",
        f"samtools depth {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.coverage"

    ]

    for c in cmds:
        process = Popen(args=c,
   	                stdout=PIPE, stderr=PIPE, shell=True)

        stdout, stderr = process.communicate()
        if stderr: print(stderr, "\n")
        #if print_out: print(stdout)


def run_ngmlr(ref, barcode, idx, print_out=False):
    cmds = [
        f"ngmlr -t 8 -x ont -r {ref} -q {barcode} -o {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.sam",
        f"samtools view -S -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.sam > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.bam",
        f"samtools sort {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.bam -o {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam",
        f"samtools index {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam",
        f"samtools depth {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.coverage",
        f"samtools stats {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam | grep '^SN' | cut -f 2- > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}-reads-chr.txt",
        f"sniffles -m {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr-10.bedpe -t 8",
        f"sniffles -m {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam -s 5 -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr-5.bedpe -t 8"
    ]

    for c in cmds:
        process = Popen(args=c,
                    stdout=PIPE, stderr=PIPE, shell=True)

        stdout, stderr = process.communicate()
        if stderr: print(stderr.decode("ascii"))
        if print_out: print(stdout)



def run_longshot(ref, idx, mapper, E=0.5):
    if mapper == "minimap2":
        cmd = f"longshot --bam {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam --ref {ref} --out {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_longshot.vcf -E {E}"
    elif mapper == "ngmlr":
        cmd = f"longshot --bam {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam --ref {ref} --out {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_longshot.vcf -E {E}"    
    else:
        print("ERROR: Please specify a correct mapper!")
        return None
    
    process = Popen(args=cmd.split(" "), 
                    stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = process.communicate()
    if stderr: print(stderr.decode("ascii"))

        
def run_canu(barcode, genome, idx):
    cmd = f"canu -p barcode{idx} -d {pathCanu}/barcode{idx} genomeSize={genome} -nanopore-raw {barcode} correctedErrorRate=0.14 stopOnLowCoverage=0"
    
    process = Popen(args=cmd.split(" "), 
                    stdout=PIPE, stderr=PIPE)
    
    stdout, stderr = process.communicate()
    if stderr: print(stderr)

      
###########################################################################################
########################             Configuration            #############################
###########################################################################################


with open(args.conf) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    #print(config)


path = config["path"]
ONTrun = config["ONTrun"]
pathRef = config["pathRef"]
num_items = config["sampleNumber"]

pathFastq = f"{path}{ONTrun}{config['pathFastq']}"
pathBarcode = f"{pathFastq}{config['pathBarcode']}"
pathCanu = f"{path}{config['pathCanu']}"
pathMinimap2 = f"{path}{config['pathMinimap2']}"
pathNgmlr = f"{path}{config['pathNgmlr']}"


###########################################################################################
########################             Referene Map             #############################
###########################################################################################

with open(args.ref) as f:
    refmap = yaml.load(f, Loader=yaml.FullLoader)

ref_files = list(refmap["references"].values())
references = [f"{pathRef}/{r}" for r in ref_files]

if len(references) != num_items:
    print("ERROR: Number of references and barcodes is not the same, make suer the reference_map.yaml file is ") 

for ref in references:
    if not os.path.isfile(ref): 
        print(f"ERROR: Reference file {ref} dosen't exist!")



genome = ["15k"]*num_items
barcode = [f"{path}{ONTrun}_S{'{0:03}'.format(i)}.fastq" for i in range(1, num_items+1)]


###########################################################################################
########################              Main Loop               #############################
###########################################################################################



create_working_dicts(idx=num_items)

run_faidx(refs=references)

run_nanoplot()


for i in range(num_items):
    # minimap2 LOOP!
    run_minimap2(ref=references[i], barcode=barcode[i], idx="{0:03}".format(i+1), print_out=False)
    run_longshot(ref=references[i], idx="{0:03}".format(i+1), mapper="minimap2", E=0.5)


for i in range(num_items):
    # ngmlr
    run_ngmlr(ref=references[i], barcode=barcode[i], idx="{0:03}".format(i+1), print_out=False)
    run_longshot(ref=references[i], idx="{0:03}".format(i+1), mapper="ngmlr", E=0.5)
    run_canu(barcode=barcode[i], genome=genome[i], idx="{0:03}".format(i+1))

