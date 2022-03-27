import os
import re
import argparse
import shutil
from subprocess import Popen, PIPE

import yaml


###########################################################################################
########################          Argument Parser             #############################
###########################################################################################

parser = argparse.ArgumentParser(
    description="Nanopore Sequencing Pipeline", formatter_class=argparse.RawTextHelpFormatter
)

## ADD Arguments
parser.add_argument("-conf", default="config.yml", required=True, type=str, help="Config file (yaml)")
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument("-all", action="store_true", help="Set flag to run the whole pipeline")
group.add_argument("-s", type=int, help="Run seperate parts of the pipeline")

args = parser.parse_args()


###########################################################################################
########################          Auxilary Function           #############################
###########################################################################################


def create_concat_fasta(non_demultiplexed, concat_path):
    files = os.listdir(non_demultiplexed)

    with open(concat_path, "wb") as fd:
        for file in files:
            with open(file, "rb") as fs:
                shutil.copyfileobj(fs, fd)


def run_minibar(path, barcode_file, concat_path):

    cmd = f'python minibar.py -T -F {barcode_file} {concat_path} -P ""'

    process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)

    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode("ascii"))

    for file in os.listdir():
        if file.split(".")[-1] == "fastq":
            if file.count("_") == 1:
                shutil.move(file, f"{path}/{ONTrun}{file}")
            else:
                shutil.move(
                    file,
                    f"{path}/misc/{ONTrun}{file}",
                )


def mkdir(folder_path):
    # Check if the folders already exists
    if not os.path.isdir(folder_path):
        os.mkdir(folder_path)


def create_working_dirs(idx):
    mkdir(os.path.join(path, "00_analysis"))
    mkdir(os.path.join(path, "00_analysis/NanoPlot"))
    mkdir(os.path.join(path, "00_analysis/canu"))
    mkdir(os.path.join(path, "00_analysis/minimap2"))
    mkdir(os.path.join(path, "00_analysis/ngmlr"))

    for i in idx:
        mkdir(os.path.join(path, f"00_analysis/canu/barcode{i}"))
        mkdir(os.path.join(path, f"00_analysis/minimap2/barcode{i}"))
        mkdir(os.path.join(path, f"00_analysis/ngmlr/barcode{i}"))


# Only run on unique
def run_faidx(run_ids):
    for rid in run_ids:
        # print(references[rid])
        process = Popen(["samtools", "faidx", references[rid]], stdout=PIPE, stderr=PIPE)
        stdout, stderr = process.communicate()
        print(stdout.decode("ascii"))
        if stderr:
            print(stderr.decode("ascii"))


def run_nanoplot():
    cmd = f"NanoPlot --summary sequencing_summary.txt --loglength -o {path}/00_analysis/NanoPlot"

    process = Popen(args=cmd.split(" "), stdout=PIPE, stderr=PIPE)

    stdout, stderr = process.communicate()

    if stderr:
        print(stderr.decode("ascii"))


def run_minimap2(ref, barcode, idx, print_out=False):
    ## Minimap2
    cmds = [
        f"minimap2 -ax map-ont {ref} {barcode} > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam",
        f"samtools view -S -b {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.sam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam",
        f"samtools sort {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.bam -o {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
        f"samtools index {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam",
        f"samtools stats {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam | grep '^SN' | cut -f 2- > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_map.txt",
        f"samtools idxstats {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}-reads-chr.txt",
        f"samtools depth {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam > {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}.coverage",
    ]

    for c in cmds:
        # print(c, "\n")
        process = Popen(args=c, stdout=PIPE, stderr=PIPE, shell=True)

        stdout, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
        if print_out:
            print(stdout)


def run_ngmlr(ref, barcode, idx, print_out=False):
    cmds = [
        f"ngmlr -t 8 -x ont -r {ref} -q {barcode} -o {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.sam",
        f"samtools view -S -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.sam > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.bam",
        f"samtools sort {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.bam -o {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam",
        f"samtools index {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam",
        f"samtools depth {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr.coverage",
        f"samtools stats {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam | grep '^SN' | cut -f 2- > {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}-reads-chr.txt",
        f"sniffles -m {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr-10.bedpe -t 8",
        f"sniffles -m {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam -s 5 -b {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr-5.bedpe -t 8",
    ]

    for c in cmds:
        process = Popen(args=c, stdout=PIPE, stderr=PIPE, shell=True)

        stdout, stderr = process.communicate()
        if stderr:
            print(stderr.decode("ascii"))
        if print_out:
            print(stdout)


def run_longshot(ref, idx, mapper, E=0.5):
    if mapper == "minimap2":
        cmd = f"longshot --bam {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_sort.bam --ref {ref} --out {pathMinimap2}/barcode{idx}/{ONTrun}S{idx}_longshot.vcf"
    elif mapper == "ngmlr":
        cmd = f"longshot --bam {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_ngmlr_sort.bam --ref {ref} --out {pathNgmlr}/barcode{idx}/{ONTrun}S{idx}_longshot.vcf"
    else:
        raise ("ERROR: Please specify a correct mapper!")

    # print(cmd, "\n")
    process = Popen(args=cmd.split(" "), stdout=PIPE, stderr=PIPE)

    stdout, stderr = process.communicate()
    if stderr:
        print(stderr.decode("ascii"))


def run_canu(barcode, genome, idx):
    cmd = f"canu -p barcode{idx} -d {pathCanu}/barcode{idx} genomeSize={genome} -nanopore-raw {barcode} correctedErrorRate=0.14 stopOnLowCoverage=0"

    print(cmd)
    process = Popen(args=cmd, stdout=PIPE, stderr=PIPE, shell=True)

    stdout, stderr = process.communicate()

    print(stdout)
    if stderr:
        print(stderr.decode("ascii"))

    """
    print(stdout.decode("ascii"))
    if stderr: print(stderr.decode("ascii"))
    """


###########################################################################################
########################             Configuration            #############################
###########################################################################################


with open(args.conf) as f:
    config = yaml.load(f, Loader=yaml.FullLoader)
    # print(config)


path = config["path"]
ONTrun = config["ONTrun"]
pathRef = config["pathRef"]
references_map = config["refMap"]
non_demultiplexed = config["nonDemultiplexed"]
bc_file = config["barcodeFile"]
concat_fasta_path = config["concat_fasta_path"]

pathFastq = f"{path}{ONTrun}{config['pathFastq']}"
pathBarcode = f"{pathFastq}{config['pathBarcode']}"
pathCanu = f"{path}{config['pathCanu']}"
pathMinimap2 = f"{path}{config['pathMinimap2']}"
pathNgmlr = f"{path}{config['pathNgmlr']}"


###########################################################################################
########################             Referene Map             #############################
###########################################################################################

with open(references_map) as f:
    refmap = yaml.load(f, Loader=yaml.FullLoader)

ref_files = list(refmap["references"].values())
references = [f"{pathRef}/{r}" for r in ref_files]

"""
if len(references) != num_items:
    print("ERROR: Number of references and barcodes is not the same, make sure the reference_map.yaml file is ") 
"""

for ref in references:
    if not os.path.isfile(ref):
        raise(f"ERROR: Reference file {ref} dosen't exist!")


# Check demultiplex output:
regex = r"(?!.*\_.*)(?<=S).*(?=\.)"  # get the running number from file_name
run_ids = []
for file_name in os.listdir(path):
    id_buffer = re.findall(regex, file_name, re.MULTILINE)
    if id_buffer:
        run_ids.extend(id_buffer)


genome = ["15k"] * len(run_ids)
barcode = [f"{path}{ONTrun}S{i}.fastq" for i in run_ids]

barcode = {}
references = {}
for rid in run_ids:
    r_buffer = refmap["references"][f"ref{rid}"]

    barcode[rid] = f"{path}{ONTrun}S{rid}.fastq"
    references[rid] = f"{pathRef}/{r_buffer}"

# barcode = sorted(buffer_barcode.items())
# references = sorted(buffer_references.items())


###########################################################################################
########################              Main Loop               #############################
###########################################################################################

def minimap_pipe():
    for rid in run_ids:
        run_minimap2(ref=references[rid], barcode=barcode[rid], idx=rid, print_out=False)
        run_longshot(ref=references[rid], idx=rid, mapper="minimap2", E=0.5)


def canu_pipe():
    for rid in run_ids:
        run_ngmlr(ref=references[rid], barcode=barcode[rid], idx=rid, print_out=False)
        run_longshot(ref=references[rid], idx=rid, mapper="ngmlr", E=0.5)
        run_canu(barcode=barcode[rid], genome="15k", idx=rid)


if concat_fasta_path:
    concat_path = concat_fasta_path
else:
    concat_path = f"{non_demultiplexed}{ONTrun}_ConcatFasta.fastq.gz"


pipeline = {
    "1": create_concat_fasta(non_demultiplexed=non_demultiplexed, concat_path=concat_path),
    "2": run_minibar(path=path, barcode_file=bc_file, concat_path=concat_path),
    "3": create_working_dirs(idx=run_ids),
    "4": run_faidx(run_ids),
    "5": run_nanoplot(),
    "6": minimap_pipe(),
    "7": canu_pipe(),
}

if args.all:
    for i in range(1, 8):
        pipeline[str(i)]

elif args.s:
    sep_pipe = [s for s in args.s]
    for sep in sep_pipe:
        pipeline[sep]

else:
    raise ("ERROR, neither --all nor --s was set")