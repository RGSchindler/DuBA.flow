import os 
import logging
import pandas as pd 

import sample
import report
import preprocessing



logging.basicConfig(format="%(asctime)s | %(levelname)8s | %(message)s", level=logging.INFO)


def validate_input() -> tuple[pd.DataFrame, pd.DataFrame]:
    ''' Check if the required files and folders exist. '''
    if not os.path.exists("input/IndexCombination.tsv"):
        raise FileNotFoundError("IndexCombination.tsv does not exist")
    elif not os.path.exists("input/references.tsv"):
        raise FileNotFoundError("references.tsv does not exist")
    elif not os.path.exists("input/input.fastq.gz"):
        raise FileNotFoundError("input.fastq.gz does not exist")
    elif not os.path.exists("input/references"):
        raise FileNotFoundError("references folder does not exist")
    logging.info("Correct files where provided.")

    idx_df = pd.read_csv("input/IndexCombination.tsv", sep="\t")
    ref_df = pd.read_csv("input/references.tsv", sep="\t")

    if not idx_df.SampleID.equals(ref_df.SampleID):
        raise ValueError("SampleID in IndexCombination.tsv and references.tsv are not the same")
    logging.info("Successfully loaded IndexCombination.tsv and references.tsv")

    return idx_df, ref_df


def main() -> None:
    idx_df, ref_df = validate_input()
    logging.info("Start preprocessing...")
    preprocessing.Preprocessing(idx_df, ref_df)

    report_data = []
    for sampleid, ref in zip(ref_df.SampleID.values, ref_df.ReferenceFile.values):
        logging.info(f"Start processing {sampleid}")
        smpl = sample.Sample(reference=ref, sample=sampleid)
        logging.info(f"Finished processing {sampleid}")
        report_data.append([
            sampleid,
            smpl.success,
            smpl.total_sequences,
            smpl.reads_mapped,
            smpl.error_rate,
            smpl.average_length,
            smpl.average_quality,
        ])

        if smpl.success:
            logging.info(f"Generating report for {sampleid}")
            rep = report.Report(
                reference=ref,
                sample_id= sampleid
            )
            del rep

        else:
            logging.info(f"Skipping report for {sampleid}")
        del smpl

    # Generate overall report tsv    
    report_df = pd.DataFrame(report_data)
    report_df.columns = ["SampleID", "Success", "Total Sequences", "Reads Mapped", "Error Rate", "Average Length", "Average Quality"]
    report_df.to_csv("input/output/report.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()