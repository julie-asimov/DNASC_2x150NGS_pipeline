#!/usr/bin/env python

#sudo docker run --memory=8g --memory-swap=8g -v /home/juliehachey/new_ngs:/code -v /home/juliehachey/analysis:/data -v $(pwd):/output -v /home/juliehachey:/meta ngs /bin/bash -c "micromamba run -n ngs  python /code/run_pipeline.py /meta/20230313_SequencingRun_060.csv /output"




import argparse
import os
import subprocess
import sys


# Import the necessary Python scripts
from get_fasta import get_fasta
from fastp import run_fastp_parallel
from qc_fastp import qc_fastp
from bwa_mapping import run_workflow
from parse_brc import run_bam
from final_metrics import write_sample_output
from final_summary import summary_file

def main():
    # Create necessary output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
     # Call the necessary Python scripts
    get_fasta(args.metadata_csv, args.master_fasta, args.output_dir)
    run_fastp_parallel(args.raw_fastq, args.output_dir)
    qc_fastp(args.output_dir)
    run_workflow( args.metadata_csv, args.output_dir)
    run_bam(args.metadata_csv, args.output_dir)
    write_sample_output(args.metadata_csv, args.output_dir)
    summary_file(args.metadata_csv, args.output_dir)

if __name__ == "__main__":
    # Define command line arguments
    parser = argparse.ArgumentParser(description='NGS Pipeline')
    parser.add_argument('raw_fastq', type=str, help='Path to input fastq')
    parser.add_argument('metadata_csv', type=str, help='Path to input samplesheet')
    parser.add_argument('master_fasta', type=str, help='Path to master fasta')
    parser.add_argument('output_dir', type=str, help='Path to output file')

    # Parse command line arguments
    args = parser.parse_args()
    # Check if input file exists
    if not os.path.isfile(args.metadata_csv):
        sys.exit('Error: Sample_sheet file not found')
    main()

