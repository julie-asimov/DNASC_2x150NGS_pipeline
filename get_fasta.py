#!/usr/bin/env python
import os
import argparse
import pandas as pd
from Bio import SeqIO
import subprocess
import sys

def get_fasta(metadata_csv, master_fasta, output_dir):
    print("get reference files started")
    os.makedirs(os.path.join(output_dir, "ref_fasta"), exist_ok=True)
    
   # Copy allSequences.fasta from Google Cloud Storage to the specified output directory. Requires bigquery access on VM.
    #catch = os.path.join(output_dir, "allSequences.fasta")
    #if not os.path.exists(catch):
     #   try:
      #      subprocess.run(
       #         ["gcloud", "storage", "cp", "gs://bios-sequences-prd-da03d51/allSequences.fasta", catch],
        #        check=True
         #   )
       # except subprocess.CalledProcessError as e:
        #    print("Error: Failed to download allSequences.fasta from Google Cloud Storage.")
         #   print(e.stderr)
          #  sys.exit(1)
           # return

    # Generate run.fa of only samples found in seq run from master fasta
    catch = os.path.join(output_dir, "ref_fasta", "run.fa")
    if not os.path.exists(catch):
        records=[]
        with open(os.path.join(output_dir, "ref_fasta", "run.fa"), "w") as out_file:
            for record in SeqIO.parse(master_fasta, "fasta"):
                if record.id in pd.read_csv(metadata_csv, usecols=['STOCK_ID'])['STOCK_ID'].unique():
                    records.append(record)
        if len(records) == 0:
            print("Warning: No records found for generating run.fa. Deleting the file.")
            if os.path.exists(catch):
                os.remove(catch)
            sys.exit(1)

        with open(catch, "w") as out_file:
            SeqIO.write(records, out_file, "fasta")

        if os.path.getsize(catch) == 0:
            print("Warning: run.fa is empty. Deleting the file.")
            os.remove(catch)
            sys.exit(1)

        if not os.path.exists(catch):
            print("Error: run.fa was not generated. Check the input data and metadata.")
            sys.exit(1)

    current_dir = os.getcwd()

    # Generate bwa-mem2 index files, samtools index file, and awk script to generate bed file for each fasta in the folder location $output_dir
    for f in os.listdir(os.getcwd()):
        if f.endswith(".fasta"):
            catch = os.path.splitext(f)[0] + ".bed"
            if not os.path.exists(os.path.join(output_dir,"ref_fasta", catch)):
                os.system("bwa-mem2 index {} &".format(f))
                os.system("samtools faidx {}".format(f))
                os.system("awk 'BEGIN {{FS=\"\\t\"}}; {{print $1 FS \"0\" FS $2}}' {}.fai > {}.bed".format(f, os.path.splitext(f)[0]))

    # Return to the original directory
    os.chdir(current_dir)

    print("get reference files finished")

if __name__ == "__main__":
    parser = argparse.ArgumentParser("requires: (1) seq run metadata.csv (2) output directory")

    parser.add_argument('metadata_csv')
    parser.add_argument('output_dir')
    args = parser.parse_args()

    get_fasta(args.metadata_csv, args.output_dir)

#    for filename in os.listdir(args.output_dir):
#        try:
#            file_path = os.path.join(args.output_dir, filename)
#            check_file_size(file_path)
#        except Exception as e:
#            print(f"Error checking file: {filename}")
#            print(e)

