#!/usr/bin/env python
import os
import argparse
import pandas as pd
from Bio import SeqIO
#from functions import check_file_size

def get_fasta(metadata_csv, output_dir):
    print("get reference files started")
    os.makedirs(os.path.join(output_dir, "ref_fasta"), exist_ok=True)
    

    # Copy allSequences.fasta from Google Cloud Storage to the specified output directory
    catch = os.path.join(output_dir, "allSequences.fasta")
    if not os.path.exists(catch):
        os.system("gcloud storage cp gs://bios-sequences-prd-da03d51/allSequences.fasta {}".format(output_dir))

    # Generate run.fa of only samples found in seq run from master fasta
    catch = os.path.join(output_dir, "ref_fasta", "run.fa")
    if not os.path.exists(catch):
        with open(os.path.join(output_dir, "ref_fasta", "run.fa"), "w") as out_file:
            for record in SeqIO.parse(output_dir+'/allSequences.fasta', "fasta"):
                if record.id in pd.read_csv(metadata_csv, usecols=['STOCK_ID'])['STOCK_ID'].unique():
                    SeqIO.write(record, out_file, "fasta")

        current_dir = os.getcwd()
        # Change the current directory to the location of run.fa
        os.chdir(os.path.join(output_dir, "ref_fasta"))

    # Split the multi-pIA-fasta into individual fasta for each subject found in run.fa
        os.system("perl -ne 'if (/^>(\S+)/) { close OUT; open OUT, \">$1.fasta\" } print OUT' run.fa")

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

