#!/usr/bin/env python

import argparse
import os
import sys
import pandas as pd
import subprocess
from multiprocessing import Pool
from subprocess import run


def bwa_mapping(stock_id, seq_name, index_id, output_dir):
    ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
    R1 = os.path.join(output_dir, "fastp", f"{index_id}_L001_R1_001.fastq.gz")
    R2 = os.path.join(output_dir, "fastp", f"{index_id}_L001_R2_001.fastq.gz")
    out_sam = os.path.join(output_dir, "bwa_mapping", f"{seq_name}.sam")
    out_lacz_sam = os.path.join(output_dir, "bwa_mapping", f"{seq_name}_lacz.sam")
    if not os.path.exists(out_sam):
    # Execute the BWA-MEM2 commands using subprocess.run
        run(["bwa-mem2", "mem", "-o", out_sam, ref_fasta, R1, R2],check=True)
        run(["bwa-mem2", "mem","-o", out_lacz_sam, os.path.join("/code/", "lacz.fasta", "lacz.fasta"), R1, R2], check=True)
 

def picard_workflow(stock_id, seq_name, output_dir):
    ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
    mapping_dir = os.path.join(output_dir, "bwa_mapping")
    in_sam = os.path.join(mapping_dir, f"{seq_name}.sam")
    in_sam_lacz = os.path.join(mapping_dir, f"{seq_name}_lacz.sam")
    name = seq_name
    name_lacz =f"{seq_name}_lacz"
    catch = os.path.join(mapping_dir, name_lacz+'.sorted.dedup.bam')
    if not os.path.exists(catch):
        #fixmate and compress bam
        cmd1 = f"samtools sort -n -O sam {in_sam} | samtools fixmate -m -O bam - {os.path.join(mapping_dir, name+'.fixmate.bam')}"
        cmd1_lacz = f"samtools sort -n -O sam {in_sam_lacz} | samtools fixmate -m -O bam - {os.path.join(mapping_dir, name_lacz+'.fixmate.bam')}"
        #samtools sort -n -O sam $PWD/$2/mapping/$name.sam | samtools fixmate -m -O bam - $PWD/$2/mapping/$name.fixmate.bam &&
        # sort
        cmd2 = f"samtools sort -O bam -o {os.path.join(mapping_dir, name+'.sorted.bam')} {os.path.join(mapping_dir, name+'.fixmate.bam')}"
        cmd2_lacz = f"samtools sort -O bam -o {os.path.join(mapping_dir, name_lacz+'.sorted.bam')} {os.path.join(mapping_dir, name_lacz+'.fixmate.bam')}"
        #samtools sort -O bam -o $PWD/$2/mapping/$name.sorted.bam $PWD/$2/mapping/$name.fixmate.bam &&
        # mark duplicates
        cmd3 = f"samtools markdup -r -S {os.path.join(mapping_dir, name+'.sorted.bam')} {os.path.join(mapping_dir, name+'.sorted.dedup.bam')}"
        cmd3_lacz = f"samtools markdup -r -S {os.path.join(mapping_dir, name_lacz+'.sorted.bam')} {os.path.join(mapping_dir, name_lacz+'.sorted.dedup.bam')}"
        #samtools markdup -r -S $PWD/$2/mapping/$name.sorted.bam $PWD/$2/mapping/$name.sorted.dedup.bam &&
        #create bai
        cmd4 = f"samtools index {os.path.join(mapping_dir, name+'.sorted.dedup.bam')}"
        #samtools index $PWD/$2/mapping/$name.sorted.dedup.bam
        subprocess.run(cmd1, shell=True)
        subprocess.run(cmd2, shell=True)
        subprocess.run(cmd3, shell=True)
        subprocess.run(cmd4, shell=True)
        subprocess.run(cmd1_lacz, shell=True)
        subprocess.run(cmd2_lacz, shell=True)
        subprocess.run(cmd3_lacz, shell=True)
    
#    if not os.path.exists(catch):
#        # filter for non standard read orientation and sort BAM file
#        cmd1 = f"samtools view -Sb -F 77 -F 141 {in_sam} | samtools sort -o {os.path.join(mapping_dir, name+'.filt.sort.bam')}"
#        cmd1_lacz =  f"samtools view -Sb -F 77 -F 141 {in_sam_lacz} | samtools sort -o {os.path.join(mapping_dir, name_lacz+'.filt.sort.bam')}"
#        # Mark duplicates
#        cmd2 = f"picard MarkDuplicates I={os.path.join(mapping_dir, name+'.filt.sort.bam')} O={os.path.join(mapping_dir, name+'.filt.sort.dedup.bam')} M={os.path.join(mapping_dir, name+'.filt.sort.dedup.metrics')} REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate"
#        cmd2_lacz = f"picard MarkDuplicates I={os.path.join(mapping_dir, name_lacz+'.filt.sort.bam')} O={os.path.join(mapping_dir, name_lacz+'.filt.sort.dedup.bam')} M={os.path.join(mapping_dir, name_lacz+'.filt.sort.dedup.metrics')} REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate"
#        # Fix mate information
#        cmd3 = f"picard FixMateInformation I={os.path.join(mapping_dir, name+'.filt.sort.dedup.bam')} O={os.path.join(mapping_dir, name+'.filt.sort.dedup.fixed.bam')} ADD_MATE_CIGAR=true IGNORE_MISSING_MATES=true TMP_DIR={mapping_dir}"
#        cmd3_lacz = f"picard FixMateInformation I={os.path.join(mapping_dir, name_lacz+'.filt.sort.dedup.bam')} O={os.path.join(mapping_dir, name_lacz+'.filt.sort.dedup.fixed.bam')} ADD_MATE_CIGAR=true IGNORE_MISSING_MATES=true TMP_DIR={mapping_dir}"
#        # Index bam file
#        cmd4 = f"samtools index {os.path.join(mapping_dir, name+'.filt.sort.dedup.fixed.bam')}"
#        subprocess.run(cmd1, shell=True)
#        subprocess.run(cmd2, shell=True)
#        subprocess.run(cmd3, shell=True)
#        subprocess.run(cmd4, shell=True)
#        subprocess.run(cmd1_lacz, shell=True)
#        subprocess.run(cmd2_lacz, shell=True)
#        subprocess.run(cmd3_lacz, shell=True)


def igv_report(stock_id, seq_name, output_dir):
    # Create the IGV report
    ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
    ref_bed = os.path.join(output_dir, "ref_fasta", f"{stock_id}.bed")
    mapping_dir = os.path.join(output_dir, "bwa_mapping")
    catch = os.path.join(output_dir, "igvreport", f"{seq_name}.html")
    if not os.path.exists(catch):
        try:
            run(["create_report",
                 ref_bed,
                 ref_fasta,
                 "--flanking=500",
                 "--tracks=" + os.path.join(mapping_dir, f"{seq_name}.sorted.dedup.bam"),
                 "--output=" + catch],
                check=True)
        except CalledProcessError as e:
            print(f"Error creating IGV report for {seq_name}: {e}")


def run_workflow(metadata_csv, output_dir):
    # Load CSV into a pandas DataFrame
    try:
        df = pd.read_csv(metadata_csv)
    except FileNotFoundError:
        print(f"Error: File {metadata_csv} not found.")
        sys.exit(1)

    # Create output directories
    os.makedirs(os.path.join(output_dir, "igvreport"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "bwa_mapping"), exist_ok=True)


    # Run bwa_mapping on all files in parallel
    with Pool(processes=4) as pool:
        pool.starmap(bwa_mapping, [(row['STOCK_ID'], row['SEQ_NAME'], row['INDEX_ID'], output_dir) for index, row in df.iterrows()])
    # Transform sam to bam
    with Pool(processes=4) as pool:
        pool.starmap(picard_workflow, [(row['STOCK_ID'], row['SEQ_NAME'], output_dir) for index, row in df.iterrows()])
    # Create the report for each sample
    with Pool(processes=4) as pool:
        pool.starmap(igv_report, [(row['STOCK_ID'], row['SEQ_NAME'], output_dir) for index, row in df.iterrows()])




if __name__ == "__main__":
    parser = argparse.ArgumentParser("requires:  (1) metadata_csv (2) output directory")
    parser.add_argument("metadata_csv")
    parser.add_argument("output_dir")
    args = parser.parse_args()

    # Run the workflow
    run_workflow(args.metadata_csv,args.output_dir )



#python bwa_mapping.py metadata.csv /path/to/output/directory --processes 4

#
#def rotate_fasta(infile, outfile):
#    # Open input and output files
#    with open(infile, "r") as inf, open(outfile, "w") as outf:
#        # Read in the entire sequence as a single string
#        seq = ""
#        for line in inf:
#            if not line.startswith(">"):
#                seq += line.strip()
#
#        # Split the sequence in the middle
#        mid = len(seq) // 2
#        left_seq = seq[:mid]
#        right_seq = seq[mid:]
#
#        # Rotate the sequence by swapping the two halves and joining in reverse order
#        new_seq = right_seq + left_seq
#
#        # Write the new sequence to the output file
#        outf.write(">rotated\n")
#        for i in range(0, len(new_seq), 60):
#            outf.write(new_seq[i:i+60] + "\n")
#
#rotate_fasta("input.fasta", "output.fasta")
#def picard_workflow(stock_id, seq_name, output_dir, bqsr=False, dbsnp_file=None):
#    ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
#    mapping_dir = os.path.join(output_dir, "bwa_mapping")
#    in_sam = os.path.join(mapping_dir, f"{seq_name}.sam")
#    name = seq_name
#    if bqsr:
#        catch = os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.recal.bam')
#    else:
#        catch = os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.bam')
#    if not os.path.exists(catch):
#        # Sort BAM file
#        cmd1 = f"samtools sort -O bam -o {os.path.join(mapping_dir, name+'.sort.bam')} -c {in_sam}"
#        # Mark duplicates
#        cmd2 = f"picard MarkDuplicates I={os.path.join(mapping_dir, name+'.sort.bam')} O={os.path.join(mapping_dir, name+'.sort.dedup.bam')} M={os.path.join(mapping_dir, name+'.sort.dedup.metrics')} REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate"
#        # Fix mate information
#        cmd3 = f"picard FixMateInformation I={os.path.join(mapping_dir, name+'.sort.dedup.bam')} O={os.path.join(mapping_dir, name+'.sort.dedup.fixed.bam')} ADD_MATE_CIGAR=true IGNORE_MISSING_MATES=true TMP_DIR={mapping_dir}"
#        # Filter for non-standard pair orientation
#        cmd4 = f"picard FilterSamReads I={os.path.join(mapping_dir, name+'.sort.dedup.fixed.bam')} O={os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.bam')} FILTER=not_primary_alignment,not_proper_pair,unmapped_mate,mate_unmapped,not_matching_ref_structure"
#        if bqsr:
#            # Base quality score recalibration
#            cmd5 = f"gatk BaseRecalibrator -I {os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.bam')} -R {ref_fasta} --known-sites {dbsnp_file} -O {os.path.join(mapping_dir, name+'.recal_data.table')}"
#            cmd6 = f"gatk ApplyBQSR -I {os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.bam')} -R {ref_fasta} --bqsr-recal-file {os.path.join(mapping_dir, name+'.recal_data.table')} -O {os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.recal.bam')}"
#        # Index bam file
#        if bqsr:
#            cmd7 = f"samtools index {os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.recal.bam')}"
#        else:
#            cmd7 = f"samtools index {os.path.join(mapping_dir, name+'.sort.dedup.fixed.filtered.bam')}"
#        subprocess.run(cmd1, shell=True)
#        subprocess.run(cmd2, shell=True)
#        subprocess.run(cmd3, shell=True)
#        subprocess.run(cmd4, shell=True)
#        if bqsr:
#            subprocess.run(cmd5, shell=True)
#            subprocess.run(cmd6, shell=True)
#        subprocess.run(cmd7, shell=True)
#Here are some scenarios where BQSR may not be necessary or appropriate:
#
#Low-quality data: BQSR requires high-quality data with a large number of aligned reads to generate accurate recalibration tables. If the sequencing data has low quality or low coverage, BQSR may not be useful or may even introduce more errors.
#Reference bias: BQSR assumes that the reference genome is correct and unbiased. If the reference genome has errors or biases, BQSR may propagate these errors or biases to the variant calls.
#Non-model organisms: BQSR relies on the availability of high-quality variant databases, such as dbSNP or 1000 Genomes Project, to identify common variants for recalibration. If the sequencing data is from a non-model organism without a well-annotated reference genome or variant database, BQSR may not be feasible or accurate.
#Downstream analysis goals: BQSR is primarily designed to improve variant calling accuracy for downstream analysis such as population genetics, phylogenetics, or clinical diagnosis. If the downstream analysis goals do not require high accuracy variant calling or if the analysis is focused on other aspects of the sequencing data, such as gene expression or epigenetics, BQSR may not be necessary.
