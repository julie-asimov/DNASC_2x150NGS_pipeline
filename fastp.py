import argparse
import concurrent.futures
import os
import subprocess
import sys
from multiprocessing import Pool
from subprocess import run
#from scripts.functions import check_file_size

def run_fastp(f,fastq_folder, output_dir):
    filename = f.rpartition('_R1_001.fastq.gz')[0]
    f2 = filename + '_R2_001.fastq.gz'
    catch = os.path.join(output_dir,"fastp", filename + '.fastp.html')
    if not os.path.exists(catch):
        inR1 = os.path.join(fastq_folder, f)
        inR2 = os.path.join(fastq_folder, f2)
        outR1 = os.path.join(output_dir, 'fastp', f)
        outR2 = os.path.join(output_dir,'fastp', f2)
        subprocess.run([
        'fastp',
        '--trim_tail1='+ str(1),
        '--detect_adapter_for_pe',
        '--adapter_sequence=CTGTCTCTTATACACATCT',
        '--cut_right',
        '--overrepresentation_analysis',
        '--html='+ os.path.join(output_dir,"fastp", filename + '.fastp.html'),
        '--json='+ os.path.join(output_dir, "fastp", filename + '.fastp.json'),
        '-i',  inR1,
        '-I',  inR2,
        '-o',  outR1,
        '-O',  outR2
        ])



def run_fastp_parallel(raw_fastq_dir, output_dir):
    os.makedirs(os.path.join(output_dir, "fastp"), exist_ok=True)
    # Create a list of R1 file names that don't have "Undetermined" in the name
    file_list = [f for f in os.listdir(raw_fastq_dir) if 'R1' in f and 'Undetermined' not in f]

    # Create a list of arguments for starmap
    arg_list = [(f, raw_fastq_dir, output_dir) for f in file_list]

    # Run the function in parallel using a process pool
    with Pool(processes=4) as pool:
        pool.starmap(run_fastp, arg_list)

if __name__ == "__main__":
    parser = argparse.ArgumentParser("requires: (1) raw fastq folder (2) output directory")
    parser.add_argument('raw_fastq_dir')
    parser.add_argument('output_dir')
    args = parser.parse_args()

    run_fastp_parallel(args.raw_fastq_dir, args.output_dir)
    

