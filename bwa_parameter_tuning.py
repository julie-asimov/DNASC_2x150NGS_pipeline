import os
import csv
import subprocess


def evaluate_aligner_params_bwa_mem2(fasta_path, output_file, num_threads=1):
    # Set the range of parameter values to test
    seed_lengths = range(17, 24)
    mismatch_penalties = range(2, 6)
    gap_open_penalties = [(5, 5), (6, 6), (7, 7)]
    results = []

    read1_path = f"{fasta_path}.read1.fq"
    read2_path = f"{fasta_path}.read2.fq"

    # Generate simulated reads using ART
    subprocess.run(["art_illumina", "-ss", "HS25", "-i", fasta_path, "-p", "-l", "150", "-f", "10", "-m", "400", "-s", "50", "-o", fasta_path])

    # Iterate over all parameter combinations
    for seed_length in seed_lengths:
        for mismatch_penalty in mismatch_penalties:
            for gap_open_penalty in gap_open_penalties:
                # Run BWA-MEM2 with the current parameter combination
                cmd = [
                    "bwa-mem2",
                    "mem",
                    "-t", str(num_threads),
                    "-k", str(seed_length),
                    "-B", str(mismatch_penalty),
                    "-O", f"{gap_open_penalty[0]},{gap_open_penalty[1]}",
                    fasta_path,
                    read1_path,
                    read2_path,
                ]
                result = subprocess.run(cmd, stdout=subprocess.PIPE)

                # Count the number of correctly aligned reads
                true_alignment_file = f"{fasta_path}.aln"
                subprocess.run(["bwa-mem2", "mem", "-x", "sr", fasta_path, read1_path, read2_path, "-R", "@RG\\tID:1\\tSM:simulated", "-o", true_alignment_file])
                true_alignments = set(open(true_alignment_file).read().splitlines())
                aligner_alignments = set(result.stdout.decode().splitlines())
                correct_alignments = len(true_alignments.intersection(aligner_alignments))

                # Calculate precision and sensitivity
                precision = correct_alignments / len(aligner_alignments)
                sensitivity = correct_alignments / len(true_alignments)

                # Save the result to the results list
                results.append({
                    "fasta_file": os.path.basename(fasta_path),
                    "seed_length": seed_length,
                    "mismatch_penalty": mismatch_penalty,
                    "gap_open_penalty": gap_open_penalty,
                    "precision": precision,
                    "sensitivity": sensitivity,
                })

    # Write the results to a CSV file
    with open(output_file, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["fasta_file", "seed_length", "mismatch_penalty", "gap_open_penalty", "precision", "sensitivity"])
        writer.writeheader()
        writer.writerows(results)

def main(input_folder, output_folder, num_threads=1):
    # Iterate over all fasta files in the input folder
   for fasta_file in os.listdir(input_folder):
        if fasta_file.endswith(".fasta") and not "." in fasta_file[len(fasta_file)-4:len(fasta_file)]:
            fasta_path = os.path.join(input_folder, fasta_file)
            output_file = os.path.join(output_folder, f"{fasta_file}.csv")
            
            # Evaluate aligner parameters for the current fasta file
            evaluate_aligner_params_bwa_mem2(fasta_path, output_file, num_threads=num_threads)

main("/Users/juliehachey/Downloads/output_230306/ref_fasta", "/Users/juliehachey/Downloads/output_230306/ref_fasta/results", num_threads=4)
