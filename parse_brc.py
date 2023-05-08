#!/usr/bin/env python

import os
import argparse
import subprocess
import pandas as pd


def run_bam_readcount(stock_id, seq_name, output_dir):
    ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
    in_bam = os.path.join(output_dir, "bwa_mapping",f"{seq_name}.sorted.dedup.bam" )
    catch = os.path.join(output_dir, "readcount_bam", f"{seq_name}.tsv")
    if not os.path.exists(catch):
        try:
            bam_readcount_cmd = [
                "bam-readcount", "-w", "0",
                "-f", ref_fasta,
                in_bam
            ]
            with open(catch, "w") as output_file:
                subprocess.run(bam_readcount_cmd, stdout=output_file)
        except subprocess.CalledProcessError:
            print(f"Error: failed to run bam-readcount for {seq_name}")
            return


def parse_brc(stock_id, seq_name, output_dir):
    bam_readcount_output = os.path.join(output_dir, "readcount_bam", f"{seq_name}.tsv")
    parse_brc_output = os.path.join(output_dir, "parse_brc", f"{seq_name}.txt")
    if not os.path.exists(parse_brc_output):
        try:
            headers = [
                'chrom', 'position', 'ref', 'base', 'vaf', 'depth', 'count', 'num_plus_strand',
                'num_minus_strand', 'avg_basequality'
            ]

            with open(parse_brc_output, 'a') as out_fh:  # Use 'a' mode to append to the file
                out_fh.write('\t'.join(headers) + '\n')

                # IMPORTANT: this relies on Python 3.6+ to maintain insertion order
                # Each field is a key with value a function to convert to the
                # appropriate data type
                base_fields = {
                    'base': str,
                    'count': int,
                    'avg_mapping_quality': float,
                    'avg_basequality': float,
                    'avg_se_mapping_quality': float,
                    'num_plus_strand': int,
                    'num_minus_strand': int,
                    'avg_pos_as_fraction': float,
                    'avg_num_mismatches_as_fraction': float,
                    'avg_sum_mismatch_qualities': float,
                    'num_q2_containing_reads': int,
                    'avg_distance_to_q2_start_in_q2_reads': float,
                    'avg_clipped_length': float,
                    'avg_distance_to_effective_3p_end': float
                }

                if os.path.getsize(bam_readcount_output) > 0:
                    with open(bam_readcount_output) as in_fh:
                        for line in in_fh:
                            line = line.strip()
                            fields = line.split('\t')
                            chrom = fields[0]
                            position = int(fields[1])
                            reference_base = fields[2].upper()
                            depth = int(fields[3])
                            for base_data_string in fields[4:]:
                                base_data = {}
                                base_values = base_data_string.split(':')
                                for i, base_field in enumerate(base_fields.keys()):
                                    base_data[base_field] = base_fields[base_field](base_values[i])
                                if base_data['count'] == 0:
                                    continue
                                vaf = base_data['count'] / depth
                                out_fh.write('\t'.join([    str(x) for x in (chrom, position, reference_base, base_data['base'],
                                                '%0.2f' % (vaf), depth, base_data['count'],
                                                base_data['num_plus_strand'],base_data['num_minus_strand'],
                                                base_data['avg_basequality'],base_data['avg_pos_as_fraction'])
                                            ]) + '\n')
        except subprocess.CalledProcessError:
            print(f"Error: failed to run parse_brc for {seq_name}")

def run_bam(metadata_csv,output_dir):
    # Load CSV into a pandas DataFrame
    try:
        df = pd.read_csv(metadata_csv)
    except FileNotFoundError:
        print(f"Error: File {metadata_csv} not found.")
        sys.exit(1)

    # Create output directories
    os.makedirs(os.path.join(output_dir, "parse_brc"), exist_ok=True)
    os.makedirs(os.path.join(output_dir, "readcount_bam"), exist_ok=True)

    # Run bamreadcount
    for index, row in df.iterrows():
        run_bam_readcount(row['STOCK_ID'], row['SEQ_NAME'],output_dir)
    # Transform sam to bam
    for index, row in df.iterrows():
        parse_brc(row['STOCK_ID'], row['SEQ_NAME'], output_dir)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Parse bam-readcount output")
    parser.add_argument('metadata_csv', help="Metadata CSV file")
    parser.add_argument('output_dir', help="Output directory")
    args = parser.parse_args()

    run_bam(args.metadata_csv,args.output_dir)
