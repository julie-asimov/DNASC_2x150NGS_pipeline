#!/usr/bin/env python
import os
import io
import argparse
import logging
import pandas as pd
import itertools
import re
import pysam
import subprocess
import math
import numpy as np
from typing import List, Tuple
from subprocess import PIPE, run
from Bio import Seq, SeqIO
from Bio.SeqUtils import gc_fraction

# bam QC metrics
def bam_depth(bam_file: str) -> str:
    stats = os.popen("samtools depth " + bam_file + "| awk '{if ($3 > 0) sum+=$3; n+=1} END { if (n > 0) print sum / n}'").read()
    if stats == '':
        return "0"
    else:
        return float(stats)


def bam_mapped(bam_file: str) -> str:
    # Run samtools flagstat and save the output to a file
    stats = os.popen("samtools flagstat " + str(bam_file) + "| awk -F "+'"'+"[(|%]"+ '"'+" 'NR== 7 {print $2}'").read()
    if "N/A" in stats:
        return "0"
    else:
        return str(round(float(stats), 2))



def bam_coverage(bam_file: str) -> str:
    stats = os.popen("samtools depth -a " + str(bam_file) +  " | awk '{c++; if($3>0) total+=1}END{if(c>0) print (total/c)*100; else print 0}'").read()
    if stats == '':
        return "0"
    else:
        return str(round(float(stats), 2))


def bam_length(bam_file: str) -> int:
    stats = os.popen( "samtools depth -a " + str(bam_file) +  "| awk '{c++}END{print c}'").read().strip()
    if stats == '':
        return "N/A"
    else:
        return int(stats)
   
def find_missing(lst):
    """find any missing number in lst """
    return [x for x in range(lst[0], lst[-1]+1)
                               if x not in lst]
def ranges(iterable):
    """get a range of missing numbers """
    iterable = sorted(set(iterable))
    for key, group in itertools.groupby(enumerate(iterable),
        lambda t: t[1] - t[0]):
        group = list(group)
        yield group[0][1], group[-1][1]
   
   
def blast_backbone_index(bbackbone_db: str, ref_fasta: str, output_dir: str) -> tuple:
    # Run BLAST command
    blast_cmd = ['blastn', '-db', bbackbone_db, '-query', ref_fasta, '-outfmt',
                 '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand',
                 '-perc_identity', '100']
    result = subprocess.run(blast_cmd, check=True, stdout=subprocess.PIPE)
    df = pd.read_csv(io.StringIO(result.stdout.decode()), sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'sstrand'])
    
    # Add a conditional statement to handle an empty DataFrame
    if df.empty:
        with open(os.path.join(output_dir, 'missing_backbone_info.txt'), 'a') as f:
            f.write(f"{os.path.basename(ref_fasta)} backbone not found")
        sseqid_prefixes = 'not_found'
        ranges_plasmid = [(0, 0)]
        ranges_backbone = [(0, 0)]
        return sseqid_prefixes, ranges_plasmid, ranges_backbone
    # Clean up df and get the length of backbone from fasta header name
 
    df[['sseqid_prefix', 'bb_len']] = df['sseqid'].str.split('.', expand=True)
    df['bb_len'] = df['bb_len'].astype(int)
    
    # Group by plasmid id and backbone and only take rows with plus strand
    grouped = df[(df['sstrand'] == 'plus')].groupby(['qseqid', 'sseqid'])
    
    # Filter each group looking for 1 in any row of column sstart and length of backbone in any row of column send
    filtered_grouped = grouped.filter(lambda x: (x['sstart'] == 1).any() & (x['send'] == x['bb_len']).any())
   
    # Remove groups with more than 2 rows e.g., a plasmid can only have one or two rows
    filtered_grouped = filtered_grouped.groupby(['qseqid', 'sseqid_prefix']).filter(lambda x: len(x) <= 2)
    
    # Create the 'range_plasmid' column
    filtered_grouped['range_plasmid'] = list(zip(filtered_grouped['qstart'], filtered_grouped['qend']))

    # Create the 'range_backbone' column
    filtered_grouped['range_backbone'] = list(zip(filtered_grouped['sstart'], filtered_grouped['send']))

    # Only keep part of the df and turn 'range_plasmid','range_backbone' into list of tuples per plasmid
    filtered = filtered_grouped[['qseqid', 'sseqid_prefix','range_plasmid','range_backbone']].copy()
    filtered = filtered.groupby(['qseqid', 'sseqid_prefix']).agg({'range_plasmid': list, 'range_backbone': list}).reset_index()

    # Initialize empty lists
    sseqid_prefixes = []
    ranges_plasmid = []
    ranges_backbone = []
    qseqid = []

    for idx, row in filtered.reset_index().iterrows():
        qseqid = row['qseqid']
        sseqid_prefixes = row['sseqid_prefix']
        ranges_plasmid = row['range_plasmid']
        ranges_backbone = row['range_backbone']
    # Check if the tuples are continuous
       # Check if the tuples are continuous
    if len(ranges_backbone) == 1:
        # Tuples are already continuous
        pass
    elif len(ranges_backbone) == 2:
        if ranges_backbone[1][1] + 1 == ranges_backbone[0][0]:
            # Tuples are continuous
            pass
        else:
            # Tuples are not continuous
            sseqid_prefixes = 'not_found'
            ranges_plasmid = [(0, 0)]
            ranges_backbone = [(0, 0)]

    # Check if output values are empty and replace with default values
    if not sseqid_prefixes:
        sseqid_prefixes = 'not_found'
    if not ranges_plasmid:
        ranges_plasmid = [(0, 0)]
    if not ranges_backbone:
        ranges_backbone = [(0, 0)]
    
    return sseqid_prefixes, ranges_plasmid, ranges_backbone
    
    
def get_homopolymer_and_snps(ref_fasta: str, listPosition: List[int], listSNP: List[str], backbone_range: List[int]) -> Tuple[List[str], str]:
    # Parse FASTA file
    fasta_sequences = SeqIO.parse(open(ref_fasta), 'fasta')
    sequence = str(next(fasta_sequences).seq).upper()

    # Find homopolymer regions
    homopolymer_str = [(
        [m.start()+1, m.end()+1],
        len(m.group()),
        m.group(1)
    ) for m in re.finditer(r'([ACGT])\1{0,}', sequence) if len(m.group()) >= 6]

    # Find positions of SNPs within homopolymer regions and not in the backbone range

    homopolymer_snp = []
    homopolymer_snp_list = []
    for ranges, amount, base in homopolymer_str:
        for num in listPosition:
            if ranges[0] <= int(num) <= ranges[1]:
                index = listPosition.index(num)
                homopolymer_snp.append(listSNP[index])
                listSNP.pop(index)
                listPosition.pop(index)
                homopolymer_snp_list.append(num)
             

    # Find positions of SNPs within homopolymer regions and in the backbone range
    homopolymer_snp_backbone = []
    hp_index= []
    for r in backbone_range:
        if len(r) == 1:
            range_start, range_end = r[0][0], r[0][1]
        else:
            range_start, range_end = r
        for num in homopolymer_snp_list:
            if int(range_start) <= int(num) <= int(range_end) and any(int(num) >= r_start and int(num) <= r_end for r_start, r_end in backbone_range):
                index = homopolymer_snp_list.index(num)
                homopolymer_snp_backbone.append(homopolymer_snp[index])
                hp_index.append(index)
    for pos in sorted(hp_index, reverse=True):
        homopolymer_snp.pop(pos)
            

    # Find positions of SNPs within the backbone range

    snp_backbone = []
    snp_index = []
    for r in backbone_range:
        if len(r) == 1:
            range_start, range_end = r[0][0], r[0][1]
        else:
            range_start, range_end = r
        for num in listPosition:
            if int(range_start) <= int(num) <= int(range_end) and any(int(num) >= r_start and int(num) <= r_end for r_start, r_end in backbone_range):
                index = listPosition.index(num)
                snp_backbone.append(listSNP[index])
                snp_index.append(index)
    for pos in sorted(snp_index, reverse=True):
        listSNP.pop(pos)
                
    return listSNP, snp_backbone, homopolymer_snp, homopolymer_snp_backbone

def parse_brc_file(parse_brc_output: str, ref_fasta: str, soft_clip_bam_file: str, backbone_name: str, plasmid_range: str, backbone_range: str ):
    
    listSNP=[]
    listPosition=[]
    known_bb_snp = []
    known_bb_index = []
    
    # get reference fasta name.fasta
    ref_fasta_id = os.path.basename(ref_fasta)
    #get name of reference without .fasta extension
    plasmid=(os.path.splitext(ref_fasta_id)[0])

    if os.path.getsize(parse_brc_output) > 0:
        with open(parse_brc_output,'r') as in_fh:
            next(in_fh)
            for line in in_fh:# Strip newline from end of line
                line = line.strip()
                # Fields are tab-separated, so split into a list on \t
                fields = line.split('\t')
                # Fields in file
                plasmid = fields[0]  # plasmid id
                if plasmid == 'chrom':
                    continue
                position =int(fields[1])  # Position (1-based)
                ref= fields[2].upper()  # Reference base
                base = fields[3]  # SNP/indel/del base(s)
                vaf=int(float(fields[4])*100)  # frequency of snp
                depth=int(fields[5])  # depth at one base
                count=int(fields[6])  # count of SNP at one base
                num_plus_strand=int(fields[7])  # #SNP plus strand
                num_minus_strand=float(fields[8])  # #SNP minus strand
                avg_basequality=float(fields[9])  # SNP base qaulity
                if base == ref or vaf <= 5 or 0.01 <= avg_basequality < 30 or count <= 5:
                    continue
              # Get reads with soft clips at the SNP position
                soft_clip_bam = pysam.AlignmentFile(soft_clip_bam_file, 'rb')
                query_names = set()
                for read in soft_clip_bam.fetch(plasmid, position - 1, position):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    if read.reference_start <= position and read.reference_end >= position:
                        query_names.add(read.query_name)
                position_soft_clip_reads_percent = round(len(query_names)/count*100)
                listSNP.append(ref+str(position)+base+"_"+str(vaf)+'%'+"(SC-"+str(position_soft_clip_reads_percent)+"%)")
                listPosition.append(position)

    if backbone_name == 'pUC-KanR_v1_mut1':
        # Create a dictionary mapping plasmid positions to backbone positions
        backbone_dict = {}
        for r in range(len(plasmid_range)):
            for i in range(plasmid_range[r][0], plasmid_range[r][1] + 1):
                backbone_position = i - plasmid_range[r][0] + backbone_range[r][0]
                backbone_dict[i] = backbone_position

        # Iterate over listPosition and check if the value is within the backbone range
        for i, pos in enumerate(listPosition):
            if pos in backbone_dict:
                bb_index = backbone_dict[pos]
                if bb_index in [317, 875]:
                    known_bb_index.append(i)
                    known_bb_snp.append(listSNP[i])

    if backbone_name == 'pL1_SpecR':
        # Create a dictionary mapping plasmid positions to backbone positions
        backbone_dict = {}
        for r in range(len(plasmid_range)):
            for i in range(plasmid_range[r][0], plasmid_range[r][1] + 1):
                backbone_position = i - plasmid_range[r][0] + backbone_range[r][0]
                backbone_dict[i] = backbone_position

        # Iterate over listPosition and check if the value is within the backbone range
        for i, pos in enumerate(listPosition):
            if pos in backbone_dict:
                bb_index = backbone_dict[pos]
                if bb_index in [1657]:
                    known_bb_index.append(i)
                    known_bb_snp.append(listSNP[i])
    # Remove known BB SNPs from listPosition_snp and listSNP_snp
    for pos in sorted(known_bb_index, reverse=True):
        listSNP.pop(pos)
        listPosition.pop(pos)
     

    return listSNP, listPosition,known_bb_snp, plasmid


def soft_clip_bam(input_bam_file, bam_readcount_output, backbone_range):
    listSoftClip=[]
    # Remove extension from input BAM file path
    input_file_name = os.path.splitext(input_bam_file)[0]
    df_counts = pd.read_table(bam_readcount_output, usecols=[1, 3],
                       names=["index", "count"])
                       
    # Generate output file paths with primary and supplementary reads
    primary_output_file = input_file_name + '.primary_soft_clipped_reads.bam'
    supplementary_output_file = input_file_name + '.supplementary_soft_clipped_reads.bam'
        
    # Open the input BAM file for reading
    bam_file = pysam.AlignmentFile(input_bam_file, "rb")

    # Open the output files for writing the soft-clipped reads
    primary_soft_clipped_reads_file = pysam.AlignmentFile(primary_output_file, "wb", template=bam_file)
    supplementary_soft_clipped_reads_file = pysam.AlignmentFile(supplementary_output_file, "wb", template=bam_file)

    # Iterate through each read in the BAM file
    for read in bam_file.fetch():
        # Check if the read has any soft-clipped bases
        if read.cigartuples:
            first_cigar = read.cigartuples[0]
            last_cigar = read.cigartuples[-1]

            if first_cigar[0] == 4 or last_cigar[0] == 4:
                # Determine if the read is primary or supplementary
                if read.is_supplementary:
                    # Write the supplementary soft-clipped read to the output file
                    supplementary_soft_clipped_reads_file.write(read)
                else:
                    # Write the primary soft-clipped read to the output file
                    primary_soft_clipped_reads_file.write(read)

    # Close the input and output files
    bam_file.close()
    primary_soft_clipped_reads_file.close()
    supplementary_soft_clipped_reads_file.close()

    # Index the primary output BAM file
    pysam.index(primary_output_file)

    # Open the output BAM file for reading
    soft_clipped_reads_file = pysam.AlignmentFile(primary_output_file, "rb")

    # Create a dictionary to store the count for each position
    position_dict = {}

    # Iterate through each soft-clipped read
    for read in soft_clipped_reads_file:
      # Check if the read has a soft-clipped cigar, a header, and a Phred quality score
        if read.cigartuples and read.query_name and read.query_qualities:
            # Determine the position of the soft-clip
            if read.cigartuples[0][0] == 4:
                position = read.reference_start + 1
            elif read.cigartuples[-1][0] == 4:
                position = read.reference_end
            else:
                continue

            # Extract the Phred quality scores from the read
            phred_scores = read.query_qualities

            # Calculate the average Phred quality score of the read
            phred_score = -10 * math.log10(sum([10**(-q/10) for q in phred_scores])/len(phred_scores))
            
            # Ignore reads with a Phred score less than 29
            if phred_score < 29:
                continue
          
            # Increment the count for the position in the dictionary
            if position in position_dict:
                position_dict[position] += 1
            else:
                position_dict[position] = 1

    # Close the input BAM file
    soft_clipped_reads_file.close()
    
    listSoftClip = []
    listBackboneSoftClip = []
    length = bam_length(input_bam_file)

    for position in [1, length]:
        position_dict.pop(position, None)
    # Check if the count for any position in position_dict is greater than 5% of the count in df_counts
    for position, count in position_dict.items():
        if position in df_counts["index"].tolist():
            idx = df_counts["index"].tolist().index(position)
            total_count = df_counts.iloc[idx]["count"]
            for r in backbone_range:
                if len(r) == 1:
                    range_start, range_end = r[0][0], r[0][1]
                else:
                    range_start, range_end = r
                if range_start <= int(position) <= range_end:
                    if int((count/total_count) * 100) >= 5:
                        listBackboneSoftClip.append(str(position) + "_" + str(int((count/total_count) * 100)) + "%")
                else:
                    perc_softclip = int((count/total_count) * 100)
                    if perc_softclip >= 5:
                        listSoftClip.append(str(position)+"_"+str(perc_softclip)+'%')

    # Return the dictionary of soft-clipped read counts by position
    return listSoftClip, listBackboneSoftClip, primary_output_file


def get_missing_positions(bam_readcount_output: str, bam_location: str) -> str:
    df = pd.read_table(bam_readcount_output, usecols=[0, 1, 3],
                       names=["plasmid", "index", "count"])
    
    # filter out indices where count is 7 or lower
    df = df[df['count'] > 7]
    
    # calculate rolling average of count
    window_size = 10
    df['rolling_avg'] = df['count'].rolling(window_size).mean()
    
    # filter out indices where percent change in rolling average is >= 90% (trying to detect recombination without softclips)
    percent_change = df['rolling_avg'].pct_change()
    df = df[percent_change.abs() < 0.9]
    
    # get list of indices to return
    col_one_list = df['index'].tolist()
    length = bam_length(bam_location)
    
    if not col_one_list:
        return {}
    
    if col_one_list[0] == 1 and col_one_list[-1] == int(length):
        missing_num = find_missing(col_one_list)
    elif col_one_list[0] != 1 and col_one_list[-1] == int(length):
        if int(col_one_list[0]) < 20:
            missing_num = find_missing(col_one_list)
        else:
            col_one_list.insert(0, int(0))
            missing_num = find_missing(col_one_list)
    elif col_one_list[0] == 1 and col_one_list[-1] != int(length):
        if col_one_list[-1] > (int(length)-20):
            missing_num = find_missing(col_one_list)
        else:
            length = int(length)+1
            col_one_list.append(length)
            missing_num = find_missing(col_one_list)
    elif col_one_list[0] != 1 and col_one_list[-1] != int(length):
        col_one_list.insert(0, int(0))
        length = int(length)+1
        col_one_list.append(length)
        missing_num = find_missing(col_one_list)

    if missing_num:
        ranges_list = list(ranges(missing_num))
        missing_position = str(list(ranges_list))[1:-1]
    else:
        missing_position = ""

    return missing_position

def find_low_coverage_regions(coverage_string, fasta_file, low_coverage_factor=0.2, min_stretch_length=6, gc_threshold=65):
    # Convert the string to a numpy array of integers
    coverage_array = np.array(list(map(int, coverage_string.split(','))))

    # Define a low coverage threshold as a fraction of the mean coverage
    low_coverage_threshold = low_coverage_factor * np.mean(coverage_array)

  # Identify stretches of low coverage
    low_coverage_stretches = []
    stretch_start = None
    for i, cov in enumerate(coverage_array):
        if cov < low_coverage_threshold:
            if stretch_start is None:
                stretch_start = i
        elif stretch_start is not None:
            stretch_end = i
            stretch_length = stretch_end - stretch_start
            if stretch_length >= min_stretch_length:
                low_coverage_stretches.append((stretch_start, stretch_end))
            stretch_start = None

 # If the final stretch extends to the end of the array, add it to the list
    if stretch_start is not None:
        stretch_end = len(coverage_array)
        stretch_length = stretch_end - stretch_start
        if stretch_length >= min_stretch_length:
            low_coverage_stretches.append((stretch_start, stretch_end))
    
    # Load the FASTA file
    records = list(SeqIO.parse(fasta_file, "fasta"))

    # Get the sequence information for each low coverage region
    low_coverage_regions = []
    for start, end in low_coverage_stretches:
        region_length = end - start
        region_seq = ""
        region_gc = 0
        for record in records:
            region_seq = str(record.seq[start:end])
            region_gc = gc_fraction(region_seq) *100
            print(region_gc)
            break
        is_homopolymer = any(region_seq[i] == region_seq[i+1] for i in range(len(region_seq)-1))
        is_high_gc = region_gc > gc_threshold

        if not is_homopolymer and not is_high_gc:
            # Identify this region as a recombination region
            low_coverage_regions.append((start, end, "unknown"))
        elif is_homopolymer:
            # Identify this region as a homopolymer
            low_coverage_regions.append((start, end, "homopolymer"))
        elif is_high_gc:
            # Identify this region as a high GC region
            low_coverage_regions.append((start, end, "high GC"))

    # Return the low coverage regions as a list of tuples
    return low_coverage_regions


    
def get_depth_from_bam(bam_file):
    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        # Get the depth at each position in the alignment
        depths = []
        for pileupcolumn in bam.pileup():
            depth = pileupcolumn.nsegments
            depths.append(depth)
        depths_string = ','.join(map(str, depths))
    return depths_string


def write_sample_output(metadata: str, output_dir: str) -> None:
    metadata_df = pd.read_csv(metadata)
    columns=['Not_Enough_Reads',
           'SEQ_NAME',
           'Plasmid_ID',
           'SNPs',
           'Missing_Positions',
           'Softclip_pileup',
           'Homopolymer_SNPs',
           'Coverage',
           'Depth',
           '%Mapped',
           'Low_coverage_regions',
           'Lacz_coverage',
           'Backbone_Found',
           'Backbone_SNPs',
           'Backbone_homopolymer_SNPs',
           'Backbone_Softclip_pileup',
           'Known_BB_SNPs'
          ]
    df = pd.DataFrame(columns=columns)
    for seq_name in metadata_df['SEQ_NAME']:
        if 'empty' in seq_name:
            continue
        # Get paths to relevant files
        stock_id = metadata_df.loc[metadata_df['SEQ_NAME'] == seq_name, 'STOCK_ID'].iloc[0]
        bam_location = os.path.join(output_dir, "bwa_mapping", f"{seq_name}.sorted.dedup.bam")
        bam_locationlz = os.path.join(output_dir, "bwa_mapping", f"{seq_name}_lacz.sorted.dedup.bam")
        bam_readcount_output = os.path.join(output_dir, "readcount_bam", f"{seq_name}.tsv")
        parse_brc_output = os.path.join(output_dir, "parse_brc", f"{seq_name}.txt")
        ref_fasta = os.path.join(output_dir, "ref_fasta", f"{stock_id}.fasta")
        print(bam_location)   
       # Get metrics for BAM file
        bam_depth_val = int(bam_depth(bam_location))
        if bam_depth_val < 10:
            new_row = {'Not_Enough_Reads' : 'True',
                       'SEQ_NAME': seq_name,
                       'Plasmid_ID': plasmid,
                       'SNPs': '',
                       'Missing_Positions': '',
                       'Softclip_pileup': '',
                       'Homopolymer_SNPs': '',
                       'Coverage': '',
                       'Depth': bam_depth_val,
                       '%Mapped': '',
                       'Low_coverage_regions': '',
                       'Lacz_coverage': '',
                       'Backbone_Found': '',
                       'Backbone_SNPs': '',
                       'Backbone_homopolymer_SNPs': '',
                       'Backbone_Softclip_pileup': '',
                       'Known_BB_SNPs': ''
                       }
        else:
            # Get metrics for BAM file
            bam_depth_val = bam_depth(bam_location)
            bam_mapped_val = bam_mapped(bam_location)
            bam_coverage_val = bam_coverage(bam_location)
            bam_LZcoverage_val = bam_coverage(bam_locationlz)

            #Get Backbone range
            bb_name, ranges_plasmid, ranges_backbone=blast_backbone_index('/code/mydb', ref_fasta, output_dir)
            # Get missing positions for BAM file
            missing_positions_val = get_missing_positions(bam_readcount_output, bam_location)
            #Get Softclips from BAM file
            listSoftClip_val, listBackboneSoftClip_val,primary_softclip_output = soft_clip_bam(bam_location, bam_readcount_output, ranges_plasmid)
            # Get list of SNPs for BAM file
            snps_list, positions_list, known_bb_snp, plasmid = parse_brc_file(parse_brc_output, ref_fasta,primary_softclip_output,bb_name, ranges_plasmid,ranges_backbone)
            snps_val, bb_snps, hp_snps, hp_bb_snps  = get_homopolymer_and_snps(ref_fasta, positions_list, snps_list, ranges_plasmid)
            # get low coverage aread and try to give them a type
            depths = get_depth_from_bam(bam_location)
            low_coverage = find_low_coverage_regions(depths, ref_fasta)

            new_row = {'Not_Enough_Reads' : '',
               'SEQ_NAME': seq_name,
               'Plasmid_ID' : plasmid,
               'SNPs': snps_val,
               'Missing_Positions': missing_positions_val,
               'Softclip_pileup': listSoftClip_val,
               'Homopolymer_SNPs' : hp_snps,
               'Coverage': bam_coverage_val,
               'Depth': bam_depth_val,
               '%Mapped': bam_mapped_val,
               'Low_coverage_regions' : low_coverage,
               'Lacz_coverage' : bam_LZcoverage_val,
               'Backbone_Found' : bb_name,
               'Backbone_SNPs' : bb_snps,
               'Backbone_homopolymer_SNPs' : hp_bb_snps,
               'Backbone_Softclip_pileup' : listBackboneSoftClip_val,
               'Known_BB_SNPs': known_bb_snp
              }
        print(new_row['SEQ_NAME'])
        df = pd.concat([df, pd.DataFrame([new_row], columns=columns)], ignore_index=True)
    
    df = df.astype(str).apply(lambda x: x.str.replace('[\[\]\{\}]', '', regex=True))
    df.to_csv(os.path.join(output_dir, "metrics.csv"), index=False)

    
if __name__ == "__main__":
    parser = argparse.ArgumentParser("requires: (1) seq run metadata.csv (2) output directory")

    parser.add_argument('metadata_csv')
    parser.add_argument('output_dir')
    args = parser.parse_args()

    write_sample_output(args.metadata_csv, args.output_dir)





#    bb_name, ranges_plasmid, ranges_backbone=blast_backbone_index('mydb', '/Users/juliehachey/Downloads/output_230306/ref_fasta/pAI-10351.fasta', '/Users/juliehachey/Downloads/output_230306/')
##    print(bb_name, ranges_plasmid, ranges_backbone)
#     #Get Softclips from BAM file
#    listSoftClip_val, listBackboneSoftClip_val,primary_softclip_output = soft_clip_bam('/Users/juliehachey/Downloads/output_230306/bwa_mapping/211635-pAI-10351-col2_S92.sorted.dedup.bam', '/Users/juliehachey/Downloads/output_230306/readcount_bam/211635-pAI-10351-col2_S92.tsv', ranges_plasmid)
#    # Get list of SNPs for BAM file
#    snps_list, positions_list, known_bb_snp, plasmid = parse_brc_file('/Users/juliehachey/Downloads/output_230306/parse_brc/211635-pAI-10351-col2_S92.txt', '/Users/juliehachey/Downloads/output_230306/ref_fasta/pAI-10351.fasta',primary_softclip_output,bb_name, ranges_plasmid,ranges_backbone)
#    print(snps_list, positions_list, known_bb_snp, plasmid)
##    snps_val, bb_snps, hp_snps, hp_bb_snps  = get_homopolymer_and_snps('/Users/juliehachey/Downloads/output_230306/ref_fasta/pAI-9890.fasta', positions_list, snps_list, ranges_plasmid)
##    print(snps_val, bb_snps, hp_snps, hp_bb_snps)

