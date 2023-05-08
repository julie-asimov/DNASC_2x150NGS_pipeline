#!/usr/bin/env python

from subprocess import PIPE, run
import itertools
import os


#functions

def out(command):
    """run packages within python script"""
    result = run(command, stdout=PIPE, stderr=PIPE, universal_newlines=True, shell=True)
    return result.stdout

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
    # for a, b in itertools.groupby(list(enumerate(i), lambda pair: pair[1] - pair[0])):
    #     numbers=list(b)
    #     yield numbers[0][1], b[-1][1]

###bam QC metric
def bam_depth(file_location):
    """average read depth"""
    command_depth="samtools depth -a " + str(file_location) +  "| awk '{c++;s+=$3}END{print s/c}'"
    avg_depth=str(out(command_depth).strip())
    return(avg_depth)

def bam_mapped(file_location):
    """% reads mapped"""
    command_mapped="samtools flagstat " + str(file_location) +  "| awk -F "+'"'+"[(|%]"+ '"'+" 'NR== 5 {print $2}'"
    perc_mapped=str(out(command_mapped).strip())
    return(perc_mapped)


def bam_coverage(file_location):
    """breatht of coverage"""
    command_coverage="samtools depth -a " + str(file_location) +  " | awk '{c++; if($3>0) total+=1}END{print (total/c)*100}'"
    coverage=str(out(command_coverage).strip())
    return(coverage)

def bam_length(file_location):
    """plasmid length"""
    command_length="samtools depth -a " + str(file_location) +  "| awk '{c++}END{print c}'"
    length=out(command_length)
    return(length)


# Define a function to check file size and delete empty files
def check_file_size(filename):
    if os.path.getsize(filename) == 0:
        os.remove(filename)
        raise Exception(f"{filename} is empty and has been deleted.")
    else:
        print(f"{filename} has a size greater than 0 bytes.")

# Example usage
filename = "example.txt"
try:
    check_file_size(filename)
except Exception as e:
    print(e)
