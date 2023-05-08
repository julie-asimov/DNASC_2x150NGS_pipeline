FROM --platform=linux/amd64 mambaorg/micromamba:latest as builder

USER root
RUN apt-get update && apt-get install -y build-essential cmake

# Create a new environment for the bioinformatics tools
RUN micromamba create -y -n ngs -c bioconda -c conda-forge \
    fastp=0.23.2 \
    pandas \
    gatk4=4.3.0.0 \
    bwa-mem2=2.2.1 \
    picard=3.0.0 \
    pysam=0.20.0 \
    igv-reports=1.7.0 \
    blast=2.13.0 \
    git \
    samtools=1.16.1

# Install merqury and Perl dependencies
SHELL ["micromamba", "run", "-n", "ngs", "/bin/bash", "-c"]
ENV PATH /opt/micromamba/envs/ngs/bin:$PATH

RUN ls
# Clone the bam-readcount repository from GitHub
RUN git clone https://github.com/genome/bam-readcount

RUN cd bam-readcount && mkdir build && cd build && cmake .. && make

# Install bam-readcount into the ngs environment
RUN cp bam-readcount/build/bin/bam-readcount /opt/conda/envs/ngs

# Install additional Python packages
RUN pip install \
    biopython==1.81

# Activate the bioinfo environment
SHELL ["/bin/bash", "-c"]
RUN echo "source activate bioinfo" >> ~/.bashrc
                                                         
