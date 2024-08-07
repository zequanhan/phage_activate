
# Phage Activation Project

## Overview
This project involves genomic analysis for phage activation. It includes scripts to analyze genome sequences using tools like MEME, Needle, and BLASTP.

## Installation

### Install Python Dependencies
To install the necessary Python dependencies, run:

```sh
pip install -r requirements.txt
```

### Install External Tools
This project requires the following external tools:
1. **MEME Suite**: For motif-based sequence analysis.
2. **Needle**: For sequence alignment.
3. **BLASTP**: For protein-protein sequence alignment.

#### Installing MEME Suite
To install MEME Suite, follow these steps:

```sh
# Download and extract MEME Suite
wget http://meme-suite.org/meme-software/5.4.1/meme-5.4.1.tar.gz
tar -zxvf meme-5.4.1.tar.gz
cd meme-5.4.1

# Configure and install
./configure --prefix=/usr/local
make
sudo make install
```

#### Installing Needle
To install Needle, follow these steps:

```sh
# On Debian/Ubuntu
sudo apt-get update
sudo apt-get install emboss

# On CentOS/RHEL
sudo yum install epel-release
sudo yum install emboss
```

#### Installing BLASTP
To install BLASTP, follow these steps:

```sh
# Download BLAST+
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.11.0+-x64-linux.tar.gz
tar -zxvf ncbi-blast-2.11.0+-x64-linux.tar.gz
cd ncbi-blast-2.11.0+

# Add BLAST+ to PATH
export PATH=$PATH:$PWD/bin
```
