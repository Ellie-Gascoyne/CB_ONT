# CB_ONT
### Workflow for creating OTU tables from FastQ files for ONT seq reads
## Introduction
### ONT
The ONT flow cell contains a vast number of small holes (nanopores) embedded in a membrane. Each nanopore has a connecting sensor chip which can measure the eletric current flowing through the nanopore. As the nucleic acid molecules pass through the nanopores, the current is disrupted and creates a signal. This raw data is recorded by the MinKNOW software and processed into reads (single strand of DNA). These reads are written as POD5 files. Basecalling algorithms, such as Guppy or Dorado, can then process these reads into basecalls, i.e they turn assign a nuceotide base to the raw signal and during this process the POD5 files are converted into FASTQ files. As CB use multiplex primers, the FASTQ files must then be demultiplexed, and this can also be done using the basecalling algorithm. 

### Coding in Bash shell
## Working in bash
Bash (Bourne Again shell) is a shell that can be used on macOS. You can also run Z Shell (Zsh) which is the default for later mac software that can autocorrect for commands, and more powerful auto-completion.

To change the shell from zsh to bash, run the following command in the terminal:

```
chsh -s /bin/bash
```
This can be changed back to zsh by using the same command but change 'bash' to 'zsh'. 

Other common commands for bash include listing all files and directories for current directory 

```
ls
```
making a new directory:
```
mkdir directory_name
```
navigating directories:
```
cd /path/to/directory
```
displaying current directory
```
pwd
```
removing files:
```
rm file_name.txt
```
Finding files:
```
find /path-name "filename"
```
## Tools
The overall steps for the workflow will be to remove the primers, filter by size, quality filter, then classify. To do this we will use the following tools:
- **Cutadapt** (primer removal)
- **Seqkit** (size filtering)
- **Chopper** (quality filtering)
- **EMU/KRAKEN** (classification)
