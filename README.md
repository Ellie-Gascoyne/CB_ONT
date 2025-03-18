# CB_ONT
Workflow for creating OTU tables from FastQ files for ONT seq reads
## Introduction and Notes on ONT
### ONT
The ONT flow cell contains a vast number of small holes (nanopores) embedded in a membrane. Each nanopore has a connecting sensor chip which can measure the eletric current flowing through the nanopore. As the nucleic acid molecules pass through the nanopores, the current is disrupted and creates a signal. This raw data is recorded by the MinKNOW software and processed into reads (single strand of DNA). These reads are written as POD5 files. Basecalling algorithms, such as Guppy or Dorado, can then process these reads into basecalls, i.e they turn assign a nuceotide base to the raw signal and during this process the POD5 files are converted into FASTQ files. As CB use multiplex primers, the FASTQ files must then be demultiplexed.
