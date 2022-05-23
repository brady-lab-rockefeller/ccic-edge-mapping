[![DOI](https://zenodo.org/badge/495478595.svg)](https://zenodo.org/badge/latestdoi/495478595)

# Edge mapping for guiding CCIC-based retrieval

## Description

A simple process to index the barcodes of a CCIC clone library to guide
retrieval of target sequences. This script expects Illumina MiSeq paired-end
reads (although it is reasonable to extend to other next-gen read technologies)
for PCR amplicons containing the vector barcode and the edge of a captured
sequence. Amplicons should cover the CCIC barcodes on the R1 (forward) reads.
See publication for details.

*NB:* parts of these scripts may use hardcoded paths or assumptions relevant
only to the data sets used for publication.

The main script is the `bc_otus.sh` BASH script. 

A helper script, `junction_mapping.pl`, is used to recover the reads which were
collapsed into OTUs. It is called by `bc_otus.sh`.

A reference genome is also expected. In this case, a hardcoded path to the S.
albidoflavus J1074 genome is used for the experiment.

## Requirements

- Read data as described above, in a single directory, fastq.gz format SeqKit >=
- 2.1.0 seqtk >= 1.3-r106 vsearch >= 2.18.0 minimap2 >= 2.24-r1122

This script has only been tested on CentOS 7, but should work on most Linux
based systems with a GNU-compatible coreutils.

A conda environment has been provided to ensure all required software is
available. You can use it in conda with:

```bash conda env create -f ccic-edge-mapping_env.yml -n ccic-edge-mapping # Or
with mamba (much faster) mamba env create -f ccic-edge-mapping_env.yml -n
ccic-edge-mapping ```

## Running

```bash # Create a workspace directory and cd into it cd workspace # Run bc_otus
script and give it the path to a directory containing all # of your read data,
in fastq.gz format. /path/to/ccic-edge-mapping/scripts/bc_otus.sh
/path/to/fastqs # Script will automatically create a symlink to the FASTQ
directory and execute # all steps. It may not stop on some steps if an error
occurs. ```

Output will be written to `workspace/minimap2/trimmed_ref_overlaps.sam`, which
can then be visualized in your preferred viewer.

