# Edge mapping for guiding CCIC-based recovery

## Description

These scripts weave together a simple process to index a clone library to
guide clone recovery. This script expects one or more sets of Illumina MiSeq
paired-end reads (although it is reasonable to extend to other next-gen read
technologies) for PCR amplicons of known sequence tags inserted by Tn5
transposase tagmentation. Amplicons should cover the CCIC barcodes on the R1
(forward) reads. See publication for details.

The main script is the `bc_otus.sh` BASH script. It should be taken as an
example or template as parts of it are hardcoded for the data sets used
for publication.

A helper script, `junction_mapping.pl`, is used to recover the reads which were
collapsed into OTUs.

A reference genome is also expected. In this case, a hardcoded path to the S.
albidoflavus J1074 genome used for the experiment.

## Requirements

- Read data as described above, in a single directory, fastq.gz format
- SeqKit >= 2.1.0
- seqtk >= 1.3-r106
- vsearch >= 2.18.0
- minimap2 >= 2.24-r1122

This script has only been tested on CentOS 7, but should work on most Linux
based systems with a GNU-compatible coreutils.

## Running

```bash
# Create a workspace directory and cd into it
cd workspace
# Run bc_otus script and give it the path to a directory containing all
# of your read data, in fastq.gz format.
/path/to/ccic-edge-mapping/scripts/bc_otus.sh /path/to/fastqs
# Script will automatically create a symlink to the FASTQ directory and execute
# all steps. It may not stop on some steps if an error occurs.
```

Output will be written to `workspace/minimap2/trimmed_ref_overlaps.sam`, which
can then be visualized in your preferred viewer.

