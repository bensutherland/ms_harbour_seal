#!/bin/bash
# Normalize all fastq files to a set number by random sampling
# Requires seqtk, run from a stacks_workflow repository


# Global variables
INPUT_FOLDER="04-all_samples"
OUTPUT_FOLDER="04-all_samples/normalized"

# User set variables
READ_COUNT=1250000

# Create directory for normalized reads
mkdir $OUTPUT_FOLDER 2>/dev/null

# Normalize reads in the folder
ls -1 $INPUT_FOLDER/*.fq.gz |
    sort -u |
    while read i
    do
        echo "Normalizing $i to $READ_COUNT reads"

        # Retain only the filename, not the path
        name=$(basename $i)

        # Use seqtk to sample the specified reads randomly from the input file
        seqtk sample $INPUT_FOLDER/$name $READ_COUNT |
            
            # Compress the output and save to the output folder
            gzip - > $OUTPUT_FOLDER/$name

    done

# Reporting
echo "Your reads are now in $OUTPUT_FOLDER. You can make a new directory for your non-normalized reads and put the normalized ones in the main folder to continue with the pipeline with normalized input." 

