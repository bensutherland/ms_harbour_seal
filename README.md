# Harbour Seal Project - Population Genetics
Manuscript analysis instruction guide and associated scripts for the Harbour Seal Population Genetics analysis. This pipeline is still in development stage, and comes with no guarantees.            
Primarily uses the repo from E. Normandeau of Labo Bernatchez for genotyping https://github.com/enormandeau/stacks_workflow, which uses Stacks v2.0.       


### Requirements    
`cutadapt` https://cutadapt.readthedocs.io/en/stable/index.html    
`fastqc` https://www.bioinformatics.babraham.ac.uk/projects/fastqc/   
`multiqc` https://multiqc.info/   
`bwa` http://bio-bwa.sourceforge.net/   
`samtools (v.1.9)` http://www.htslib.org/    
`stacks (v2.3e)` http://catchenlab.life.illinois.edu/stacks/     
`fineRADstructure` http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html     
`ape` (install within R)     
`stacks_workflow` https://github.com/enormandeau/stacks_workflow         


## 0. Obtain repo for analysis
Clone https://github.com/enormandeau/stacks_workflow and change directory into the repo. All scripts will be run from here.   

## 1. Preparing Data
### a. Set up 
1. Put all raw data in `02-raw` using cp or cp -l    
2. Prepare the sample info file (see template in repo sample_information.csv). Note: tab-delimited, even though name is .csv.    
3. Download reference genome (oyster_v9): https://www.ncbi.nlm.nih.gov/assembly/GCF_000297895.1/      
note: source citation is Zhang et al. 2012, Nature. https://www.ncbi.nlm.nih.gov/pubmed/22992520/       


### b. Clean data
View raw data with fastqc and multiqc:    
#### FastQC raw data reads
```
mkdir 02-raw/fastqc_raw
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw
```

#### Prepare lane info
```
./00-scripts/00_prepare_lane_info.sh
```

#### Run cutadapt in order to trim off adapters, and remove any reads less than 50 bp
```
./00-scripts/01_cutadapt.sh 12
```

#### Demultiplex with two rxn enzymes in parallel over multiple CPUs
```
./00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 nsiI mspI 14
```

#### Rename the files
```
./00-scripts/03_rename_samples.sh
```

#### FastQC de-multiplexed data
```
mkdir 04-all_samples/fastqc_demulti
fastqc 04-all_samples/*.fq.gz -o 04-all_samples/fastqc_demulti/ -t 14
multiqc -o 04-all_samples/fastqc_demulti/ 04-all_samples/fastqc_demulti
```

#### #TODO Here is where I will add some detail to remove low sequence samples ###
#### #TODO First need to determine a threshold ###

### c. Map reads against the reference genome
#### Index the genome
```
# (Only needed once) Change directory into the genome folder and index the genome 
bwa index -p GCA_004348235.1_GSC_HSeal_1.0_genomic ./GCA_004348235.1_GSC_HSeal_1.0_genomic.fna.gz
```

#### Align individual files against the genome
#### Edit variables within the following script then launch:    
```
# First update the script below to point towards the directory containing your genome
GENOMEFOLDER="/home/ben/Documents/genomes"
GENOME="GCA_004348235.1_GSC_HSeal_1.0_genomic"
# also comment out the module load commands

# Launch
./00-scripts/bwa_mem_align_reads.sh 14
```

### d. Genotype
#### Prepare and run gstacks
```
# Prepare the population map
./00-scripts/04_prepare_population_map.sh

# Edit and run gstacks
# Only update the NUM_CPU variable and run
./00-scripts/stacks2_gstacks_reference.sh

# Edit and run the ./00-scripts/stacks2_populations_reference.sh script
#  for haplotypes:
#populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
#    -t "$NUM_CPU" -p 5 -r 0.7 \
#    --renz nsiI --merge-sites \
#    --min-maf 0.01 \
#    --ordered-export --genepop \
#    --radpainter

# OR 

# Edit and run populations module for single SNP and fasta export, and hwe
#populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
#    -t "$NUM_CPU" -p 5 -r 0.7 \
#    --renz nsiI --merge-sites \
#    --min-maf 0.01 \
#    --ordered-export --genepop \
#    # --radpainter
#    --write-single-snp
#    --hwe --fasta-loc \



./00-scripts/stacks2_populations_reference.sh
```
