# Harbour Seal Project - Population Genetics
Manuscript analysis for the harbour seal population genetics analysis. This pipeline is still in development stage, and comes with no guarantees.            
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
3. Download GenBank version reference genome: https://www.ncbi.nlm.nih.gov/genome/?term=Phoca+vitulina      


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

## 2. Analyze data
### a. Map reads against the reference genome
#### Index the genome
```
# (Only needed once) Change directory into the genome folder and index the genome 
bwa index -p GCA_004348235.1_GSC_HSeal_1.0_genomic ./GCA_004348235.1_GSC_HSeal_1.0_genomic.fna.gz
```

#### Align individual files against the genome
Edit variables within the following script then launch:    
```
# First update the script below to point towards the directory containing your genome
GENOMEFOLDER="/home/ben/Documents/genomes"
GENOME="GCA_004348235.1_GSC_HSeal_1.0_genomic"
# also comment out the module load commands

# Launch
./00-scripts/bwa_mem_align_reads.sh 14
```

### b. Inspect alignment results
Compare per-sample reads and alignments, and per-sample reads and number of aligned scaffolds:      
```
./../ms_harbour_seal/01_scripts/assess_results.sh    # also uses plotting Rscript 'assess_reads_mappings.R'   
./../ms_harbour_seal/01_scripts/determine_number_unique_scaff_mapped.sh
# These scripts will produce: 
# 04-all_samples/reads_per_sample_table.txt 
# 04-all_samples/mappings_per_sample.txt
# graph: number reads and aligned reads
# graph: number reads and scaffolds mapped

```

Other optional calculations:    
```
# Total reads in all samples:     
awk '{ print $2 } ' 04-all_samples/reads_per_sample_table.txt | paste -sd+ - | bc
# Total reads before de-multiplexing (note: divide by 4 due to fastqc):   
for i in $(ls 02-raw/*.fastq.gz) ; do echo $i ; gunzip -c $i | wc -l ; done
```

#### Remove low sequence depth samples 
In general, this requires running the below first to see if there are any odd effects of low coverage, although one could also look to the mean coverage above and remove any clear outliers.     
For this project, we will remove any individuals with fewer than 1 M reads. There are two that have almost no reads, so these are an easy choice to remove. A few have between 1 M - 1.5 M, we will keep an eye on those.    
```
# make a directory to store removed samples
mkdir 04-all_samples/removed_samples

# move problematic samples to the directory
mv 04-all_samples/NBC_110* 04-all_samples/removed_samples/
mv 04-all_samples/ORE_104.* 04-all_samples/removed_samples/
```  
Now you must also correct `01-info_files/sample_information.csv` by dropping those two samples as well, before building the population map below.       
```
# First create an all-sample backup
cp 01-info_files/sample_information.csv 01-info_files/sample_information_all_samples.csv  

# Then remove the offending samples (-P allows use of \t to identify tabs), although you cannot write to the same file, so a move is required
grep -vP 'NBC\t110|ORE\t104' 01-info_files/sample_information.csv > 01-info_files/sample_information_outlier_rem.csv
mv 01-info_files/sample_information_outlier_rem.csv 01-info_files/sample_information.csv 

# Check to make sure the correct number of samples were removed
wc -l 01-info_files/sample_information*

# now you have a backup file, and a trimmed sample_information file. Proceed below to rebuild the pop map
```

Redo the genotyping below to make sure you get all of the proper outputs of stacks with the problematic samples removed.     

### d. Genotype
#### Prepare and run gstacks
```
# Prepare the population map
./00-scripts/04_prepare_population_map.sh

# Edit and run gstacks
# Only update the NUM_CPU variable and run
./00-scripts/stacks2_gstacks_reference.sh

# Edit and run the ./00-scripts/stacks2_populations_reference.sh script
#  for single SNP and fasta export (and per locus hwe) 
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 7 -r 0.7 \
    --renz nsiI --merge-sites \
    --min-maf 0.01 \
    --ordered-export --genepop \
    --write-single-snp --hwe --fasta-loc

#  for haplotypes:
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 7 -r 0.7 \
    --renz nsiI --merge-sites \
    --min-maf 0.01 \
    --ordered-export --genepop \
    --radpainter

# Suggest to make a different folder containing everything populations.* after each of those runs in order to preserve the data from multiple populations runs.     
e.g., populations_out_single_snp and populations_out_microhaplotypes     

```

## 3. Analysis of results
Clone simple_pop_stats.       

Copy in the harbour seal colour file from an earlier version of simple_pop_stats. The file is entitled 'harbour_seal_pops_colours.csv'       

Copy out the single-variant-per-locus genepop to simple_pop_stats:      
`cp 05-stacks/populations_out_single_snp/populations.snps.genepop ../simple_pop_stats/02_input_data/`

Rename your file
`bhs_p<X>_r<0.X>_maf<0.0X>_20<XX-XX-XX>.gen` (note: customize as per true values for variables)         

Analyze via `ms_harbour_seal/01_scripts/hs_popn_analysis.R`     

## 4. Analysis of results - relatedness
##### This section not run yet #####
Note: must run differentiation analysis first.    
#### i. Shared coancestry by related
Input: single SNP, HWE filtered from adegenet output (adegenet_output.RData).      
Open Rscript `01_scripts/relatedeness.R` to translate the genlight obj to related format, and calculate coancestry per population, plotting per-population relatedness metrics.       
Output: See `09-diversity_stats/relatedness_*.pdf`   
##### /end/ This section not run yet #####

#### ii. Shared haplotypes by fineRADstructure
Input: haplotype output from stacks populations as input to fineRADstructure         
Run all from directory `ms_harbour_seal`.     

Copy out radpainter input (SimpleMatrix) generated by stacks population module (above)         
`cp ../stacks_workflow/05-stacks/populations_out_microhaplotypes/populations.haps.radpainter ./04_relatedness/populations.haps.radpainter`      

In general, follow instructions from the fineRADstructure tutorial (http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html), but in brief:  
1) Calculate co-ancestry matrix:     
`RADpainter paint 04_relatedness/populations.haps.radpainter`      
2) Assign individuals to populations:     
`finestructure -x 100000 -y 100000 04_relatedness/populations.haps_chunks.out 04_relatedness/populations.haps_chunks.out.mcmc.xml`
Uses Markov chain Monte Carlo (mcmc) without a tree, assumes data is from one pop (-I 1).    
3) Build tree:        
`finestructure -m T -x 10000 04_relatedness/populations.haps_chunks.out 04_relatedness/populations.haps_chunks.out.mcmc.xml 04_relatedness/populations.haps_chunks.out.mcmcTree.xml`      

Then plot using the Rscripts adapted from the fineRADstructure site (see above)   
`01_scripts/fineRADstructurePlot.R` (follow instructions here, running this will source the libraries below)      
`01_scripts/FinestructureLibrary.R`     
This will produce plots in the working directory.  


## 5. Re-analysis of previous data
Requires a fresh `stacks_workflow` repo.       
Put fastq files in `02-raw`, as they have already been demultiplexed, but still need to be trimmed (contain adapters).        
Compress files:     
`gzip *.fastq`

Run fastqc on the files:       
`mkdir 02-raw/fastqc/`       
`fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc/ -t 2`       
`multiqc -o 02-raw/fastqc/ 02-raw/fastqc`       

Remove adapters and trim:      
`00-scripts/01_cutadapt.sh 3`        

Run fastqc on the files:       
`mkdir 02-raw/trimmed/fastqc`      
`fastqc 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc/ -t 2`
`multiqc -o 02-raw/trimmed/fastqc 02-raw/trimmed/fastqc`      

Move the trimmed files to `04-all_samples`, which skips script `00-scripts/03_rename_samples.sh`      
`mv 02-raw/trimmed/*.fastq.gz 04-all_samples/`      

After editing the script, run the alignment       
`./00-scripts/bwa_mem_align_reads.sh 3`       

Aligned data will be in `04-all_samples`.       

Prepare a renaming script (this is required because the renaming script of stacks workflow was not used, nor was demultiplexing)         
```
grep -vE '^#' 01-info_files/sample_information.csv | awk '{ print "mv " $1".sorted.bam " $3"_"$4".sorted.bam" }' - > 00-scripts/rename_aligned_bams.sh 

# Manually add the shebang to the file (#!/bin/bash), then make it executable using chmod.      

# Run the renaming
cd 04-all_samples
./../00-scripts/rename_aligned_bams.sh    
cd ..
# Files are ready, and should all be in the format of the population map file
```

#### Genotyping
Prepare and run gstacks        
```
# Prepare the population map
./00-scripts/04_prepare_population_map.sh

# Edit and run gstacks
# Only update the NUM_CPU variable and run
./00-scripts/stacks2_gstacks_reference.sh

# Since our main goal here is to observe diversity, keep multiple SNPs per locus to capture full variation, therefore use the following edits in the script, then run it:       

#populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
#    -t "$NUM_CPU" -p 6 -r 0.7 \
#    --min-maf 0.01 \
#    --ordered-export --genepop 

./00-scripts/stacks2_populations_reference.sh

```

### de novo approach
Copy all raw data to 02-raw.     
Copy sample information to 01-info_files
Create script for renaming raw data (at fastq.gz rather than bam):      
`grep -vE '^#' 01-info_files/sample_information.csv | awk '{ print "mv " $1 " " $3"_"$4".fq.gz" }' > 00-scripts/rename_raw_fastq_for_pop_map.sh`        
...then add the shebang and chmod to make executable.       

```
cd 02-raw
./../00-scripts/rename_raw_fastq_for_pop_map.sh
cd ..

```

Run fastqc on the files:       
`mkdir 02-raw/fastqc/`       
`fastqc 02-raw/*.fq.gz -o 02-raw/fastqc/ -t 2`      
`multiqc -o 02-raw/fastqc/ 02-raw/fastqc`       

Remove adapters and trim:      
(note: first edit script to identify .fq.gz instead of .fastq.gz)      
`00-scripts/01_cutadapt.sh 3`        

Run fastqc on the files:       
`mkdir 02-raw/trimmed/fastqc`      
`fastqc 02-raw/trimmed/*.fq.gz -o 02-raw/trimmed/fastqc/ -t 2`
`multiqc -o 02-raw/trimmed/fastqc 02-raw/trimmed/fastqc`      

If running de novo, need to truncate to a constant read length. Use the following custom script to truncate all reads to 70 bp and only keep the 70 bp reads:     
```
mkdir 02-raw/standardized
./00-scripts/standardize_reads.sh
```
The output will be in the new standardized folder. 

fastqc again:      
```
mkdir 02-raw/standardized/fastqc 
fastqc 02-raw/standardized/*.fq.gz -o 02-raw/standardized/fastqc/ -t 2
multiqc -o 02-raw/standardized/fastqc 02-raw/standardized/fastqc  
```


Clone three copies of stacks workflow, for each denovo analysis. Copy in the sample info file:    
```
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_wc/01-info_files/
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_ec/01-info_files/
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_eu/01-info_files/
```

rsync in the raw data into the 04-all_data folder:     
```
rsync -avzP 02-raw/standardized/END* ./../stacks_workflow_denovo_wc/04-all_samples/
rsync -avzP 02-raw/standardized/KOD* ./../stacks_workflow_denovo_wc/04-all_samples/
rsync -avzP 02-raw/standardized/BIC* ./../stacks_workflow_denovo_ec/04-all_samples/
rsync -avzP 02-raw/standardized/NFL* ./../stacks_workflow_denovo_ec/04-all_samples/
rsync -avzP 02-raw/standardized/ORK* ./../stacks_workflow_denovo_eu/04-all_samples/
rsync -avzP 02-raw/standardized/WNL* ./../stacks_workflow_denovo_eu/04-all_samples/
```

For each:     
Make population map
`./00-scripts/04_prepare_population_map.sh`       

```
# note: for all of the following, update the NUM_CPU variable before launching
./00-scripts/stacks2_ustacks.sh
./00-scripts/stacks2_cstacks.sh
./00-scripts/stacks2_sstacks.sh
./00-scripts/stacks2_tsv2bam.sh
./00-scripts/stacks2_gstacks.sh # note: this should be gstacks_denovo
./00-scripts/stacks2_populations.sh # note: this should be populations_denovo

```

## Supplemental Information 
Note: if you are struggling with multiqc, use the following approach:      
```
# install miniconda by package installer
conda create -y -n my-conda-env
source activate my-conda-env
conda install -c bioconda -c conda-forge multiqc
# then run your multiqc commands with your conda env active
```

