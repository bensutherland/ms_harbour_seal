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
`01_scripts/fineRADstructurePlot.R` (follow instructions here)      
`01_scripts/FinestructureLibrary.R`     
This will produce plots in the working directory.  

