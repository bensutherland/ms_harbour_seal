# Harbour Seal Project - Population Genetics
Analysis for harbour seal population genomics study. The pipeline is designed for the uses described in the manuscript and comes with no guarantees outside of that analysis.    
The pipeline primarily depends on the [stacks_workflow](https://github.com/enormandeau/stacks_workflow) repo of E. Normandeau (Labo Bernatchez) for genotyping, and [simple_pop_stats](https://github.com/bensutherland/simple_pop_stats) for population genetic analyses.         


### Requirements    
[Stacks2 (v2.3e)](http://catchenlab.life.illinois.edu/stacks/)     
[stacks_workflow](https://github.com/enormandeau/stacks_workflow)         
[simple_pop_stats](https://github.com/bensutherland/simple_pop_stats)      
[cutadapt](https://cutadapt.readthedocs.io/en/stable/index.html)    
[fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)      
[multiqc](https://multiqc.info/)     
[bwa](http://bio-bwa.sourceforge.net/)    
[samtools (v.1.9)](http://www.htslib.org/)    
[plink2](https://www.cog-genomics.org/plink/2.0/)           
[fineRADstructure](http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html)    


## 0. Obtain repos for analysis
Clone this repository, `stacks_workflow`, and `simple_pop_stats`, all at the same directory level. Change into the `stacks_workflow` directory for all commands, and use R for `simple_pop_stats`.       


## 1. Preparing Data
### a. Set up 
1. Put all raw data in `02-raw` using cp or cp -l    
2. Prepare the sample info file (see template in repo `sample_information.csv`). Note: tab-delimited, even though name is .csv.    
3. Download GenBank version reference genome: https://www.ncbi.nlm.nih.gov/genome/?term=Phoca+vitulina      

### b. Clean data
View raw data with fastqc and multiqc:    
```
mkdir 02-raw/fastqc_raw
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc_raw/ -t 5
multiqc -o 02-raw/fastqc_raw/ 02-raw/fastqc_raw
```

Prepare lane info
```
./00-scripts/00_prepare_lane_info.sh
```

Run cutadapt in order to trim off adapters, and remove any reads less than 50 bp
```
./00-scripts/01_cutadapt.sh 12
```

Demultiplex with two rxn enzymes in parallel over multiple CPUs
```
./00-scripts/02_process_radtags_2_enzymes_parallel.sh 80 nsiI mspI 14
```

Rename the files
```
./00-scripts/03_rename_samples.sh
```

FastQC de-multiplexed data
```
mkdir 04-all_samples/fastqc_demulti
fastqc 04-all_samples/*.fq.gz -o 04-all_samples/fastqc_demulti/ -t 14
multiqc -o 04-all_samples/fastqc_demulti/ 04-all_samples/fastqc_demulti
```


## 2. Alignment and genotyping
### a. Map reads against the reference genome
Index the genome
```
# (Only needed once) Change directory into the genome folder and index the genome 
bwa index -p GCA_004348235.1_GSC_HSeal_1.0_genomic ./GCA_004348235.1_GSC_HSeal_1.0_genomic.fna.gz
```

Align individual files against the genome:     
```
# First update the script below to point towards the directory containing your genome
GENOMEFOLDER="/home/ben/Documents/genomes"
GENOME="GCA_004348235.1_GSC_HSeal_1.0_genomic"
# note: also comment out the module load commands

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

#### Optional: remove low sequence depth samples 
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


### d. Genotype using reference genome-guided approach
Prepare and run gstacks      
```
# Prepare the population map
./00-scripts/04_prepare_population_map.sh

# Edit and run gstacks
# Update the NUM_CPU variable and run
./00-scripts/stacks2_gstacks_reference.sh

# Edit and run the ./00-scripts/stacks2_populations_reference.sh script
# For summary statistics (e.g., HOBS) and microhaplotypes
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 7 -r 0.7 \
    --min-maf 0.01 --ordered-export --radpainter \
    --vcf --hwe

# Make directory and move output 
mkdir 05-stacks/popn_out_mhaps
mv 05-stacks/populations* 05-stacks/popn_out_mhaps/

# Run populations again, for population genetic analysis (i.e., single-SNP per tag)
populations -P "$STACKS_FOLDER" -M "$INFO_FILES_FOLDER"/"$POP_MAP" \
    -t "$NUM_CPU" -p 7 -r 0.7 \
    --min-maf 0.01 \
    --ordered-export --genepop \
    --write-single-snp --hwe --fasta-loc

# Make directory and move output 
mkdir 05-stacks/popn_out_single_snp
mv 05-stacks/populations* 05-stacks/popn_out_single_snp/
```

## 3. Analysis of results
### Data setup
Change directory into `simple_pop_stats`.      

Copy in the colour file included with the harbour seal repo:       
`cp ../ms_harbour_seal/00_archive/harbour_seal_pops_colours.csv ./00_archive`       

Copy in the single-variant per-locus genepop:       
`cp ../stacks_workflow/05-stacks/popn_out_single_snp/populations.snps.genepop ./02_input_data/`     

Rename your file
`bhs_p<X>_r<0.X>_maf<0.0X>_20<XX-XX-XX>.gen` (note: customize as per true values for variables)         

Optional: copy in the multiple SNP per locus VCF (for individual inbreeding stat):      
`cp 05-stacks/popn_out_mhaps/populations.snps.vcf ../ms_harbour_seal/02_input_data/populations.snps_multi_per_locus.vcf`      


### Analysis
Analyze via `ms_harbour_seal/01_scripts/hs_popn_analysis_part1.R`, `*_part2.R` and `*_part3.R` sequentially.      

Note: if analyzing full dataset, make sure that the variable in the top of the first script is `dataset <- "all"` (default) instead of 'balanced'. The balanced dataset is to evaluate effect of sample size and read depth.     

In brief, this analysis will:       
1. Load data from the genepop; 
2. Evaluate missing data per individual; 
3. Global analysis to generate PCA and population-level FST;  
4. Atlantic coast only, filters (i.e., MAF, HWE, excess HOBS), then PCA, FST, per locus allele freq (myFreq.atl), per locus statistics, inter-individual relatedness;     
5. Pacific coast only, filters (i.e., MAF, HWE, excess HOBS), then PCA, FST, per locus allele freq (myFreq.atl), per locus statistics, inter-individual relatedness;     
6. Private alleles in Burrard Inlet outlier group; 
7. Per locus FST within individual contrasts filtered for low MAF (NBC vs. SOG, no outliers; ORE vs. SOG, no outliers; and SOG vs. outliers), and direct comparison of FST vals between dataset analyses;
8. Plotting of MAF distribution of markers in either Pacific or Atlantic datasets. 

Output is also saved as an image (`03_results/completed_analysis.RData`).       

There is another script to analyze MAF distributions using the _de novo_ analysis, and this is a short script `01_scripts/hs_popn_analysis_denovo.R`.       

## 4. Inter-individual relatedness 
#### i. Shared coancestry by related
This is now implemented in the main script (see above), and will generate PDFs of within-population relatedness, but also will generate text file of pairs that are the most highly related (Ritland > 0.1).      

#### ii. Shared haplotypes by fineRADstructure
Input: haplotype output from stacks populations as input to fineRADstructure         

Copy out radpainter input (SimpleMatrix) generated by stacks population module (above)         
`cp 05-stacks/popn_out_mhaps/populations.haps.radpainter ./../ms_harbour_seal/04_relatedness/`      

Change directory and run all commands from `ms_harbour_seal`       

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

## 5. Inbreeding coefficient per individual
Uses the vcf generated by the reference-guided genotyping, multiple SNP per locus, moved to the current repository.       
`vcftools --vcf 02_input_data/populations.snps_multi_per_locus.vcf --het --stdout > 03_results/per_ind_het.txt`        



## 5. Re-analysis of Liu et al. 2022
This section re-analyzes the data within [Liu et al. 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.16365) to compare with the present data, as well as to compare the results when using a reference genome or a _de novo_ approach.       

Clone a new `stacks_workflow` repository. Obtain raw data from NCBI SRA (e.g., via RunSelector and SRA Toolkit).      

### a. Data preparation and cleaning
Raw files have been demultiplexed, but still contain adapters. They also need to be trimmed to a constant length as per program specifications of `Stacks2` _de novo_. Put fastq data in `02-raw`, and compress using:      
`gzip 02-raw/*.fastq`      

Rename files to fit the population map:       
```
# create script to rename raw data      
grep -vE '^#' 01-info_files/sample_information.csv | awk '{ print "mv " $1 " " $3"_"$4".fq.gz" }' > 00-scripts/rename_raw_fastq_for_pop_map.sh        
# then add the shebang and chmod to make executable.       

# run the renaming script:     
cd 02-raw
./../00-scripts/rename_raw_fastq_for_pop_map.sh
cd ..
```

Inspect raw data:       
```
mkdir 02-raw/fastqc/       
fastqc 02-raw/*.fastq.gz -o 02-raw/fastqc/ -t 2       
multiqc -o 02-raw/fastqc/ 02-raw/fastqc       
```

Remove adapters and trim:      
`00-scripts/01_cutadapt.sh 3`        

Inspect trimmed data:      
```
mkdir 02-raw/trimmed/fastqc      
fastqc 02-raw/trimmed/*.fastq.gz -o 02-raw/trimmed/fastqc/ -t 2
multiqc -o 02-raw/trimmed/fastqc 02-raw/trimmed/fastqc      
```

Truncate the reads to require all reads to be 70 bp in constant length (as per de novo stacks requirements):      
```
mkdir 02-raw/standardized
./00-scripts/standardize_reads.sh
```
The output will be in the new standardized folder. 

Inspect standardized data:      
```
mkdir 02-raw/standardized/fastqc 
fastqc 02-raw/standardized/*.fq.gz -o 02-raw/standardized/fastqc/ -t 2
multiqc -o 02-raw/standardized/fastqc 02-raw/standardized/fastqc  
```

Move prepared files to `04-all_samples` (note: skipping renaming script)       
`mv 02-raw/trimmed/*.fastq.gz 04-all_samples/`      

#### Align and genotype
After editing the script, run the alignment       
`./00-scripts/bwa_mem_align_reads.sh 3`       

Aligned data will be in `04-all_samples`.       

```
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

## Genotyping de novo, coast-specific
#### Set up
Clone and rename as `stacks_workflow_trimming`         
Copy all raw data to `02-raw`.     
Copy sample information to `01-info_files`      

#### Region-specific analysis
Clone three copies of stacks workflow, for each denovo analysis. Copy in the sample info file:    
```
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_wc/01-info_files/
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_ec/01-info_files/
cp stacks_workflow_trimming/01-info_files/sample_information.csv  stacks_workflow_denovo_eu/01-info_files/
```

rsync in the raw data into the `04-all_data` folder:     
```
rsync -avzP 02-raw/standardized/END* ./../stacks_workflow_denovo_wc/04-all_samples/
rsync -avzP 02-raw/standardized/KOD* ./../stacks_workflow_denovo_wc/04-all_samples/
rsync -avzP 02-raw/standardized/BIC* ./../stacks_workflow_denovo_ec/04-all_samples/
rsync -avzP 02-raw/standardized/NFL* ./../stacks_workflow_denovo_ec/04-all_samples/
rsync -avzP 02-raw/standardized/ORK* ./../stacks_workflow_denovo_eu/04-all_samples/
rsync -avzP 02-raw/standardized/WNL* ./../stacks_workflow_denovo_eu/04-all_samples/
```

#### For each of denovo wc, ec, and eu:     
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
# for populations, use -p 2, -r 0.7, and --min-maf 0.01 as well as --hwe

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


### Admixture analysis ###
Change directory into `ms_harbour_seal` (this repository).    

Get the original single-SNP per locus VCF file output by Stacks2 populations module:     
`cp ../stacks_workflow_2023-02-27/05-stacks/popn_out_single_snp/populations.snps.vcf ./02_input_data/`       
i

As part of the population analysis, whitelists were made for individuals and loci for the Pacific, Atlantic, and both coasts together. Copy these whitelists into the admixture folder:    
`cp ../simple_pop_stats/03_results/*whitelist* 05_admixture/`     


#### Prepare the input files for admixture
Subset the input VCF to only contain the target individuals and loci for each of the three datasets:     

```
# Subset individuals
bcftools view -S 05_admixture/both_coasts_whitelist_inds.txt 02_input_data/populations.snps.vcf -Ov -o 05_admixture/populations.snps_filt_inds.vcf

# (Optional) view the inds remaining: 
bcftools query -l 05_admixture/populations.snps_filt_inds.vcf | wc -l

# Prepare a bed file to subset SNPs by using the following Rscript interactively:   
#  note: select which dataset is the goal
01_scripts/generate_ranges_file.R
#  note: this will output 05_admixture/*whitelist_loci_bed.txt

# Subset loci
intersectBed -a 05_admixture/populations.snps_filt_inds.vcf -b 05_admixture/both_coasts_whitelist_loci_bed.txt -header > 05_admixture/populations.snps_filt_inds_loci.vcf

# Convert to plink format
~/programs/plink2 --make-bed --vcf 05_admixture/populations.snps_filt_inds_loci.vcf --allow-no-sex --no-sex --no-parents --no-fid --no-pheno --allow-extra-chr --out 05_admixture/pvit_all

# Fix chr names issue
awk '{$1=0;print $0}' 05_admixture/pvit_all.bim > 05_admixture/pvit_all_fixed.bim
#awk '{$1=0;print $0}' 02_input_data/harbour_seal.bim > 02_input_data/harbour_seal_rem_chr.bim

# Overwrite
mv 05_admixture/pvit_all_fixed.bim 05_admixture/pvit_all.bim

```

#### Run admixture analysis

```
# Change into admixture directory
cd 05_admixture

# Run admixture to evaluate if everything is working with a single k:     
#admixture --cv 02_input_data/harbour_seal.bed 5 | tee 02_input_data/admixture_out.log

# Run in a loop with multiple k to test out different k levels:      
for K in 1 2 3 4 5 6 7 8; do admixture --cv pvit_all.bed $K | tee pvit_all_log${K}.out; done 

# View cv results
less pvit_all_log

# Collect CV error results into a file
cat 05_admixture/pvit_all_log* | grep 'CV error' > 05_admixture/pvit_all_CV_error.txt
```

Interactively use `01_scripts/admixture.R` to plot output.    
note: will need the GPS coordinate file.   

Copy `additional_file_s1_sample_metadata_2023-05-11.xlsx` from the additional files to `02_input_data`. Open it and save as a tab-delimited text file, of the same name but with the .txt suffix.  This will be used by the plotting script above.       

