## Harbour seal simple pop stats and loci selection
# C Rycroft

## 00_setup ----

# additional libraries 
library(tidyverse)

# source simple_pop_stats.R
# species prompt: [8] Harbour seal
 
#  create custom output folder
#set_output_folder(result_path = '~/02_r_code/01_phvi/04_phvi_AID20210047_output/')

# note: how are we going to deal with this as some files are piped into the custom file, 
#       but many are stored in the simple pop stats '03_results'

# note: '## '##_description' correlates to headers in simple pop stats pipeline

 # **end section ----


## 01_loading_data ----
# load genepop data file
load_genepop(datatype = 'SNP')
# note: baseline harbour seal; pops = 7; maf = 0.01

# create vector of unique populations
pops_unique <- unique(pop(obj))
# note: script edit to make 'locus_stats' a global object

# **end section ----


## 03_characterize_data ----
# find the number of samples, markers, alleles and population sample size. store in custom output directory
characterize_genepop(df = obj 
                     , pdf_width = 30 
                     , pdf_height = 6 
                     , cex_names = 0.3 
                     , N = 30 # note: draws horizontal line across plot
                     ) 


# calculate per-locus fst, per-locus stats, and write to 03_results (labelled with date)
per_locus_stats(data = obj)

# **end section ----


## 04_drop_loci ----
# required: this needs changes obj to obj_filt
drop_loci(drop_monomorphic = TRUE
          #, drop_file = ''
          ) 

# **end section ----


## 05_drop_pops ----
#create drop file
 # drop_pop_file <- c('EQB_120'
 #                    , 'LAB_107'
 #                    , 'NFL_124'
 #                    )
 
 # write.table(drop_pop_file
 #             , file = 'C:/Users/RycroftC/Documents/02_r_sandbox/01_phvi/00_simple_pop_stats/drop_pop_file_west_coast.txt'          
 #             #, sep = ''
 #             )
 # note: need to go into txt file and remove the row numbers and header that gets written. 

# remove populations based on sample size
drop_pops(df = obj
          , drop_by_pop_size = FALSE 
          , min_indiv = NULL
          , drop_file = 'C:/Users/RycroftC/Documents/02_r_sandbox/01_phvi/00_simple_pop_stats/drop_pop_file_west_coast.txt'
          , drop_unused = FALSE
          , names.dat.FN = NULL 
          , drop_indivs = NULL
          )
# action item: change the code to use a vector instead of a file. 
# note: this script is a pain. when writing from the vector to the tab delimited file, need to go in and edit. 

# **end section ----


## 06_genetic_differentiation ----
calculate_FST(format = 'genind'
              , dat = obj_pop_filt
              , separated = FALSE
              )

# **end section ----


## 07_build_tree # note: not needed for analyses ----
# make_tree(bootstrap = TRUE
#           , boot_obj = obj_pop_filt
#           , nboots = 10000
#           , dist_metric = 'edwards.dist'
#           , separated = FALSE
#          )

# **end section ----


## 08_run_multidimensional_scaling_techniques ----
# conduct PCA
pca_from_genind(data = obj_pop_filt
                , PCs_ret = 3
                , plot_eigen = TRUE
                , plot_allele_loadings = TRUE
                #, colour_file = NULL
                ) # note: manually set colours for each pop; see lab book pg 20 - 21


# conduct DAPC
dapc_from_genind(data = obj_pop_filt
                 , plot_allele_loadings = TRUE
                 #, colour_file = NULL
                 ) # note: manually set colours for each pop; see lab book pg 20 - 21
# **end section ----


## 09_calculate_relatedness ----
# convert data from genind to relatedness
relatedness_calc(data = obj_pop_filt
                 , datatype = 'SNP'
                 )


# plot result
# relatedness_plot(file = '03_results/kinship_analysis_2021-11-05.Rdata'
#                  , same_pops = TRUE
#                  , plot_by = 'codes'
#                  )


# population marker HWE evaluation summary
hwe_eval(data = obj_pop_filt
         , alpha = 0.01
         )

# **end section ----


## 12_generate_allele_frequency_table ----
# calculate allele frequencies; split by pop
calculate_allele_freq(data = obj_pop_filt)


# calculate allele frequencies; global
# create new object and alter for one population code
global_pop_filt <- obj_pop_filt
pop(global_pop_filt) <- rep('all', times = length(pop(global_pop_filt)))

# convert to genpop, calculate af
global.genpop <- genind2genpop(x = global_pop_filt)
global_freq <- makefreq(x = global.genpop)

# transpose
global_freq <- t(global_freq)

# make df
global_freq <- as.data.frame(global_freq, stringsAsFactors = F)

global_freq <- global_freq %>% 
    rownames_to_column('allele.id') %>% 
    rename(global_allele_freq = all)


# **end section ----


### leaving simple pop stats script; development of harbour seal loci selection


## relatedness ----
# create relatedness data frame from 
head(output$relatedness)

relatedness.df <- output$relatedness
unique(relatedness.df$group)


# identify related individuals within each population
relatedness.df.pops <- as.data.frame(relatedness.df) %>% 
        filter(group == 'SOSO' & ritland > 0.09
               | group == 'NBNB' & ritland > 0.09
               | group == 'CACA' & ritland > 0.09
               | group == 'OROR' & ritland > 0.09)
# note: bjgs originally had this at 0.1
# result: OROR samples were retained at a value of 0.09, but were lost at 0.1

unique(relatedness.df.pops$group)


# split samples into populations and identify relatives by applying more stringent filter
# note: filtered relatedness.df.pops by 0.1, but having looked at relatedness_plot output, 0.09 seems better to use.
relatedness.NBC.sibs <- relatedness.df %>% 
        filter(group == 'NBNB')

relatedness.SOG.sibs <- relatedness.df %>% 
        filter(group == 'SOSO')

relatedness.ORE.sibs <- relatedness.df %>% 
        filter(group == 'OROR')

relatedness.CAL.sibs <- relatedness.df %>% 
        filter(group == 'CACA')


# create sibling objects
sib_pairs <- relatedness.df.pops %>%
    dplyr::select('ind1.id'
                  , 'ind2.id')
    
removed_sib <- sib_pairs %>%
         dplyr::select(ind1.id) %>% 
         distinct()


# remove sibling 1 from related pair 
relatedness.plot.filt <- relatedness.df %>% 
        filter(!(ind1.id %in% removed_sib$ind1.id)) %>% 
        filter(group == 'NBNB' |
                       group == 'SOSO' |
                       group == 'OROR' |
                       group == 'CACA'
                       )

unique(relatedness.plot.filt$ind1.id) # note: 93 animals unfiltered; 59 filtered


# plot relatedness results
pop.plot <- ggplot(relatedness.plot.filt, 
                   aes(x = group
                       , y = ritland
                       , fill = group
                       )) +
        geom_boxplot() +
        labs(x = 'Population'
             , y = 'Ritland relatedness estimate'
             ) +
        scale_fill_discrete(name = 'Population'
                            , breaks = c('CACA'
                                         , 'NBNB'
                                         , 'OROR'
                                         , 'SOSO'
                                         )
                            , labels = c('California'
                                         , 'Northern\nBritish Columbia'
                                         , 'Oregon'
                                         , 'Strait of\nGeorgia'
                                         )) +
        theme_minimal()

#pdf('C:/Users/RycroftC/Documents/02_r_sandbox/01_phvi/00_simple_pop_stats/03_results/01_west_coast_panel/relatedness_ritland_no_sib_2021-11-10.pdf')
pop.plot
#dev.off()

# **end section ----


## hobs filter ----
# identify top 300 differentiating loci with Hobs less than 0.5
# create fixed variable for heterozygousity 
het_level <- 0.5

hobs_loci <- as.data.frame(locus_stats) %>% 
    rownames_to_column('locus') %>% 
    filter(Hobs < het_level) %>% # note: may need to do this differently and filter out Hobs later
    arrange(desc(Fst))
# note: removed 'top_loci' as more filtering required before top diff loci can be identified

removed_hobs_loci <- as.data.frame(locus_stats) %>% 
    rownames_to_column('locus') %>% 
    filter(Hobs > het_level)

# **end section ----    


## hardy-weinberg ----
# create fixed variable for chi^2
prob <- 0.05 # note: corresponds to p = 0.05

# split populations, remove monomorphic loci (?) and assign hwe status
# note: work toward putting this into a loop
hwe_NBC_df <- as.data.frame(hwe.dat$NBC_115) %>% 
    rownames_to_column('locus') %>% 
    filter(!(`chi^2` == 0 & df == 0) # note: remove monomorphic sites
           , (locus %in% hobs_loci$locus)) %>% # note: remove loci Hobs> 0.5
    mutate(hwe_NBC = case_when(`Pr(chi^2 >)` > prob ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`) 

hwe_SOG_df <- as.data.frame(hwe.dat$SOG_141) %>% 
    rownames_to_column('locus') %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (locus %in% hobs_loci$locus)) %>% 
    mutate(hwe_SOG = case_when(`Pr(chi^2 >)` > prob ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)

hwe_ORE_df <- as.data.frame(hwe.dat$ORE_120) %>% 
    rownames_to_column('locus') %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (locus %in% hobs_loci$locus)) %>% 
    mutate(hwe_ORE = case_when(`Pr(chi^2 >)` > prob ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)

hwe_CAL_df <- as.data.frame(hwe.dat$CAL_130) %>% 
    rownames_to_column('locus') %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (locus %in% hobs_loci$locus)) %>% 
    mutate(hwe_CAL = case_when(`Pr(chi^2 >)` > prob ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)


# **end section ----


## allele_frequency_table ----
allele_freq <- as.data.frame(freq.df) %>% 
        rownames_to_column('allele.id') %>% 
        mutate(locus = gsub('.{3}$', '', allele.id))


# **end section ----


## analysis data frame ----
# combine hobs_filtered_loci, allele_frequencies, hwe 
all_loci_df<- full_join(as.data.frame(hobs_loci)
                          , allele_freq
                          , by = 'locus') %>%
    left_join(.
                  , hwe_NBC_df %>% 
                          dplyr::select('locus'
                                        , 'hwe_NBC')
                  , by = 'locus') %>%
    left_join(.
                  , hwe_SOG_df %>% 
                          dplyr::select('locus'
                                        , 'hwe_SOG')
                  , by = 'locus') %>% 
    left_join(.
                  , hwe_ORE_df %>% 
                          dplyr::select('locus'
                                        , 'hwe_ORE')
                  , by = 'locus') %>% 
    left_join(.
                  , hwe_CAL_df %>% 
                          dplyr::select('locus'
                                        , 'hwe_CAL')
                  , by = 'locus') %>% 
        dplyr::select(locus
                      , allele.id
                      , everything()) 


# remove any loci with NA in hwe_pop column or in hwe dis-equil
# note: depending on how loci selection goes, consider removing the dis-equil filter
all_loci_filt <- all_loci_df %>% 
    filter_at(vars(hwe_NBC
                   , hwe_SOG
                   , hwe_ORE
                   , hwe_CAL)
              , all_vars(!is.na(.))) %>% 
    filter(!hwe_NBC == 'non-equil'
           , !hwe_SOG == 'non-equil' 
           , !hwe_ORE == 'non-equil'
           , !hwe_CAL == 'non-equil') 

# **end section ----


## minor allele frequency ----

# set fixed variables
allele_freq_coeff <- 0.5
maf_filter <- 0.01

# join hwe filtered loci with global maf
maf_loci <- all_loci_filt %>% 
    group_by(locus) %>% 
    rowwise() %>% 
    left_join(.
              , global_freq
              , 'allele.id') %>% 
    mutate(global_maf = case_when(global_allele_freq >= allele_freq_coeff ~ 'ref'
                              , TRUE ~ 'alt' )) %>% 
    mutate(maf_NBC = case_when(NBC_115 >= 0.5 ~ 'ref' # note: identify pop specific ref and alt alleles
                               , TRUE ~ 'alt') 
           , maf_SOG = case_when(SOG_141 >= 0.5 ~ 'ref'
                                 , TRUE ~ 'alt')
           , maf_ORE = case_when(ORE_120 >= 0.5 ~ 'ref'
                                 , TRUE ~ 'alt')
           , maf_CAL = case_when(CAL_130 >= 0.5 ~ 'ref'
                                 , TRUE ~ 'alt')
    ) %>% 
    mutate(alt_allele = case_when(maf_NBC == 'alt' # note:create Y/N column of allele stat
                                  & maf_SOG == 'alt'
                                  & maf_ORE == 'alt'
                                  & maf_CAL == 'alt' ~ 'Y'
                                  , TRUE ~ 'N'
                                  )) %>% 
    filter(!(maf_NBC <= maf_filter | # note: noting removed as imported data filtered at 0.01
                 maf_SOG <= maf_filter |
                 maf_ORE <= maf_filter |
                 maf_CAL <= maf_filter )) %>% 
    dplyr::select(locus
                  , allele.id
                  , Fit
                  , Fst
                  , Fis
                  , Hobs
                  , Hexp
                  , alt_allele
                  , NBC_115
                  , hwe_NBC
                  , maf_NBC
                  , SOG_141
                  , hwe_SOG
                  , maf_SOG
                  , ORE_120
                  , hwe_ORE
                  , maf_ORE
                  , CAL_130
                  , hwe_CAL
                  , maf_CAL
                  , everything())

# identify top loci fst
top_fst_df <- maf_loci %>% 
    ungroup() %>% 
    distinct(locus
             , .keep_all = TRUE) %>%
    arrange(desc(Fst)) %>% 
    slice_max(Fst
              , n = 300)

# identify top loci heterozygous 
top_hobs_df <- maf_loci %>% 
    ungroup() %>% 
    arrange(desc(Hobs)) %>% 
    distinct(locus, .keep_all = TRUE) %>%
    slice_max(Hobs
              , n = 300)

# join top fst and hobs designation to main df
maf_loci_df <- maf_loci %>% 
    mutate(top_fst_loci = 'N' # note: create columns designating top loci
           , top_hobs_loci = 'N') %>% 
    mutate(top_fst_loci = case_when(locus %in% top_fst_df$locus 
                                    & allele.id %in% top_fst_df$allele.id ~ 'Y'
                                   , TRUE ~ top_fst_loci)
           , top_hobs_loci = case_when(locus %in% top_hobs_df$locus 
                                       & allele.id %in% top_hobs_df$allele.id ~ 'Y'
                                       , TRUE ~ top_hobs_loci))
    
# ** end section ----


## removed loci ----
# hobs 
removed_hobs_loci 
# note: full code above


# hwe 
removed_hwe_loci <- rbind(as.data.frame(hwe.dat$NBC_115)
                          , as.data.frame(hwe.dat$SOG_141)
                          , as.data.frame(hwe.dat$ORE_120)
                          , as.data.frame(hwe.dat$CAL_130)) %>%
    rownames_to_column('locus') %>% 
    filter((`chi^2` == 0 & df == 0)
           | (locus %in% removed_hobs_locus$locus))

removed_mono_loci <- all_loci_df %>% 
    filter_at(vars(hwe_NBC
                   , hwe_SOG
                   , hwe_ORE
                   , hwe_CAL)
              , all_vars(is.na(.)))

removed_non_hwe_loci <- all_loci_df %>% 
    filter(hwe_NBC == 'non-equil'
           | hwe_SOG == 'non-equil' 
           | hwe_ORE == 'non-equil'
           | hwe_CAL == 'non-equil') 


# maf
removed_maf_loci <- all_maf_loci_pop %>% 
    filter(maf_NBC == 'ref'
           | maf_SOG == 'ref'
           | maf_ORE == 'ref'
           | maf_CAL == 'ref') 


# **end selection ----

## write to csv ----
write.csv (relatedness.df
           , '03_results/01_west_coast_panel/relatedness_west_coast_df_2021-11-05.csv')

write.csv (relatedness.df.pops
           , '03_results/01_west_coast_panel/relatedness_west_coast_df_pops_2021-11-05.csv')

write.csv (relatedness.NBC.sibs
           , '03_results/01_west_coast_panel/relatedness_NBC_sibs_2021-11-05.csv')

write.csv (relatedness.SOG.sibs
           , '03_results/01_west_coast_panel/relatedness_SOG_sibs_2021-11-05.csv')

write.csv (relatedness.ORE.sibs
           , '03_results/01_west_coast_panel/relatedness_ORE_sibs_2021-11-05.csv')

write.csv (relatedness.CAL.sibs
           , '03_results/01_west_coast_panel/relatedness_CAL_sibs_2021-11-05.csv')

write.csv(global_freq, '03_results/01_west_coast_panel/global_freq_2021-11-18.csv')

write.csv(sib_pairs, '03_results/01_west_coast_panel/sib_pairs_2021-11-18.csv')

write.csv(top_fst_df, '03_results/01_west_coast_panel/top_fst_2021-11-24.csv')

write.csv(top_hobs_df, '03_results/01_west_coast_panel/top_hobs_2021-11-24.csv')


# **end section ----

## data checks ----
# check to see if the same loci are being selected b/t chi^2 and Pr(chi^2) ---- 
# note: code is now variable based above, but not reflected here
hwe.NBC_pr <- as.data.frame(hwe.dat$NBC_115) %>% 
    rownames_to_column('locus') %>% 
    filter(!(`chi^2` == 0 & df == 0)) %>%
    mutate(hwe_NBC = case_when(`Pr(chi^2 >)` > 0.05 ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(desc(`chi^2`))

loci.selection.method <- anti_join(hwe.NBC
                                   , hwe.NBC_pr
                                   , by = 'locus')

# hwe using chi^2 ----
chi_coeff <- 3.84 # note: corresponds to p = 0.05

# split populations, remove monomorphic loci (?) and assign hwe status
# note: work toward putting this into a loop
hwe_NBC_df <- as.data.frame(hwe.dat$NBC_115) %>% 
    rownames_to_column('loci') %>% 
    filter(!(`chi^2` == 0 & df == 0) # note: remove monomorphic sites
           , (loci %in% hobs_loci$loci)) %>% # note: remove loci Hobs> 0.5
    mutate(hwe_NBC = case_when(`chi^2` < chi_coeff ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`) 

hwe_SOG_df <- as.data.frame(hwe.dat$SOG_141) %>% 
    rownames_to_column('loci') %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (loci %in% hobs_loci$loci)) %>% 
    mutate(hwe_SOG = case_when(`chi^2` < chi_coeff ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)

hwe_ORE_df <- as.data.frame(hwe.dat$ORE_120) %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (loci %in% hobs_loci$loci)) %>% 
    mutate(hwe_ORE = case_when(`chi^2` < chi_coeff ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)

hwe_CAL_df <- as.data.frame(hwe.dat$CAL_130) %>% 
    rownames_to_column('loci') %>% 
    filter(!(`chi^2` == 0 & df == 0) 
           , (loci %in% hobs_loci$loci)) %>% 
    mutate(hwe_CAL = case_when(`chi^2` < chi_coeff ~ 'equil'
                               , TRUE ~ 'non-equil')) %>%
    arrange(`chi^2`)

# * end section ----

