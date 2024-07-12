#__________________________________________________________________________________________________________________________#
## Script:            simple cancer cohort study analysis of PXD047839 raw data
## raw data origin:   https://pubs.acs.org/doi/10.1021/acs.jproteome.3c00646
## raw data access:   https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD047839
## short description: In the study, plasma samples from 20 late-stage lung cancer patients and 20 control individuals
##                    were analyzed using the timsTOF HT mass spectrometer. The samples were processed with the Proteograph
##                    Product Suite, which employs five distinct nanoparticles (NP1−5) to compress the dynamic range of the
##                    plasma proteome, enabling deeper proteome profiling. Specifically, the NP2 nanoparticle was used for 
##                    SEER enrichment of plasma. Here the Evosep/timsTOF HT measurements of NP2 were used.
#__________________________________________________________________________________________________________________________#



# load libraries ----------------------------------------------------------
library(SpectroPipeR)
library(gprofiler2)
library(tidyverse)

# setup parameters --------------------------------------------------------
params <- list(output_folder = "PXD047839_cancer_cohort_NP2_SEER_SN19")


# define Spectronaut report file location ---------------------------------
file_path <- "PXD047839_cancer_cohort_NP2_SEER_SN19_Report_SpectroPipeR (Normal).tsv"

# perform the analysis
SpectroPipeR_analysis <- SpectroPipeR(file = file_path,
                                      parameter = params,
                                      condition_comparisons = cbind(c("cancer","control"))
)



# candidates XIC plot -----------------------------------------------------
# due to the bigger size the XIC SQLite files will not be deposit on github

XIC_plot_module(Spectronaut_report_path = file_path,
                Spectronaut_xicDB_path = "PXD047839_cancer_cohort_NP2_SEER_SN19_XIC-DBs",
                protein_groups = "P0DJI8",
                number_of_cores = 2,
                output_path = "PXD047839_cancer_cohort_NP2_SEER_SN19/XIC_plots",
                export_csv_files = F)

# simple Gprofiler2 analysis using the statistical analysis of SpectroPipeR as basis
# https://stemicha.github.io/SpectroPipeR/articles/a07_Gprofiler2_code_suggestion.html

# load SpectroPipeR statistics --------------------------------------------

#path to statistical_analysis.csv
stats <- read_csv("PXD047839_cancer_cohort_NP2_SEER_SN19/06_statistics/40_sample_analysis/statistical_analysis.csv")

# filter statistics -------------------------------------------------------
stats_filtered <- stats %>%
  # filter for ≥2 peptides
  filter(n>=2) %>%
  # filter for abs. FC & q-value
  filter(fold_change_absolute >= 1.5 & p.fdr <= 0.05)


# get statistical comparisons ---------------------------------------------
comparison <- unique(stats_filtered$slr_ratio_meta)


# generate output folder --------------------------------------------------
gprofiler_output_dir <- "PXD047839_cancer_cohort_NP2_SEER_SN19/gprofiler_query_Results"
dir.create(gprofiler_output_dir,showWarnings = F)


# perform Gprofiler analysis in a for loop over comparisons ---------------
for(i in comparison){
  
  # filter data for comparison
  grpofiler_input <- stats_filtered %>%
    filter(slr_ratio_meta %in% i)
  
  #
  # run gprofiler
  # query         - character vector of protein IDs
  # significant   - whether all or only statistically significant results should be returned.
  # organism      - Organism names are constructed by concatenating the first letter of the name and
  #                 the family name. Example: human - 'hsapiens', mouse - 'mmusculus'
  
  grpofiler_res<- gost(query = grpofiler_input$PG.ProteinGroups,
                       significant = T, # only statistically significant results should be returned
                       organism = "hsapiens" # select the right organism
  )
  
  # write grpofiler results container
  write_rds(x = grpofiler_res,file = paste0(gprofiler_output_dir,
                                            str_replace_all(i,"/","_vs_"),".rds"))
  # write grpofiler results table
  write_csv(x = as_tibble(grpofiler_res$result),
            file = paste0(gprofiler_output_dir,str_replace_all(i,"/","_vs_"),".csv"))
  
  # write grpofiler results plot
  publish_gostplot(
    p = gostplot(grpofiler_res, interactive = FALSE),
    highlight_terms = NULL,
    filename = paste0(gprofiler_output_dir,str_replace_all(i,"/","_vs_"),"__plot.png"),
    width = NA,
    height = NA
  )
  
  # write grpofiler results table
  publish_gosttable(
    gostres = grpofiler_res,
    filename = paste0(gprofiler_output_dir,str_replace_all(i,"/","_vs_"),"__table.png"),
    highlight_terms = grpofiler_res$result$term_id[which(grpofiler_res$result$p_value<1E-16)]
  )
}

