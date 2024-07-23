#__________________________________________________________________________________________________________________________#
## Script:            simple HYE mix analysis
## raw data origin:   doi:10.1002/pmic.202300294
## raw data access:   https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=0f33717d84fd45b1a318ad40670022cc
## short description: HYE mix analysis using Exploris480 measurements
#__________________________________________________________________________________________________________________________#

# load libraries ----------------------------------------------------------
library(SpectroPipeR)
library(ggrepel)
library(tidyverse)
library(ggh4x)
library(scales)
library(ggrepel)
library(ggtext)
library(patchwork)

# setup parameters --------------------------------------------------------
params <- list(output_folder = "species_mix_analysis")


# define Spectronaut report file location ---------------------------------
file_path <- "HYE_Exploris480_SN19_Report_SpectroPipeR (Normal).tsv"

# perform the analysis
SpectroPipeR_analysis <- SpectroPipeR(file = file_path,
                                      parameter = params,
                                      condition_comparisons = cbind(c("HYE mix A","HYE mix B"))
)

# candidates XIC plot -----------------------------------------------------
# due to the bigger size the XIC SQLite files will not be deposit on github

XIC_plot_module(Spectronaut_report_path = file_path,
                Spectronaut_xicDB_path = "HYE_Exploris480_SN19_XIC-DBs",
                protein_groups = "P0A9X9",
                number_of_cores = 2,
                output_path = "species_mix_analysis/XIC_plots",
                export_csv_files = F)


# plot species mix analysis -----------------------------------------------

organisms_vector <- c(
  "Escherichia coli (strain K12)",
  "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)",
  "Homo sapiens"
)

# load statistics
stats_global <- read_delim("species_mix_analysis/06_statistics/8_sample_analysis/statistical_analysis.csv")


# filter for 2 peptides ---------------------------------------------------

stats_global <- stats_global %>% filter(n>=2)


# add orgnism to protein groups in statistics
organism <- read_delim(file = file_path) %>% 
  distinct(PG.ProteinGroups, PG.Organisms)

stats_global <- left_join(stats_global,
                          organism,
                          by = join_by(PG.ProteinGroups))

# filter for species in mix only
stats_global <- stats_global %>%
  filter(PG.Organisms %in% organisms_vector)

# get iBAQ int.
iBAQ <- read_csv("species_mix_analysis/05_processed_data/8_sample_analysis/iBAQ_protein_intensity_data_extracted_from_Spectronaut_summary.csv")

# get peptide int.
pep <- read_delim("species_mix_analysis/05_processed_data/8_sample_analysis/peptide_intensities_final_zero_values_replaced_with_half_minimal_intensity.csv")

# peptide CV
CV_pep <- pep %>% 
  group_by(R.Condition,PG.ProteinGroups) %>% 
  filter(n_distinct(EG.ModifiedPeptide)>=2) %>% # 2 peptide count filter
  group_by(R.Condition,PG.ProteinGroups,EG.ModifiedPeptide) %>% 
  summarise(CV = sd(peptide_intensity,na.rm=T)/mean(peptide_intensity,na.rm=T)) %>% 
  ungroup()

# filter for species in mix only
CV_pep <- left_join(CV_pep,organism,by = join_by(PG.ProteinGroups)) %>% 
  filter(PG.Organisms %in% organisms_vector)

# merge iBAQ mixB and stats 
stats_global <- left_join(stats_global,
                          iBAQ %>% filter(R.Condition=="HYE mix B"),
                          by = join_by(PG.ProteinGroups))


# benchmarking plot -------------------------------------------------------

Ecoli_threshold <- log2(1/3)
yeast_threshold <- log2(1.5)
human_threshold <- log2(1.2)
q_value_threshold <- 0.05

passed_threshold_Ecoli <- stats_global %>% 
  filter(PG.Organisms == organisms_vector[1]) %>% 
  mutate(passed_threshold = ifelse(slr<=Ecoli_threshold & p.fdr<0.05,
                                   yes = 1,
                                   no = 0))
passed_threshold_yeast <- stats_global %>% 
  filter(PG.Organisms == organisms_vector[2]) %>% 
  mutate(passed_threshold = ifelse(slr>=yeast_threshold & p.fdr<0.05,
                                   yes = 1,
                                   no = 0))
passed_threshold_human <- stats_global %>% 
  filter(PG.Organisms == organisms_vector[3]) %>% 
  mutate(passed_threshold = ifelse(abs(slr) >= human_threshold & p.fdr<=0.05,
                                   yes = 0,
                                   no = 1))


stats_global_final <- bind_rows(passed_threshold_Ecoli,passed_threshold_yeast,passed_threshold_human)



# summarize threshold passing PG counts 
stats_global_final_counts <- stats_global_final %>% 
  group_by(PG.Organisms,
           passed_threshold) %>% 
  summarise(count = n_distinct(PG.ProteinGroups)) %>% 
  ungroup()

# valid protein quantity count
valid_PG_quantity_count <- stats_global_final %>%
  group_by(
    PG.Organisms
  ) %>%
  summarise(distinct_PG = n_distinct(PG.ProteinGroups))

# max density estimate 
# denisty max calculate ratio of density max == mode
density_max <- stats_global_final %>%
  group_by( PG.Organisms) %>%
  summarise(
    density_max = max(density(slr, na.rm = T)$y),
    density_mode = 2^density(slr, na.rm = T)$x[which(density(slr, na.rm = T)$y == max(density(slr, na.rm = T)$y))]
  )

# median ------------------------------------------------------------------
org_median <- stats_global_final %>%
  group_by(PG.Organisms) %>%
  summarise(median = median(slr, na.rm = T))

org_median <- left_join(org_median, density_max)

# plot parameter
coords <- c(-4, 4)
protein_id_max <- 12000

# text parameter
axis_text_param <- 18
plot_tag_param <- 30
plot_stripped_text_param <- 14

facets_ncol <- 1


# scatter plot
plot_out <- ggplot(stats_global_final, aes(mean_iBAQ_intensities,slr,color = PG.Organisms)) +
  geom_point(alpha = 0.1) +
  geom_hline(yintercept = c(0, 1, -2), color = "darkgrey") +
  scale_x_log10() +
  facet_nested(~slr_ratio_meta) +
  annotation_logticks(sides = "l") +
  scale_color_brewer(palette = "Dark2") +
  stat_smooth(
    method = "loess",
    mapping = aes(group = PG.Organisms),
    color = "black",
    linetype = "dashed", se = F
  ) +
  labs(
    title = "ratio vs. iBAQ plot",
    x = "iBAQ (Mix B)",
    y = expression(log[2] ~ "protein quantity ratio")
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.title = element_markdown(size = 16, face ="bold"),
    strip.text = element_text(face = "bold", size = plot_stripped_text_param),
    axis.text = element_text(size = axis_text_param),
    axis.title.y = element_text(size = axis_text_param + 2),
    axis.title.x  = element_blank(),
    plot.tag = element_text(size = plot_tag_param, face = "bold")
  ) +
  coord_flip(ylim = coords) +
  guides(color = "none")

# density plot
dens_out <- ggplot(
  stats_global_final,
  aes(slr, fill = PG.Organisms)
) +
  geom_density(linewidth = 0, alpha = 0.8) +
  # coord_flip(xlim = coords)+
  coord_cartesian(xlim = coords) +
  theme_void() +
  scale_fill_brewer(palette = "Dark2") +
  scale_color_brewer(palette = "Dark2") +
  geom_vline(xintercept = c(0, 1, -2), color = "darkgrey") +
  guides(color = "none") +
  facet_wrap(~ slr_ratio_meta, ncol = facets_ncol) +
  geom_vline(data = org_median, mapping = aes(
    xintercept = median,
    color = PG.Organisms
  ), linetype = "dashed") +
  geom_label_repel(
    data = org_median,
    mapping = aes(
      x = log2(density_mode),
      y = density_max,
      color = PG.Organisms,
      label = round(density_mode, digits = 2)
    ),
    min.segment.length = 0.01, fill = "white", segment.color = "black"
  ) +
  guides(fill = "none") +
  theme(
    strip.text = element_blank(),
    axis.text = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_text(size = axis_text_param + 2),
    plot.tag = element_text(size = plot_tag_param, face = "bold")
  ) +
  scale_y_reverse()+
  labs(x = expression(log[2] ~ "ratio (A/B)"))

# CV plot
cv_plot_out <- ggplot(CV_pep,
                      aes(interaction(R.Condition,PG.Organisms), CV, fill = PG.Organisms)) +
  geom_violin() +
  geom_boxplot(outlier.colour = NA, width = 0.15, fill = "white") +
  facet_nested(R.Condition~ ., scales = "free_y") +
  scale_fill_brewer(palette = "Dark2") +
  coord_flip(ylim = c(0, 1.5)) +
  labs(
    title = "peptide CVs",
    tag = "",
    subtitle = "",
    y = "coefficient of variation of peptide intensity", x = ""
  ) +
  geom_hline(yintercept = 0.2) +
  geom_hline(yintercept = 0.1, linetype = "dashed") +
  theme_light() +
  guides(fill = "none") +
  theme(
    strip.background.x = element_blank(),
    strip.text.x = element_blank(),
    plot.title = element_markdown(size = 16, face ="bold"),
    strip.text.y = element_text(face = "bold", size = plot_stripped_text_param),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = axis_text_param, angle = 90,hjust = 1,vjust = 0.5),
    axis.title = element_text(size = axis_text_param + 2),
    plot.tag = element_text(size = plot_tag_param, face = "bold")
  )

# ID plot complex
ID_plot_out <- ggplot(valid_PG_quantity_count,
                      aes(1, distinct_PG, fill = PG.Organisms)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(
    palette = "Dark2",
    limits = c("Escherichia coli (strain K12)", "Homo sapiens", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"),
    labels = c("*Escherichia coli*", "*Homo sapiens*", "*Saccharomyces cerevisiae*")
  ) +
  labs(
    title = "protein groups with ≥ 2 peptides",
    tag = "",
    subtitle = "",
    y = "valid PG quantity count", x = "", fill = "organism"
  ) +
  # species specific ID count
  geom_text(
    position = position_stack(vjust = 0.5),
    angle = 0,
    aes(label = distinct_PG),
    size = 5
  ) +
  # global ID count
  geom_label(
    data = valid_PG_quantity_count %>%
      summarise(distinct_PG = sum(distinct_PG)),
    mapping = aes(1, distinct_PG, label = distinct_PG),
    size = 5, inherit.aes = F, hjust = -0.2
  ) +
  coord_flip(ylim = c(0, protein_id_max)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.text = element_markdown(size = axis_text_param),
    legend.title = element_text(size = axis_text_param + 5),
    plot.title = element_markdown(size = 16, face ="bold"),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    axis.text = element_text(size = axis_text_param),
    axis.text.y = element_blank(),
    axis.title = element_text(size = axis_text_param + 2),
    plot.tag = element_text(size = plot_tag_param, face = "bold")
  )+
  guides(fill="none")



significant_count <- ggplot(stats_global_final_counts %>% 
                              filter(passed_threshold==1), 
                            aes(x = 1, 
                                y = count, 
                                fill = PG.Organisms)) +
  geom_bar(stat= "identity") +
  scale_fill_brewer(
    palette = "Dark2",
    limits = c("Escherichia coli (strain K12)", "Homo sapiens", "Saccharomyces cerevisiae (strain ATCC 204508 / S288c)"),
    labels = c("*Escherichia coli*", "*Homo sapiens*", "*Saccharomyces cerevisiae*")
  ) +
  # species specific ID count
  geom_text(
    position = position_stack(vjust = 0.5),
    angle = 0,
    aes(label = count),
    size = 5
  ) +
  # global ID count
  geom_label(
    data = stats_global_final_counts %>% 
      filter(passed_threshold==1) %>% 
      summarise(count = sum(count)),
    mapping = aes(x = 1, count, label = count),
    size = 5, inherit.aes = F, hjust = -0.2
  ) +
  coord_flip(ylim = c(0, protein_id_max)) +
  theme_light() +
  theme(
    legend.position = "bottom",
    plot.subtitle = element_markdown(),
    plot.title = element_markdown(size = 16, face ="bold"),
    strip.text = element_markdown(size = 16),
    legend.text = element_markdown(size = axis_text_param),
    legend.title = element_text(size = axis_text_param + 5),
    axis.text = element_text(size = axis_text_param),
    axis.text.y = element_blank(),
    axis.title = element_text(size = axis_text_param + 2),
    plot.tag = element_text(size = plot_tag_param, face = "bold")
  )+
  labs(title = "count of protein groups passing thresholds",
       #subtitle = paste("**thresholds:**<br>
       #                 *E.coli* (theo. FC = -4): FC ≤", 1/(2^Ecoli_threshold)*-1," & q-value ≤",q_value_threshold,"<br>",
       #                 "*H.sapiens* (theo. FC = 1): abs. FC <", (2^human_threshold)," & q-value >",q_value_threshold,"<br>",
       #                 "*S.cerevisiea* (theo. FC = 2): abs. FC ≥", (2^yeast_threshold)," & q-value ≤",q_value_threshold,
       #                 ""),
       y = "PG count passing thresholds",
       fill = "",
       x = "",
       tag = "")+
  guides(fill=guide_legend(ncol=1))


# prepare final plot
final_plot <- plot_out + dens_out +cv_plot_out + ID_plot_out + significant_count+
  plot_layout(heights = c(4, 1.5,3, 1,1)) +
  plot_annotation(
    theme = theme(
      plot.title = element_text(size = 30, face = "bold"),
      plot.subtitle = element_text(size = 18)
    ),
    title = paste0("Exploris480 - HYE mixtures"),
    subtitle = "protein level - FDR = 0.01, Q-value < 0.01"
  )
ggsave(
  filename = paste0("species_mix_analysis/HYE_benchmarking_detailed_plot_panel.png"),
  plot = final_plot, device = "png", height = 20, width = 7
)

