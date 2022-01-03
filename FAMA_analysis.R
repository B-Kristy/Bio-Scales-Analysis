# Load Libraries
library(readxl)
library(reshape2)
library(ggplot2)
library(dplyr)
library(viridis)
library(cowplot)
library(vegan)
library(tidyverse)
library(ggpubr)

setwd("~/Dropbox/BioScales/FAMA_data/")
# Import FAMA category Data
FAMA_category <- read_excel("FAMA_category.xlsx")
FAMA_functions <- read_excel("FAMA_functions.xlsx")
# Rhizosphere vs. Soil Analysis ----------------------------------------------
# Develop distance matrix
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype", "pH",
                   "Ca (mg/kg (ppm))", "K (mg/kg (ppm))", "Mg (mg/kg (ppm))",
                   "Mn (mg/kg (ppm))", "P (mg/kg (ppm))", "Zn (mg/kg (ppm))", "NH4-N (mg/kg)",
                   "NO3-N (mg/kg)", "C (%)", "N (%)", "MC (%)", "Compartment")

idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`,
                                           pH, `Ca (mg/kg (ppm))`, `K (mg/kg (ppm))`, `Mg (mg/kg (ppm))`,
                                           `Mn (mg/kg (ppm))`, `P (mg/kg (ppm))`, `Zn (mg/kg (ppm))`, `NH4-N (mg/kg)`,
                                           `NO2-N (mg/kg)`, `NO3-N (mg/kg)`, `C (%)`, `N (%)`, `MC (%)`, `Compartment`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")
perm = 9999
set.seed(309125)

# Plot dbRDA
rda_FAMA <- dbrda(dist_FAMA ~ Compartment + `Site Location`+ Compartment*`Site Location`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, shape = `Site Location`, fill = `Compartment`), 
             size = 3, alpha = .65) +
  ggtitle("Nitrogen Gene Composition") +
  scale_fill_manual(values = c("darkorchid1", "goldenrod2"), name = "Compartment", 
                    labels = c(expression("Rhizosphere"),
                               expression("Soil"))) +
  scale_shape_manual(values = c(21,24), name = "Site Location",
                     labels = c(expression("Clatskanie"), 
                                expression("Corvallis"))) +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12))+
        guides(fill = guide_legend(override.aes = list(pch = 21), title = "Microbiome Compartment")) -> p_dbrda_FAMA

ggsave(filename=paste("dbRDA_soil_vs_rhizo.tiff", sep=""), plot=p_dbrda_FAMA, 
       width=5, height=4, dpi=600)   

# Measure dispersion

dispersion <- betadisper(dist_FAMA, FAMA_metadata$Compartment, type = c("centroid"), bias.adjust = TRUE, sqrt.dist = FALSE, add = FALSE)
anova(dispersion)

set.seed(1221990)
permutest(dispersion, pairwise = TRUE, permutations = 9999)
(dispersion_16S_HSD <- TukeyHSD(dispersion_16S))

# Plot Dispersion

data <- data.frame(Distance_to_centroid=dispersion$distances, Group=dispersion$group)
data$Compartment <- FAMA_metadata$Compartment
data$`Site Location` <- FAMA_metadata$`Site Location`
my_comparisons <- list(c("Rhizosphere", "Soil"))

groups <- dispersion$group
ggplot(data=data, aes(x = Group, y=Distance_to_centroid, color = `Site Location`)) + 
  geom_jitter(size =2.25, position = position_jitter(0.1)) +
  scale_y_continuous(name = "Distance to Centroid") +
  geom_vline(xintercept = c(1.5, 2.5), linetype = "dotted") +
  stat_compare_means(comparisons = my_comparisons, method ="t.test") +
  panel_border() +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  theme(strip.text.y = element_text(size = 12), legend.text = element_text(hjust = 0),
        plot.title = element_text(size=12, face="bold"),
        legend.title = element_text(size = 12)) -> disper_plot

ggsave(filename=paste("disper_soil_vs_rhizo.tiff", sep=""), plot=disper_plot, 
       width=5, height=4, dpi=600)

# Plot Heatmap
FAMA_functions_melt <- melt(FAMA_functions)
FAMA_functions_ammonification <- subset(FAMA_functions_melt, variable == "NrfA"|
                                          variable == "NrfB"| variable == "NrfC"|
                                          variable == "NrfD"| variable == "NrfH")
FAMA_functions_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "AmoA_PmoA"|
                                              variable == "AmoB_PmoB"| variable == "AmoC_PmoC")
FAMA_functions_anaerobic_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "HAO"|
                                                        variable == "Hzo"| variable =="HzsA"|
                                                        variable == "HzsB"| variable =="HzsC")
FAMA_functions_denitrification <- subset(FAMA_functions_melt, variable == "NirB3"|
                                           variable == "NirC"| variable == "NirK"|
                                           variable == "NirM"| variable == "NirN"|
                                           variable == "NirS"| variable == "NirT"|
                                           variable == "NosZ"| variable == "cNor-C"|
                                           variable == "cNorB_qNor")
FAMA_functions_nitrate_assimilatory_reduction <- subset(FAMA_functions_melt, variable == "NasA"|
                                                          variable == "NasB")
FAMA_functions_nitrate_dissimilatory_reduction <- subset(FAMA_functions_melt, variable == "NapA"|
                                                           variable == "NapB"| variable == "NapC"|
                                                           variable == "NapD"| variable == "NapE"|
                                                           variable == "NapF"| variable == "NapG"|
                                                           variable == "NapH"| variable == "NapK"|
                                                           variable == "NapL"| variable == "NarC"|
                                                           variable == "NarG_NxrA"| variable == "NarH_NxrB"|
                                                           variable == "NarI"| variable == "NarJ")
FAMA_functions_nitrite_assimilation <- subset(FAMA_functions_melt, variable == "NasI"|
                                                variable == "NasJ"| variable == "NirA"|
                                                variable == "NirB"| variable == "NirD"|
                                                variable == "NirU")
FAMA_functions_nitrogen_fixation <- subset(FAMA_functions_melt, variable == "AnfG_VnfG"|
                                             variable == "NifB"| variable == "NifD_AnfD_VnfD"|
                                             variable == "NifH_AnfH_VnfH"| variable == "NifK_AnfK_VnfK")
FAMA_functions_urease <- subset(FAMA_functions_melt, variable == "UreA"|
                                  variable == "UreB"| variable == "UreC")

FAMA_functions_ammonification$category <- "Ammonification"
FAMA_functions_ammonium_oxidation$category <- "Ammonium oxidation"
FAMA_functions_anaerobic_ammonium_oxidation$category <- "Anaerobic ammonium oxidation"
FAMA_functions_denitrification$category <- "Denitrification"
FAMA_functions_nitrate_assimilatory_reduction$category <- "Nitrate assimilatory reduction"
FAMA_functions_nitrate_dissimilatory_reduction$category <- "Nitrate dissimilatory reduction"
FAMA_functions_nitrite_assimilation$category <- "Nitrite assimilation"
FAMA_functions_nitrogen_fixation$category <- "Nitrogen Fixation"
FAMA_functions_urease$category <- "Urease"

FAMA_functions_melt <- rbind(FAMA_functions_ammonification,
                             FAMA_functions_ammonium_oxidation,
                             FAMA_functions_anaerobic_ammonium_oxidation,
                             FAMA_functions_denitrification,
                             FAMA_functions_nitrate_assimilatory_reduction,
                             FAMA_functions_nitrate_dissimilatory_reduction,
                             FAMA_functions_nitrite_assimilation,
                             FAMA_functions_nitrogen_fixation,
                             FAMA_functions_urease)


# Total heatmap (Same Scale)
FAMA_functions_melt %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_compartment

# Subsetted heatmaps (Different scale for each category)
FAMA_functions_ammonification %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonification_heatmap

FAMA_functions_ammonium_oxidation %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonium_oxidation_heatmap

FAMA_functions_anaerobic_ammonium_oxidation %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Anaerobic ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_anaerobic_ammonium_oxidation_heatmap

FAMA_functions_denitrification %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Denitrification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_denitrification_heatmap

FAMA_functions_nitrate_assimilatory_reduction %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate assimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_assimilatory_reduction_heatmap

FAMA_functions_nitrate_dissimilatory_reduction %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate dissimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_dissimilatory_reduction_heatmap

FAMA_functions_nitrite_assimilation %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrite assimilation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrite_assimilation_heatmap

FAMA_functions_nitrogen_fixation %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrogen fixation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrogen_fixation_heatmap

FAMA_functions_urease %>%
  ggplot(aes(`Compartment`, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Urease") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_urease_heatmap

ggsave(filename=paste("ammonification_heatmap.tiff", sep = ""), plot = p_ammonification_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("ammonium_oxidation_heatmap.tiff", sep = ""), plot = p_ammonium_oxidation_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("anaerobic_ammonium_oxidation_heatmap.tiff", sep = ""), plot = p_anaerobic_ammonium_oxidation_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("nitrate_assimilatory_reduction_heatmap.tiff", sep = ""), plot = p_nitrate_assimilatory_reduction_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("nitrate_dissimilatory_reduction_heatmap.tiff", sep = ""), plot = p_nitrate_dissimilatory_reduction_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("nitrite_assimilation_heatmap.tiff", sep = ""), plot = p_nitrite_assimilation_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("nitrogen_fixation_heatmap.tiff", sep = ""), plot = p_nitrogen_fixation_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("urease_heatmap.tiff", sep = ""), plot = p_urease_heatmap,
       width=7, height=5, dpi=600)
ggsave(filename=paste("denitrification_heatmap.tiff", sep = ""), plot = p_denitrification_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrogen_heatmap_compartment.tiff", sep=""), plot=p_functions_heatmap_compartment,
       width=11, height=8.5, dpi=600) 
# Rhizosphere: Metadata Pairwise Correlation ----------------------------------------------
# Remove samples that lack metadata values 
FAMA_category_corr <- FAMA_category[-grep("Not Found", FAMA_category$pH),]

# Corvallis Site
FAMA_category_corr_melt <- melt(FAMA_category_corr)
FAMA_category_corr_melt$pH <- as.numeric(FAMA_category_corr_melt$pH)
FAMA_category_corr_melt$`Ca (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Ca (mg/kg (ppm))`)
FAMA_category_corr_melt$`K (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`K (mg/kg (ppm))`)
FAMA_category_corr_melt$`Mg (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Mg (mg/kg (ppm))`)
FAMA_category_corr_melt$`Mn (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Mn (mg/kg (ppm))`)
FAMA_category_corr_melt$`P (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`P (mg/kg (ppm))`)
FAMA_category_corr_melt$`Zn (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Zn (mg/kg (ppm))`)
FAMA_category_corr_melt$`NH4-N (mg/kg)` <- as.numeric(FAMA_category_corr_melt$`NH4-N (mg/kg)`)
FAMA_category_corr_melt$`NO3-N (mg/kg)` <- as.numeric(FAMA_category_corr_melt$`NO3-N (mg/kg)`)
FAMA_category_corr_melt$`C (%)` <- as.numeric(FAMA_category_corr_melt$`C (%)`)
FAMA_category_corr_melt$`N (%)` <- as.numeric(FAMA_category_corr_melt$`N (%)`)
FAMA_category_corr_melt$`MC (%)` <- as.numeric(FAMA_category_corr_melt$`MC (%)`)

FAMA_corvallis_corr_melt <- subset(FAMA_category_corr_melt, `Site Location` == "Corvallis")

# pH
ggscatter(FAMA_corvallis_corr_melt, x = "pH", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_pH

ggsave(filename=paste("corvallis_pH.tiff", sep=""), plot=corvallis_pH, 
       width=9.5, height=8, dpi=600)

# Ca
ggscatter(FAMA_corvallis_corr_melt, x = "Ca (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_Ca

ggsave(filename=paste("corvallis_Ca.tiff", sep=""), plot=corvallis_Ca, 
       width=9.5, height=8, dpi=600)

# K
ggscatter(FAMA_corvallis_corr_melt, x = "K (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_K

ggsave(filename=paste("corvallis_K.tiff", sep=""), plot=corvallis_K, 
       width=9.5, height=8, dpi=600)
# Mg
ggscatter(FAMA_corvallis_corr_melt, x = "Mg (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_Mg

ggsave(filename=paste("corvallis_Mg.tiff", sep=""), plot=corvallis_Mg, 
       width=9.5, height=8, dpi=600)
# Mn
ggscatter(FAMA_corvallis_corr_melt, x = "Mn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
         facet_wrap(~variable, scales = "free") -> corvallis_Mn

ggsave(filename=paste("corvallis_K.tiff", sep=""), plot=corvallis_Mn, 
       width=9.5, height=8, dpi=600)
# P
ggscatter(FAMA_corvallis_corr_melt, x = "P (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_P

ggsave(filename=paste("corvallis_P.tiff", sep=""), plot=corvallis_P, 
       width=9.5, height=8, dpi=600)
# Zn 
ggscatter(FAMA_corvallis_corr_melt, x = "Zn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_Zn

ggsave(filename=paste("corvallis_Zn.tiff", sep=""), plot=corvallis_Zn, 
       width=9.5, height=8, dpi=600)
# NH4-N
ggscatter(FAMA_corvallis_corr_melt, x = "NH4-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
         facet_wrap(~variable, scales = "free") -> corvallis_NH4_N

ggsave(filename=paste("corvallis_NH4-N.tiff", sep=""), plot=corvallis_NH4_N, 
       width=9.5, height=8, dpi=600)
# NO3-N
ggscatter(FAMA_corvallis_corr_melt, x = "NO3-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_NO3_N

ggsave(filename=paste("corvallis_NO3-N.tiff", sep=""), plot=corvallis_NO3_N, 
       width=9.5, height=8, dpi=600)
# C
ggscatter(FAMA_corvallis_corr_melt, x = "C (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_C

ggsave(filename=paste("corvallis_C.tiff", sep=""), plot=corvallis_C, 
       width=9.5, height=8, dpi=600)
# N
ggscatter(FAMA_corvallis_corr_melt, x = "N (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_N

ggsave(filename=paste("corvallis_N.tiff", sep=""), plot=corvallis_N, 
       width=9.5, height=8, dpi=600)
# MC
ggscatter(FAMA_corvallis_corr_melt, x = "MC (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          xlab("Soil Moisture Content (%)") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> corvallis_MC

ggsave(filename=paste("corvallis_MC.tiff", sep=""), plot=corvallis_MC, 
       width=9.5, height=8, dpi=600)

# Now do Clatskanie
FAMA_clatskanie_corr_melt <- subset(FAMA_category_corr_melt, `Site Location` == "Clatskanie")

# pH
ggscatter(FAMA_clatskanie_corr_melt, x = "pH", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_pH

ggsave(filename=paste("clatskanie_pH.tiff", sep=""), plot=clatskanie_pH, 
       width=9.5, height=8, dpi=600)

# Ca
ggscatter(FAMA_clatskanie_corr_melt, x = "Ca (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_Ca

ggsave(filename=paste("clatskanie_Ca.tiff", sep=""), plot=clatskanie_Ca, 
       width=9.5, height=8, dpi=600)

# K
ggscatter(FAMA_clatskanie_corr_melt, x = "K (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_K

ggsave(filename=paste("clatskanie_K.tiff", sep=""), plot=clatskanie_K, 
       width=9.5, height=8, dpi=600)
# Mg
ggscatter(FAMA_clatskanie_corr_melt, x = "Mg (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_Mg

ggsave(filename=paste("clatskanie_Mg.tiff", sep=""), plot=clatskanie_Mg, 
       width=9.5, height=8, dpi=600)
# Mn
ggscatter(FAMA_clatskanie_corr_melt, x = "Mn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_Mn

ggsave(filename=paste("clatskanie_K.tiff", sep=""), plot=clatskanie_Mn, 
       width=9.5, height=8, dpi=600)
# P
ggscatter(FAMA_clatskanie_corr_melt, x = "P (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_P

ggsave(filename=paste("clatskanie_P.tiff", sep=""), plot=clatskanie_P, 
       width=9.5, height=8, dpi=600)
# Zn 
ggscatter(FAMA_clatskanie_corr_melt, x = "Zn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_Zn

ggsave(filename=paste("clatskanie_Zn.tiff", sep=""), plot=clatskanie_Zn, 
       width=9.5, height=8, dpi=600)
# NH4-N
ggscatter(FAMA_clatskanie_corr_melt, x = "NH4-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_NH4_N

ggsave(filename=paste("clatskanie_NH4-N.tiff", sep=""), plot=clatskanie_NH4_N, 
       width=9.5, height=8, dpi=600)
# NO3-N
ggscatter(FAMA_clatskanie_corr_melt, x = "NO3-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_NO3_N

ggsave(filename=paste("clatskanie_NO3-N.tiff", sep=""), plot=clatskanie_NO3_N, 
       width=9.5, height=8, dpi=600)
# C
ggscatter(FAMA_clatskanie_corr_melt, x = "C (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_C

ggsave(filename=paste("clatskanie_C.tiff", sep=""), plot=clatskanie_C, 
       width=9.5, height=8, dpi=600)
# N
ggscatter(FAMA_clatskanie_corr_melt, x = "N (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          facet_wrap(~variable, scales = "free") -> clatskanie_N

ggsave(filename=paste("clatskanie_N.tiff", sep=""), plot=clatskanie_N, 
       width=9.5, height=8, dpi=600)
# MC
ggscatter(FAMA_corvallis_corr_melt, x = "MC (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
          ylab("Normalized Gene Abundance (EFPKG)") +
          xlab("Soil Moisture Content (%)") +
          facet_wrap(~variable, scales = "free") -> corvallis_MC

ggsave(filename=paste("corvallis_MC.tiff", sep=""), plot=corvallis_MC, 
       width=9.5, height=8, dpi=600)

# Rhizosphere: Statistics -------------------------------------------------------------------

#Subset FAMA_categories by nitrogen cycle process
FAMA_category_melt <- melt(FAMA_category) # melt FAMA category dataframe
FAMA_ammonification <- subset(FAMA_category_melt, variable == "Ammonification")
FAMA_ammonium_oxidation <- subset(FAMA_category_melt, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation <- subset(FAMA_category_melt, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification <- subset(FAMA_category_melt, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction <- subset(FAMA_category_melt, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction <- subset(FAMA_category_melt, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation <- subset(FAMA_category_melt, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation <- subset(FAMA_category_melt, variable == "Nitrogen fixation")
FAMA_urease <- subset(FAMA_category_melt, variable == "Urease")

FAMA_corvallis <- subset(FAMA_category_melt, `Site Location` == "Corvallis")
FAMA_ammonification_co <- subset(FAMA_corvallis, variable == "Ammonification")
FAMA_ammonium_oxidation_co <- subset(FAMA_corvallis, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation_co <- subset(FAMA_corvallis, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification_co <- subset(FAMA_corvallis, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction_co <- subset(FAMA_corvallis, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction_co <- subset(FAMA_corvallis, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation_co <- subset(FAMA_corvallis, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation_co <- subset(FAMA_corvallis, variable == "Nitrogen fixation")
FAMA_urease_co <- subset(FAMA_corvallis, variable == "Urease")

FAMA_clatskanie <- subset(FAMA_category_melt, `Site Location` == "Clatskanie")
FAMA_ammonification_cl <- subset(FAMA_clatskanie, variable == "Ammonification")
FAMA_ammonium_oxidation_cl <- subset(FAMA_clatskanie, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation_cl <- subset(FAMA_clatskanie, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification_cl <- subset(FAMA_clatskanie, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction_cl <- subset(FAMA_clatskanie, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction_cl <- subset(FAMA_clatskanie, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation_cl <- subset(FAMA_clatskanie, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation_cl <- subset(FAMA_clatskanie, variable == "Nitrogen fixation")
FAMA_urease_cl <- subset(FAMA_clatskanie, variable == "Urease")


# Perform two-way anova to determine if chemotype influences cycle gene abundance across both sites
anova_ammonification <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_ammonification)
anova_ammonium_oxidation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_ammonium_oxidation)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation)
anova_denitrification <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_denitrification)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction)
anova_nitrite_assimilation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrite_assimilation)
anova_nitrogen_fixation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrogen_fixation)
anova_urease <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_urease)
anova_nifB <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nifB)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction)
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation)
summary(anova_nitrogen_fixation)
summary(anova_urease)

# Perform one-way anova to determine if chemotype influnces cycle gene abundance for both sites independently
# Clatskanie
anova_ammonification <- aov(value ~ `Metabolite Level`, data = FAMA_ammonification_cl)
anova_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_ammonium_oxidation_cl)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation_cl)
anova_denitrification <- aov(value ~ `Metabolite Level`, data = FAMA_denitrification_cl)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction_cl)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction_cl)
anova_nitrite_assimilation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrite_assimilation_cl)
anova_nitrogen_fixation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrogen_fixation_cl)
anova_urease <- aov(value ~ `Metabolite Level`, data = FAMA_urease_cl)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction) # Clatskanie site is significant (p = 0.0453)
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation) #Clatskanie site is significant (p = 0.0394)
summary(anova_nitrogen_fixation)
summary(anova_urease)

# Identify which chemotypes are significantly different within Clatskanie dataset
anova_nitrate_assimilatory_reduction <- aov(value ~ `Chemotype`, data = FAMA_nitrate_assimilatory_reduction_cl)
anova_nitrite_assimilation <- aov(value ~ `Chemotype`, data = FAMA_nitrite_assimilation_cl)

summary(anova_nitrate_assimilatory_reduction)
summary(anova_nitrite_assimilation) 

# Corvallis
anova_ammonification <- aov(value ~ `Metabolite Level`, data = FAMA_ammonification_co)
anova_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_ammonium_oxidation_co)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation_co)
anova_denitrification <- aov(value ~ `Metabolite Level`, data = FAMA_denitrification_co)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction_co)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction_co)
anova_nitrite_assimilation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrite_assimilation_co)
anova_nitrogen_fixation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrogen_fixation_co)
anova_urease <- aov(value ~ `Metabolite Level`, data = FAMA_urease_co)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction) 
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation)
summary(anova_nitrogen_fixation)
summary(anova_urease)


anova_nitrogen_fixation <- aov(value ~ `Chemotype`, data = FAMA_nitrogen_fixation_co) # Check chemotype influence on nitrogen fixation
summary(anova_nitrogen_fixation)

# Rhizosphere: Box Plots --------------------------------------------------------------------

# Boxplot of all categories together
# Plot to compare each site
my_comparisons <- list(c("Clatskanie", "Corvallis"))
FAMA_category_melt %>%
  ggplot(aes(x=`Site Location`, y=value, fill = `Metabolite Level`)) +
  geom_boxplot(alpha=0.65) + 
  facet_wrap(~variable,
             scales = "free") +
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  stat_compare_means(comparisons = my_comparisons, method="t.test",
                     label.x = 1.5) +
  theme_bw() -> p_FAMA_category

ggsave(filename=paste("nitrogen_gene_boxplot_site.tiff", sep=""), plot=p_FAMA_category, 
       width=8.5, height=12, dpi=600)

# Plot each individual chemotype out
FAMA_category_melt %>%
  ggplot(aes(x=`Site Location`, y=value, fill = Chemotype)) +
  geom_boxplot(alpha=0.65) + 
  facet_wrap(~variable,
             scales = "free") +
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  theme_bw() -> p_FAMA_category

ggsave(filename=paste("nitrogen_gene_boxplot_chemotype.tiff", sep=""), plot=p_FAMA_category, 
       width=8.5, height=6, dpi=600)

# Subset boxplot and group by genotype within each fild site
FAMA_corvallis %>%
  ggplot(aes(x = `Chemotype`, y=value, fill = Genotype)) +
  geom_boxplot(alpha = 0.65) +
  facet_wrap(~variable,
             scales = 'free') + 
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  -> p_genotype_corvallis

ggsave(filename=paste("corvallis_genotype_boxplot.tiff", sep=""), plot = p_genotype_corvallis,
       width =10.5, height=8, dpi=600)

FAMA_clatskanie %>%
  ggplot(aes(x = `Chemotype`, y=value, fill = Genotype)) +
  geom_boxplot(alpha = 0.65) +
  facet_wrap(~variable,
             scales = 'free') + 
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  -> p_genotype_clatskanie

ggsave(filename=paste("clatskanie_genotype_boxplot.tiff", sep=""), plot = p_genotype_clatskanie,
       width =10.5, height=8, dpi=600)
# Subset boxplot for categories of interest
my_comparisons <- list( c("High", "Low"))
FAMA_clatskanie_plot <- subset(FAMA_clatskanie,
                               variable == "Nitrate assimilatory reduction" |
                               variable == "Nitrite assimilation")

FAMA_clatskanie_plot %>%
  ggplot(aes(x=`Metabolite Level`, y=value, fill = `Metabolite Level`)) +
  geom_boxplot(alpha=0.65, outlier.shape = NA) + 
  facet_wrap(~variable,
             scales = "free") +
  geom_jitter(shape = 16, alpha=0.75, size =2, position = position_jitter(0.08)) +
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +  
  ggtitle("Clatskanie Site Comparison") +
  stat_compare_means(comparisons = my_comparisons, method ="t.test") +
  theme_bw() -> p_FAMA_clatskanie_plot

ggsave(filename=paste("clatskanie_significant_plots.tiff", sep=""), plot = p_FAMA_clatskanie_plot,
       width =6, height=5, dpi=300)
# Rhizosphere: Heatmaps ---------------------------------------------------------------------
# Heatmap of all categories together

FAMA_functions_melt <- melt(FAMA_functions)
FAMA_functions_ammonification <- subset(FAMA_functions_melt, variable == "NrfA"|
                                        variable == "NrfB"| variable == "NrfC"|
                                        variable == "NrfD"| variable == "NrfH")
FAMA_functions_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "AmoA_PmoA"|
                                            variable == "AmoB_PmoB"| variable == "AmoC_PmoC")
FAMA_functions_anaerobic_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "HAO"|
                                                      variable == "Hzo"| variable =="HzsA"|
                                                      variable == "HzsB"| variable =="HzsC")
FAMA_functions_denitrification <- subset(FAMA_functions_melt, variable == "NirB3"|
                                         variable == "NirC"| variable == "NirK"|
                                         variable == "NirM"| variable == "NirN"|
                                         variable == "NirS"| variable == "NirT"|
                                         variable == "NosZ"| variable == "cNor-C"|
                                         variable == "cNorB_qNor")
FAMA_functions_nitrate_assimilatory_reduction <- subset(FAMA_functions_melt, variable == "NasA"|
                                                        variable == "NasB")
FAMA_functions_nitrate_dissimilatory_reduction <- subset(FAMA_functions_melt, variable == "NapA"|
                                                         variable == "NapB"| variable == "NapC"|
                                                         variable == "NapD"| variable == "NapE"|
                                                         variable == "NapF"| variable == "NapG"|
                                                         variable == "NapH"| variable == "NapK"|
                                                         variable == "NapL"| variable == "NarC"|
                                                         variable == "NarG_NxrA"| variable == "NarH_NxrB"|
                                                         variable == "NarI"| variable == "NarJ")
FAMA_functions_nitrite_assimilation <- subset(FAMA_functions_melt, variable == "NasI"|
                                              variable == "NasJ"| variable == "NirA"|
                                              variable == "NirB"| variable == "NirD"|
                                              variable == "NirU")
FAMA_functions_nitrogen_fixation <- subset(FAMA_functions_melt, variable == "AnfG_VnfG"|
                                           variable == "NifB"| variable == "NifD_AnfD_VnfD"|
                                           variable == "NifH_AnfH_VnfH"| variable == "NifK_AnfK_VnfK")
FAMA_functions_urease <- subset(FAMA_functions_melt, variable == "UreA"|
                                variable == "UreB"| variable == "UreC")

FAMA_functions_ammonification$category <- "Ammonification"
FAMA_functions_ammonium_oxidation$category <- "Ammonium oxidation"
FAMA_functions_anaerobic_ammonium_oxidation$category <- "Anaerobic ammonium oxidation"
FAMA_functions_denitrification$category <- "Denitrification"
FAMA_functions_nitrate_assimilatory_reduction$category <- "Nitrate assimilatory reduction"
FAMA_functions_nitrate_dissimilatory_reduction$category <- "Nitrate dissimilatory reduction"
FAMA_functions_nitrite_assimilation$category <- "Nitrite assimilation"
FAMA_functions_nitrogen_fixation$category <- "Nitrogen Fixation"
FAMA_functions_urease$category <- "Urease"

FAMA_functions_melt <- rbind(FAMA_functions_ammonification,
                             FAMA_functions_ammonium_oxidation,
                             FAMA_functions_anaerobic_ammonium_oxidation,
                             FAMA_functions_denitrification,
                             FAMA_functions_nitrate_assimilatory_reduction,
                             FAMA_functions_nitrate_dissimilatory_reduction,
                             FAMA_functions_nitrite_assimilation,
                             FAMA_functions_nitrogen_fixation,
                             FAMA_functions_urease)

FAMA_functions_corvallis <- subset(FAMA_functions_melt, `Site Location` == "Corvallis")

FAMA_functions_clatskanie <- subset(FAMA_functions_melt, `Site Location` == "Clatskanie")

# Total heatmap (Same Scale)
FAMA_functions_melt %>%
  ggplot(aes(`Site Location`, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_site

# Heatmap for each site
#By Chemotype
FAMA_functions_corvallis %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_corvallis

FAMA_functions_clatskanie %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_clatskanie

# Now do Genotype
FAMA_functions_corvallis %>%
  ggplot(aes(Genotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_corvallis_genotype

FAMA_functions_clatskanie %>%
  ggplot(aes(Genotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_clatskanie_genotype

# Subsetted heatmaps (Different scale for each category)
FAMA_functions_ammonification %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonification_heatmap

FAMA_functions_ammonium_oxidation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonium_oxidation_heatmap

FAMA_functions_anaerobic_ammonium_oxidation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Anaerobic ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_anaerobic_ammonium_oxidation_heatmap

FAMA_functions_denitrification %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Denitrification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_denitrification_heatmap

FAMA_functions_nitrate_assimilatory_reduction %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate assimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_assimilatory_reduction_heatmap

FAMA_functions_nitrate_dissimilatory_reduction %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate dissimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_dissimilatory_reduction_heatmap

FAMA_functions_nitrite_assimilation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrite assimilation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrite_assimilation_heatmap

FAMA_functions_nitrogen_fixation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrogen fixation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrogen_fixation_heatmap

FAMA_functions_urease %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Urease") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_urease_heatmap

ggsave(filename=paste("nitrogen_fixation_heatmap.tiff", sep = ""), plot = p_nitrogen_fixation_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrate_assimilatory_reduction_heatmap.tiff", sep=""), plot = p_nitrate_assimilatory_reduction_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrite_assimilation_heatmap.tiff", sep=""), plot = p_nitrite_assimilation_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("urease_heatmap.tiff", sep=""), plot = p_urease_heatmap,
       width=7, height=5, dpi=600)
                  
ggsave(filename=paste("nitrogen_heatmap_total.tiff", sep=""), plot=p_functions_heatmap_site,
       width=11, height=8.5, dpi=600)    

ggsave(filename=paste("nitrogen_heatmap_corvallis.tiff", sep=""), plot=p_functions_heatmap_corvallis,
       width=11, height=8.5, dpi=600)   

ggsave(filename=paste("nitrogen_genotype_heatmap_corvallis.tiff", sep=""), plot=p_functions_heatmap_corvallis_genotype,
       width=12.5, height=8.5, dpi=600)   

ggsave(filename=paste("nitrogen_genotype_heatmap_clatskanie.tiff", sep=""), plot=p_functions_heatmap_clatskanie_genotype,
       width=12.5, height=8.5, dpi=600)   


# Rhizosphere: PCA --------------------------------------------------------------------------
# Do a full PCA with all factors
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype", "pH",
                   "Ca (mg/kg (ppm))", "K (mg/kg (ppm))", "Mg (mg/kg (ppm))",
                   "Mn (mg/kg (ppm))", "P (mg/kg (ppm))", "Zn (mg/kg (ppm))", "NH4-N (mg/kg)",
                   "NO3-N (mg/kg)", "C (%)", "N (%)", "MC (%)")
                   
idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`,
                                           pH, `Ca (mg/kg (ppm))`, `K (mg/kg (ppm))`, `Mg (mg/kg (ppm))`,
                                           `Mn (mg/kg (ppm))`, `P (mg/kg (ppm))`, `Zn (mg/kg (ppm))`, `NH4-N (mg/kg)`,
                                           `NO2-N (mg/kg)`, `NO3-N (mg/kg)`, `C (%)`, `N (%)`, `MC (%)`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")
perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Site Location` + `Metabolite Level` + `Chemotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = `Metabolite Level`, fill = `Metabolite Level`, shape = `Site Location`), 
             size = 3, alpha = 0.65) +
  ggtitle("Nitrogen Gene Composition") +
  scale_shape_manual(values = c(21,24), name = "Site Location",
                     labels = c(expression("Clatskanie"), 
                                expression("Corvallis"))) +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
               #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA

ggsave(filename=paste("nitrogen_dbRDA_total.tiff", sep=""), plot=p_dbrda_FAMA, 
       width=5, height=4, dpi=400)   

# Now do just Corvallis 
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype")
idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- subset(FAMA_metadata, `Site Location` == "Corvallis")
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- subset(FAMA_table, `Site Location` == "Corvallis")
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")

perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Metabolite Level` + `Chemotype` + `Genotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = `Genotype`, fill = `Genotype`, shape = `Metabolite Level`), 
             size = 4.5, alpha = 0.65) +
  ggtitle("Corvallis Gene Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA_corvallis

# Now do just Clatskanie
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype")
idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- subset(FAMA_metadata, `Site Location` == "Clatskanie")
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- subset(FAMA_table, `Site Location` == "Clatskanie")
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")

perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Metabolite Level` + `Chemotype` + `Genotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = `Genotype`, fill = `Genotype`, shape = `Metabolite Level`), 
             size = 4.5, alpha = 0.65) +
  ggtitle("Clatskanie Gene Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA_clatskanie

leg <- get_legend(p_dbrda_FAMA_corvallis)
prow <- plot_grid(p_dbrda_FAMA_clatskanie + theme(legend.position = "none"), p_dbrda_FAMA_corvallis + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave(filename=paste("clatskanie_genotype_dbrda.tiff", sep=""), plot=p_dbrda_FAMA_clatskanie, 
       width=7.5, height=7.5, dpi=600) 

ggsave(filename=paste("corvallis_genotype_dbrda.tiff", sep=""), plot=p_dbrda_FAMA_corvallis,
       width=7.5, height=7.5, dpi=600)   





# Subset Data -- Soils only
FAMA_category_melt <- melt(FAMA_category) # melt FAMA category dataframe
FAMA_category_melt <- subset(FAMA_category_melt, `Compartment` == "Soil")

FAMA_category_corr <- FAMA_category[-grep("Not Found", FAMA_category$pH),]
FAMA_category_corr <- subset(FAMA_category_corr, `Compartment` == "Soil")

# Soil: Metadata Pairwise Correlation ----------------------------------------------
# Remove samples that lack metadata values 

# Corvallis Site
FAMA_category_corr_melt <- melt(FAMA_category_corr)
FAMA_category_corr_melt$pH <- as.numeric(FAMA_category_corr_melt$pH)
FAMA_category_corr_melt$`Ca (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Ca (mg/kg (ppm))`)
FAMA_category_corr_melt$`K (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`K (mg/kg (ppm))`)
FAMA_category_corr_melt$`Mg (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Mg (mg/kg (ppm))`)
FAMA_category_corr_melt$`Mn (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Mn (mg/kg (ppm))`)
FAMA_category_corr_melt$`P (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`P (mg/kg (ppm))`)
FAMA_category_corr_melt$`Zn (mg/kg (ppm))` <- as.numeric(FAMA_category_corr_melt$`Zn (mg/kg (ppm))`)
FAMA_category_corr_melt$`NH4-N (mg/kg)` <- as.numeric(FAMA_category_corr_melt$`NH4-N (mg/kg)`)
FAMA_category_corr_melt$`NO3-N (mg/kg)` <- as.numeric(FAMA_category_corr_melt$`NO3-N (mg/kg)`)
FAMA_category_corr_melt$`C (%)` <- as.numeric(FAMA_category_corr_melt$`C (%)`)
FAMA_category_corr_melt$`N (%)` <- as.numeric(FAMA_category_corr_melt$`N (%)`)
FAMA_category_corr_melt$`MC (%)` <- as.numeric(FAMA_category_corr_melt$`MC (%)`)

FAMA_corvallis_corr_melt <- subset(FAMA_category_corr_melt, `Site Location` == "Corvallis")

# pH
ggscatter(FAMA_corvallis_corr_melt, x = "pH", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_pH

ggsave(filename=paste("corvallis_pH.tiff", sep=""), plot=corvallis_pH, 
       width=9.5, height=8, dpi=600)

# Ca
ggscatter(FAMA_corvallis_corr_melt, x = "Ca (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_Ca

ggsave(filename=paste("corvallis_Ca.tiff", sep=""), plot=corvallis_Ca, 
       width=9.5, height=8, dpi=600)

# K
ggscatter(FAMA_corvallis_corr_melt, x = "K (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_K

ggsave(filename=paste("corvallis_K.tiff", sep=""), plot=corvallis_K, 
       width=9.5, height=8, dpi=600)
# Mg
ggscatter(FAMA_corvallis_corr_melt, x = "Mg (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_Mg

ggsave(filename=paste("corvallis_Mg.tiff", sep=""), plot=corvallis_Mg, 
       width=9.5, height=8, dpi=600)
# Mn
ggscatter(FAMA_corvallis_corr_melt, x = "Mn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_Mn

ggsave(filename=paste("corvallis_K.tiff", sep=""), plot=corvallis_Mn, 
       width=9.5, height=8, dpi=600)
# P
ggscatter(FAMA_corvallis_corr_melt, x = "P (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_P

ggsave(filename=paste("corvallis_P.tiff", sep=""), plot=corvallis_P, 
       width=9.5, height=8, dpi=600)
# Zn 
ggscatter(FAMA_corvallis_corr_melt, x = "Zn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_Zn

ggsave(filename=paste("corvallis_Zn.tiff", sep=""), plot=corvallis_Zn, 
       width=9.5, height=8, dpi=600)
# NH4-N
ggscatter(FAMA_corvallis_corr_melt, x = "NH4-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_NH4_N

ggsave(filename=paste("corvallis_NH4-N.tiff", sep=""), plot=corvallis_NH4_N, 
       width=9.5, height=8, dpi=600)
# NO3-N
ggscatter(FAMA_corvallis_corr_melt, x = "NO3-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_NO3_N

ggsave(filename=paste("corvallis_NO3-N.tiff", sep=""), plot=corvallis_NO3_N, 
       width=9.5, height=8, dpi=600)
# C
ggscatter(FAMA_corvallis_corr_melt, x = "C (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_C

ggsave(filename=paste("corvallis_C.tiff", sep=""), plot=corvallis_C, 
       width=9.5, height=8, dpi=600)
# N
ggscatter(FAMA_corvallis_corr_melt, x = "N (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_N

ggsave(filename=paste("corvallis_N.tiff", sep=""), plot=corvallis_N, 
       width=9.5, height=8, dpi=600)
# MC
ggscatter(FAMA_corvallis_corr_melt, x = "MC (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  xlab("Soil Moisture Content (%)") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> corvallis_MC

ggsave(filename=paste("corvallis_MC.tiff", sep=""), plot=corvallis_MC, 
       width=9.5, height=8, dpi=600)

# Now do Clatskanie
FAMA_clatskanie_corr_melt <- subset(FAMA_category_corr_melt, `Site Location` == "Clatskanie")

# pH
ggscatter(FAMA_clatskanie_corr_melt, x = "pH", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_pH

ggsave(filename=paste("clatskanie_pH.tiff", sep=""), plot=clatskanie_pH, 
       width=9.5, height=8, dpi=600)

# Ca
ggscatter(FAMA_clatskanie_corr_melt, x = "Ca (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_Ca

ggsave(filename=paste("clatskanie_Ca.tiff", sep=""), plot=clatskanie_Ca, 
       width=9.5, height=8, dpi=600)

# K
ggscatter(FAMA_clatskanie_corr_melt, x = "K (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_K

ggsave(filename=paste("clatskanie_K.tiff", sep=""), plot=clatskanie_K, 
       width=9.5, height=8, dpi=600)
# Mg
ggscatter(FAMA_clatskanie_corr_melt, x = "Mg (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_Mg

ggsave(filename=paste("clatskanie_Mg.tiff", sep=""), plot=clatskanie_Mg, 
       width=9.5, height=8, dpi=600)
# Mn
ggscatter(FAMA_clatskanie_corr_melt, x = "Mn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_Mn

ggsave(filename=paste("clatskanie_K.tiff", sep=""), plot=clatskanie_Mn, 
       width=9.5, height=8, dpi=600)
# P
ggscatter(FAMA_clatskanie_corr_melt, x = "P (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_P

ggsave(filename=paste("clatskanie_P.tiff", sep=""), plot=clatskanie_P, 
       width=9.5, height=8, dpi=600)
# Zn 
ggscatter(FAMA_clatskanie_corr_melt, x = "Zn (mg/kg (ppm))", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_Zn

ggsave(filename=paste("clatskanie_Zn.tiff", sep=""), plot=clatskanie_Zn, 
       width=9.5, height=8, dpi=600)
# NH4-N
ggscatter(FAMA_clatskanie_corr_melt, x = "NH4-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_NH4_N

ggsave(filename=paste("clatskanie_NH4-N.tiff", sep=""), plot=clatskanie_NH4_N, 
       width=9.5, height=8, dpi=600)
# NO3-N
ggscatter(FAMA_clatskanie_corr_melt, x = "NO3-N (mg/kg)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_NO3_N

ggsave(filename=paste("clatskanie_NO3-N.tiff", sep=""), plot=clatskanie_NO3_N, 
       width=9.5, height=8, dpi=600)
# C
ggscatter(FAMA_clatskanie_corr_melt, x = "C (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_C

ggsave(filename=paste("clatskanie_C.tiff", sep=""), plot=clatskanie_C, 
       width=9.5, height=8, dpi=600)
# N
ggscatter(FAMA_clatskanie_corr_melt, x = "N (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_N

ggsave(filename=paste("clatskanie_N.tiff", sep=""), plot=clatskanie_N, 
       width=9.5, height=8, dpi=600)
# MC
ggscatter(FAMA_corvallis_corr_melt, x = "MC (%)", y = "value",
          add = "reg.line", conf.int = TRUE,
          cor.coef = TRUE, cor.method = "spearman") +
  ylab("Normalized Gene Abundance (EFPKG)") +
  xlab("Soil Moisture Content (%)") +
  facet_wrap(~variable, scales = "free") -> clatskanie_MC

ggsave(filename=paste("clatskanie_MC.tiff", sep=""), plot=clatskanie_MC, 
       width=9.5, height=8, dpi=600)

# Soil: Statistics -------------------------------------------------------------------

#Subset FAMA_categories by nitrogen cycle process
FAMA_ammonification <- subset(FAMA_category_melt, variable == "Ammonification")
FAMA_ammonium_oxidation <- subset(FAMA_category_melt, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation <- subset(FAMA_category_melt, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification <- subset(FAMA_category_melt, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction <- subset(FAMA_category_melt, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction <- subset(FAMA_category_melt, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation <- subset(FAMA_category_melt, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation <- subset(FAMA_category_melt, variable == "Nitrogen fixation")
FAMA_urease <- subset(FAMA_category_melt, variable == "Urease")

FAMA_corvallis <- subset(FAMA_category_melt, `Site Location` == "Corvallis")
FAMA_ammonification_co <- subset(FAMA_corvallis, variable == "Ammonification")
FAMA_ammonium_oxidation_co <- subset(FAMA_corvallis, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation_co <- subset(FAMA_corvallis, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification_co <- subset(FAMA_corvallis, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction_co <- subset(FAMA_corvallis, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction_co <- subset(FAMA_corvallis, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation_co <- subset(FAMA_corvallis, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation_co <- subset(FAMA_corvallis, variable == "Nitrogen fixation")
FAMA_urease_co <- subset(FAMA_corvallis, variable == "Urease")

FAMA_clatskanie <- subset(FAMA_category_melt, `Site Location` == "Clatskanie")
FAMA_ammonification_cl <- subset(FAMA_clatskanie, variable == "Ammonification")
FAMA_ammonium_oxidation_cl <- subset(FAMA_clatskanie, variable == "Ammonium oxidation")
FAMA_anaerobic_ammonium_oxidation_cl <- subset(FAMA_clatskanie, variable == "Anaerobic ammonium oxidation")
FAMA_denitrification_cl <- subset(FAMA_clatskanie, variable == "Denitrification")
FAMA_nitrate_assimilatory_reduction_cl <- subset(FAMA_clatskanie, variable == "Nitrate assimilatory reduction")
FAMA_nitrate_dissimilatory_reduction_cl <- subset(FAMA_clatskanie, variable == "Nitrate dissimilatory reduction")
FAMA_nitrite_assimilation_cl <- subset(FAMA_clatskanie, variable == "Nitrite assimilation")
FAMA_nitrogen_fixation_cl <- subset(FAMA_clatskanie, variable == "Nitrogen fixation")
FAMA_urease_cl <- subset(FAMA_clatskanie, variable == "Urease")


# Perform two-way anova to determine if chemotype influences cycle gene abundance across both sites
anova_ammonification <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_ammonification)
anova_ammonium_oxidation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_ammonium_oxidation)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation)
anova_denitrification <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_denitrification)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction)
anova_nitrite_assimilation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrite_assimilation)
anova_nitrogen_fixation <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nitrogen_fixation)
anova_urease <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_urease)
anova_nifB <- aov(value ~ `Site Location` + `Metabolite Level`, data = FAMA_nifB)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction)
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation)
summary(anova_nitrogen_fixation)
summary(anova_urease)

# Perform one-way anova to determine if chemotype influnces cycle gene abundance for both sites independently
# Clatskanie
anova_ammonification <- aov(value ~ `Metabolite Level`, data = FAMA_ammonification_cl)
anova_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_ammonium_oxidation_cl)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation_cl)
anova_denitrification <- aov(value ~ `Metabolite Level`, data = FAMA_denitrification_cl)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction_cl)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction_cl)
anova_nitrite_assimilation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrite_assimilation_cl)
anova_nitrogen_fixation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrogen_fixation_cl)
anova_urease <- aov(value ~ `Metabolite Level`, data = FAMA_urease_cl)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction) # Clatskanie site is significant (p = 0.0453)
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation) #Clatskanie site is significant (p = 0.0394)
summary(anova_nitrogen_fixation)
summary(anova_urease)

# Corvallis
anova_ammonification <- aov(value ~ `Metabolite Level`, data = FAMA_ammonification_co)
anova_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_ammonium_oxidation_co)
anova_anaerobic_ammonium_oxidation <- aov(value ~ `Metabolite Level`, data = FAMA_anaerobic_ammonium_oxidation_co)
anova_denitrification <- aov(value ~ `Metabolite Level`, data = FAMA_denitrification_co)
anova_nitrate_assimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_assimilatory_reduction_co)
anova_nitrate_dissimilatory_reduction <- aov(value ~ `Metabolite Level`, data = FAMA_nitrate_dissimilatory_reduction_co)
anova_nitrite_assimilation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrite_assimilation_co)
anova_nitrogen_fixation <- aov(value ~ `Metabolite Level`, data = FAMA_nitrogen_fixation_co)
anova_urease <- aov(value ~ `Metabolite Level`, data = FAMA_urease_co)

summary(anova_ammonification)
summary(anova_ammonium_oxidation)
summary(anova_anaerobic_ammonium_oxidation)
summary(anova_denitrification)
summary(anova_nitrate_assimilatory_reduction) 
summary(anova_nitrate_dissimilatory_reduction)
summary(anova_nitrite_assimilation)
summary(anova_nitrogen_fixation)
summary(anova_urease)




# Soil: Box Plots --------------------------------------------------------------------

# Boxplot of all categories together
# Plot to compare each site
my_comparisons <- list(c("Clatskanie", "Corvallis"))
FAMA_category_melt %>%
  ggplot(aes(x=`Site Location`, y=value, fill = `Metabolite Level`)) +
  geom_boxplot(alpha=0.65) + 
  facet_wrap(~variable,
             scales = "free") +
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  stat_compare_means(comparisons = my_comparisons, method="t.test",
                     label.x = 1.5) +
  theme_bw() -> p_FAMA_category

ggsave(filename=paste("nitrogen_gene_boxplot_site.tiff", sep=""), plot=p_FAMA_category, 
       width=8.5, height=12, dpi=600)

# Plot each individual chemotype out
FAMA_category_melt %>%
  ggplot(aes(x=`Site Location`, y=value, fill = Chemotype)) +
  geom_boxplot(alpha=0.65) + 
  facet_wrap(~variable,
             scales = "free") +
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5), linetype = "dotted") +
  theme_bw() -> p_FAMA_category

ggsave(filename=paste("nitrogen_gene_boxplot_chemotype.tiff", sep=""), plot=p_FAMA_category, 
       width=8.5, height=6, dpi=600)

# Subset boxplot and group by genotype within each fild site
FAMA_corvallis %>%
  ggplot(aes(x = `Chemotype`, y=value, fill = Genotype)) +
  geom_boxplot(alpha = 0.65) +
  facet_wrap(~variable,
             scales = 'free') + 
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  -> p_genotype_corvallis

ggsave(filename=paste("corvallis_genotype_boxplot.tiff", sep=""), plot = p_genotype_corvallis,
       width =10.5, height=8, dpi=600)

FAMA_clatskanie %>%
  ggplot(aes(x = `Chemotype`, y=value, fill = Genotype)) +
  geom_boxplot(alpha = 0.65) +
  facet_wrap(~variable,
             scales = 'free') + 
  scale_y_continuous(name = "Normalized Gene Abundance (EFPKG)") +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), linetype = "dotted") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))  -> p_genotype_clatskanie

ggsave(filename=paste("clatskanie_genotype_boxplot.tiff", sep=""), plot = p_genotype_clatskanie,
       width =10.5, height=8, dpi=600)

# Soil: Heatmaps ---------------------------------------------------------------------
# Heatmap of all categories together

FAMA_functions_melt <- melt(FAMA_functions)
FAMA_functions_ammonification <- subset(FAMA_functions_melt, variable == "NrfA"|
                                          variable == "NrfB"| variable == "NrfC"|
                                          variable == "NrfD"| variable == "NrfH")
FAMA_functions_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "AmoA_PmoA"|
                                              variable == "AmoB_PmoB"| variable == "AmoC_PmoC")
FAMA_functions_anaerobic_ammonium_oxidation <- subset(FAMA_functions_melt, variable == "HAO"|
                                                        variable == "Hzo"| variable =="HzsA"|
                                                        variable == "HzsB"| variable =="HzsC")
FAMA_functions_denitrification <- subset(FAMA_functions_melt, variable == "NirB3"|
                                           variable == "NirC"| variable == "NirK"|
                                           variable == "NirM"| variable == "NirN"|
                                           variable == "NirS"| variable == "NirT"|
                                           variable == "NosZ"| variable == "cNor-C"|
                                           variable == "cNorB_qNor")
FAMA_functions_nitrate_assimilatory_reduction <- subset(FAMA_functions_melt, variable == "NasA"|
                                                          variable == "NasB")
FAMA_functions_nitrate_dissimilatory_reduction <- subset(FAMA_functions_melt, variable == "NapA"|
                                                           variable == "NapB"| variable == "NapC"|
                                                           variable == "NapD"| variable == "NapE"|
                                                           variable == "NapF"| variable == "NapG"|
                                                           variable == "NapH"| variable == "NapK"|
                                                           variable == "NapL"| variable == "NarC"|
                                                           variable == "NarG_NxrA"| variable == "NarH_NxrB"|
                                                           variable == "NarI"| variable == "NarJ")
FAMA_functions_nitrite_assimilation <- subset(FAMA_functions_melt, variable == "NasI"|
                                                variable == "NasJ"| variable == "NirA"|
                                                variable == "NirB"| variable == "NirD"|
                                                variable == "NirU")
FAMA_functions_nitrogen_fixation <- subset(FAMA_functions_melt, variable == "AnfG_VnfG"|
                                             variable == "NifB"| variable == "NifD_AnfD_VnfD"|
                                             variable == "NifH_AnfH_VnfH"| variable == "NifK_AnfK_VnfK")
FAMA_functions_urease <- subset(FAMA_functions_melt, variable == "UreA"|
                                  variable == "UreB"| variable == "UreC")

FAMA_functions_ammonification$category <- "Ammonification"
FAMA_functions_ammonium_oxidation$category <- "Ammonium oxidation"
FAMA_functions_anaerobic_ammonium_oxidation$category <- "Anaerobic ammonium oxidation"
FAMA_functions_denitrification$category <- "Denitrification"
FAMA_functions_nitrate_assimilatory_reduction$category <- "Nitrate assimilatory reduction"
FAMA_functions_nitrate_dissimilatory_reduction$category <- "Nitrate dissimilatory reduction"
FAMA_functions_nitrite_assimilation$category <- "Nitrite assimilation"
FAMA_functions_nitrogen_fixation$category <- "Nitrogen Fixation"
FAMA_functions_urease$category <- "Urease"

FAMA_functions_melt <- rbind(FAMA_functions_ammonification,
                             FAMA_functions_ammonium_oxidation,
                             FAMA_functions_anaerobic_ammonium_oxidation,
                             FAMA_functions_denitrification,
                             FAMA_functions_nitrate_assimilatory_reduction,
                             FAMA_functions_nitrate_dissimilatory_reduction,
                             FAMA_functions_nitrite_assimilation,
                             FAMA_functions_nitrogen_fixation,
                             FAMA_functions_urease)

FAMA_functions_corvallis <- subset(FAMA_functions_melt, `Site Location` == "Corvallis")

FAMA_functions_clatskanie <- subset(FAMA_functions_melt, `Site Location` == "Clatskanie")

# Total heatmap (Same Scale)
FAMA_functions_melt %>%
  ggplot(aes(`Site Location`, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_site

# Heatmap for each site
#By Chemotype
FAMA_functions_corvallis %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_corvallis

FAMA_functions_clatskanie %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_clatskanie

# Now do Genotype
FAMA_functions_corvallis %>%
  ggplot(aes(Genotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_corvallis_genotype

FAMA_functions_clatskanie %>%
  ggplot(aes(Genotype, variable, fill = value)) +
  facet_wrap(~category,
             scales = "free") +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_functions_heatmap_clatskanie_genotype

# Subsetted heatmaps (Different scale for each category)
FAMA_functions_ammonification %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonification_heatmap

FAMA_functions_ammonium_oxidation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_ammonium_oxidation_heatmap

FAMA_functions_anaerobic_ammonium_oxidation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Anaerobic ammonium oxidation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_anaerobic_ammonium_oxidation_heatmap

FAMA_functions_denitrification %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Denitrification") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_denitrification_heatmap

FAMA_functions_nitrate_assimilatory_reduction %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate assimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_assimilatory_reduction_heatmap

FAMA_functions_nitrate_dissimilatory_reduction %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrate dissimilatory reduction") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrate_dissimilatory_reduction_heatmap

FAMA_functions_nitrite_assimilation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrite assimilation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrite_assimilation_heatmap

FAMA_functions_nitrogen_fixation %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Nitrogen fixation") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_nitrogen_fixation_heatmap

FAMA_functions_urease %>%
  ggplot(aes(Chemotype, variable, fill = value)) +
  geom_tile() +
  labs(fill = "Normalized Gene Abundance (EFPKG)") +
  ylab("Gene Name") +
  labs(title = "Urease") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_fill_viridis(discrete=FALSE) -> p_urease_heatmap

ggsave(filename=paste("nitrogen_fixation_heatmap.tiff", sep = ""), plot = p_nitrogen_fixation_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrate_assimilatory_reduction_heatmap.tiff", sep=""), plot = p_nitrate_assimilatory_reduction_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrite_assimilation_heatmap.tiff", sep=""), plot = p_nitrite_assimilation_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("urease_heatmap.tiff", sep=""), plot = p_urease_heatmap,
       width=7, height=5, dpi=600)

ggsave(filename=paste("nitrogen_heatmap_total.tiff", sep=""), plot=p_functions_heatmap_site,
       width=11, height=8.5, dpi=600)    

ggsave(filename=paste("nitrogen_heatmap_corvallis.tiff", sep=""), plot=p_functions_heatmap_corvallis,
       width=11, height=8.5, dpi=600)   

ggsave(filename=paste("nitrogen_heatmap_clatskanie.tiff", sep=""), plot=p_functions_heatmap_clatskanie,
       width=11, height=8.5, dpi=600) 

ggsave(filename=paste("nitrogen_genotype_heatmap_corvallis.tiff", sep=""), plot=p_functions_heatmap_corvallis_genotype,
       width=12.5, height=8.5, dpi=600)   

ggsave(filename=paste("nitrogen_genotype_heatmap_clatskanie.tiff", sep=""), plot=p_functions_heatmap_clatskanie_genotype,
       width=12.5, height=8.5, dpi=600)   


# Soil: PCA --------------------------------------------------------------------------
# Do a full PCA with all factors
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype", "pH",
                   "Ca (mg/kg (ppm))", "K (mg/kg (ppm))", "Mg (mg/kg (ppm))",
                   "Mn (mg/kg (ppm))", "P (mg/kg (ppm))", "Zn (mg/kg (ppm))", "NH4-N (mg/kg)",
                   "NO3-N (mg/kg)", "C (%)", "N (%)", "MC (%)", "Compartment")

idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`,
                                           pH, `Ca (mg/kg (ppm))`, `K (mg/kg (ppm))`, `Mg (mg/kg (ppm))`,
                                           `Mn (mg/kg (ppm))`, `P (mg/kg (ppm))`, `Zn (mg/kg (ppm))`, `NH4-N (mg/kg)`,
                                           `NO2-N (mg/kg)`, `NO3-N (mg/kg)`, `C (%)`, `N (%)`, `MC (%)`, `Compartment`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")
perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Site Location` + `Metabolite Level` + `Chemotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = `Metabolite Level`, fill = `Metabolite Level`, shape = `Site Location`), 
             size = 3, alpha = 0.65) +
  ggtitle("Nitrogen Gene Composition") +
  scale_shape_manual(values = c(21,24), name = "Site Location",
                     labels = c(expression("Clatskanie"), 
                                expression("Corvallis"))) +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA

ggsave(filename=paste("nitrogen_dbRDA_total.tiff", sep=""), plot=p_dbrda_FAMA, 
       width=5, height=4, dpi=400)   

# Now do just Corvallis 
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype", "pH",
                   "Ca (mg/kg (ppm))", "K (mg/kg (ppm))", "Mg (mg/kg (ppm))",
                   "Mn (mg/kg (ppm))", "P (mg/kg (ppm))", "Zn (mg/kg (ppm))", "NH4-N (mg/kg)",
                   "NO3-N (mg/kg)", "C (%)", "N (%)", "MC (%)", "Compartment")
idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- subset(FAMA_metadata, `Site Location` == "Corvallis")
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- subset(FAMA_table, `Site Location` == "Corvallis")
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`,
                                           pH, `Ca (mg/kg (ppm))`, `K (mg/kg (ppm))`, `Mg (mg/kg (ppm))`,
                                           `Mn (mg/kg (ppm))`, `P (mg/kg (ppm))`, `Zn (mg/kg (ppm))`, `NH4-N (mg/kg)`,
                                           `NO2-N (mg/kg)`, `NO3-N (mg/kg)`, `C (%)`, `N (%)`, `MC (%)`, `Compartment`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")

perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Metabolite Level` + `Chemotype` + `Genotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = Chemotype, fill = Chemotype, shape = `Metabolite Level`), 
             size = 4.5, alpha = 0.65) +
  ggtitle("Corvallis Gene Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA_corvallis

# Now do just Clatskanie
Metadata_list <- c("Sample_ID", "Chemotype", "Metabolite Level", "Site Location", "Genotype", "pH",
                   "Ca (mg/kg (ppm))", "K (mg/kg (ppm))", "Mg (mg/kg (ppm))",
                   "Mn (mg/kg (ppm))", "P (mg/kg (ppm))", "Zn (mg/kg (ppm))", "NH4-N (mg/kg)",
                   "NO3-N (mg/kg)", "C (%)", "N (%)", "MC (%)", "Compartment")
idx <- match(Metadata_list, names(FAMA_functions))
idx <- sort(c(idx))
FAMA_metadata <- FAMA_functions[,idx]
FAMA_metadata <- subset(FAMA_metadata, `Site Location` == "Clatskanie")
FAMA_metadata <- FAMA_metadata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")

FAMA_table <- FAMA_functions
FAMA_table <- subset(FAMA_table, `Site Location` == "Clatskanie")
FAMA_table <- FAMA_table %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
FAMA_table <- dplyr::select(FAMA_table, -c(Chemotype, `Metabolite Level`, `Site Location`, `Genotype`,
                                           pH, `Ca (mg/kg (ppm))`, `K (mg/kg (ppm))`, `Mg (mg/kg (ppm))`,
                                           `Mn (mg/kg (ppm))`, `P (mg/kg (ppm))`, `Zn (mg/kg (ppm))`, `NH4-N (mg/kg)`,
                                           `NO2-N (mg/kg)`, `NO3-N (mg/kg)`, `C (%)`, `N (%)`, `MC (%)`, `Compartment`))
FAMA_table <- as.matrix(FAMA_table)
FAMA_table <- as.data.frame(FAMA_table)

dist_FAMA <- vegdist(FAMA_table, "jaccard")

perm = 9999
set.seed(309125)
rda_FAMA <- dbrda(dist_FAMA ~ `Metabolite Level` + `Chemotype` + `Genotype`, FAMA_metadata, permutations = perm) 
anova(rda_FAMA) 
print(anova(rda_FAMA, by="terms", permutation = 9999))

plot.data <- merge(summary(rda_FAMA)$sites, FAMA_metadata, by = "row.names") 
ggplot() +
  geom_point(data = plot.data, aes(x = dbRDA1, y = dbRDA2, colour = Chemotype, fill = Chemotype, shape = `Metabolite Level`), 
             size = 4.5, alpha = 0.65) +
  ggtitle("Clatskanie Gene Composition") +
  scale_x_continuous(name = paste0(paste0("dbRDA 1 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,1], 1)), "%)")) +
  scale_y_continuous(name = paste0(paste0("dbRDA 2 (", round(100*as.data.frame(summary(rda_FAMA)$cont)[2,2], 1)), "%)")) +
  #geom_segment(data = as.data.frame(summary(rda_FAMA)$biplot), aes(x = 0, y = 0, xend = dbRDA1, yend = dbRDA2),
  #size = 0.25, color = "black", arrow = arrow()) +
  #geom_text(data = as.data.frame(summary(rda_FAMA)$biplot) %>% rownames_to_column(), aes(x = dbRDA1, y = dbRDA2, label = rowname)) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), axis.text = element_text(size = 10),
        plot.title = element_text(size=12, face="bold"),
        legend.text = element_text(hjust = 0), title = element_text(size = 12),
        legend.title = element_text(size = 12)) -> p_dbrda_FAMA_clatskanie

leg <- get_legend(p_dbrda_FAMA_corvallis)
prow <- plot_grid(p_dbrda_FAMA_clatskanie + theme(legend.position = "none"), p_dbrda_FAMA_corvallis + theme(legend.position = "none"),
                  ncol = 2, labels = "AUTO", align = "vh")
prow <- plot_grid(prow, leg, ncol = 2, rel_widths = c(0.8, 0.2))

ggsave(filename=paste("clatskanie_genotype_dbrda.tiff", sep=""), plot=p_dbrda_FAMA_clatskanie, 
       width=5, height=7.5, dpi=600) 

ggsave(filename=paste("corvallis_genotype_dbrda.tiff", sep=""), plot=p_dbrda_FAMA_corvallis,
       width=5, height=7.5, dpi=600)   


