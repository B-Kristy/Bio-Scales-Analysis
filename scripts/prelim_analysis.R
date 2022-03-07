# Load libraries and set working directory 
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(viridis)
setwd("~/Dropbox/BioScales/prelim_results/")

# Load and process input data 
rhizo_denitrification <- read.table("Rhizosphere_denitrification.txt", header = T, sep = "\t")
all_na <- function(x) any(!is.na(x)) # remove empty columns
rhizo_denitrification <- rhizo_denitrification %>% select_if(all_na)
rhizo_denitrification$cycle <- "Denitrification" # add cycle column

rhizo_nitrification <- read.table("Rhizosphere_nitrification.txt", header =T, sep = "\t")
rhizo_nitrification$cycle <- "Nitrification" # add cycle column

rhizo_nitrogen_fixation <- read.table("rhizosphere_nitrogen_fixation.txt", header = T, sep = "\t")
rhizo_nitrogen_fixation <- rhizo_nitrogen_fixation %>% select_if(all_na) # remove empty columns
rhizo_nitrogen_fixation$cycle <- "Nitrogen Fixation"

rhizo_nitrogen_cycle <- rbind(rhizo_denitrification, rhizo_nitrification, rhizo_nitrogen_fixation)
write.table(rhizo_nitrogen_cycle, "rhizosphere_nitrogen_cycle.txt", sep = "\t", row.names =T, col.names=T) # table is exported and reformatted in excel (matrix)

rhizo_kegg <- read.table("Rhizosphere_KEGG.txt", header = T, sep = "\t")
# Visualize function set data (Heat map)
rhizo_nitrogen_cycle <- read.table ("Rhizosphere_nitrogen_function_set.txt", header = T, sep = "\t")

rhizo_denitrification <- rhizo_nitrogen_cycle %>% filter((Cycle == "Denitrification")) # split into separate heat maps
rhizo_nitrification <- rhizo_nitrogen_cycle %>% filter((Cycle == "Nitrification"))
rhizo_nitrogen_fiation <- rhizo_nitrogen_cycle %>% filter((Cycle == "Nitrogen Fixation"))

p1 <- ggplot(rhizo_denitrification, aes(Chemotype, Function.Name, fill = Abundance )) +
        geom_tile() +
        labs(title = "Denitrification Genes", x = "", y ="") +
        theme(axis.text.x = element_blank()) +
        scale_fill_viridis(discrete=FALSE) 

p2 <- ggplot(rhizo_nitrification, aes(Chemotype, Function.Name, fill = Abundance )) +
        geom_tile() +
        labs(title = "Nitrification Genes", x = "", y ="") +
        theme(axis.text.x = element_blank()) +
        scale_fill_viridis(discrete=FALSE) 

p3 <- ggplot(rhizo_nitrogen_fiation, aes(Chemotype, Function.Name, fill = Abundance )) +
        geom_tile() +
        labs(title = "Nitrogen Fixation Genes", x = "Chemotype", y ="") +
        theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
        scale_fill_viridis(discrete=FALSE) 

prow_heatmaps <- ggarrange(p1, p2, p3, align = "hv", ncol =1, labels = c("A", "B", "C"))

ggsave(filename=paste("heatmap.tiff",sep=""), plot=prow_heatmaps, width=9, height=9, dpi=600)


# Visualize boxplot of significant KEGG differential abundance pathway (boxplot)
rhizo_nitrate_reduction <- rhizo_kegg %>% filter ((Feature == "M00530"))



# Visualize taxonomic data (phylogenetic tree; heatmap of significantly abundant microbes)



