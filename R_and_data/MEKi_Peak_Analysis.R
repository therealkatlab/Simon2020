#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                 Peak call in Control vs MEKi Treatments                      #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          MEKi_data.R
#          smoothpks_PD03.csv output from Matlab script PeakCall_smooth_PD03.m
#          smoothdat_PD03.csv output from Matlab script PeakCall_smooth_PD03.m
# Outputs: Supplemental Figure 5A
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('MEKi_Data.R')



## SAVE
write_csv(good_spots, "PD03_peaks.csv")


#------------------------------------------------------------------------------# 
# Run Matlab FindPeaks Script ---- PeakCall_smooth_PD03.m
#------------------------------------------------------------------------------# 

#------------------------------------------------------------------------------# 
# Read Matlab FindPeaks Output ----
#------------------------------------------------------------------------------# 

pksmatlab <- read.csv("smoothpks_PD03.csv")
embryos_CN_PD03 <- read.csv("smoothdat_PD03.csv") 
# Note Matlab CSV which removes . so TE.ICM now TEICM

#------------------------------------------------------------------------------# 
# Run Matlab Threshold  Script ----
#------------------------------------------------------------------------------# 

pksmatlab <- pksmatlab %>% 
  filter (pks > 1.2142)
# Threshold from Matlab threshold script

colnames(pksmatlab)[colnames(pksmatlab)=="ID"] <- "MtUniqueID"

Treatment <- embryos_CN_PD03 %>% 
  group_by (MtUniqueID, TEICM, Treatment, EmbryoId) %>% 
  summarise(freq=n()) 

min <- embryos_CN_PD03 %>% 
  group_by (MtUniqueID) %>% 
  summarise (min = min(smoothed), med = median(smoothed))

pks <- dplyr::full_join(pksmatlab, Treatment, by = "MtUniqueID")

pks <- dplyr::left_join(pks, min, by = "MtUniqueID")

pks$h <- pks$pk - pks$min

pks$hmed <- pks$pk - pks$med

pks$wmin <- pks$w * 5

pks_y <- pks %>% 
  filter(!is.na(pks)) %>% 
  group_by(EmbryoId, TEICM,  MtUniqueID, Treatment, freq) %>% 
  summarise(numPk = n())

pks_n <- pks %>% 
  filter(is.na(pks)) %>% 
  group_by(EmbryoId, TEICM, MtUniqueID, Treatment, freq) %>% 
  summarise(numPk = 0)

summary <- rbind(pks_y, pks_n)

summary$freqPk <- summary$numPk / (summary$freq * 5 / 60) 

## Counts for numPk ####

counts <- summary %>% group_by(Treatment, TEICM, numPk) %>% summarise(freq=n())


#------------------------------------------------------------------------------# 
# Supplemental Figure 5A ----

g <- summary %>% 
  filter (TEICM == "ICM") %>%
  ggplot (aes (x = numPk)) 
g <- g + geom_histogram (aes (y = stat (width * density), 
  fill = fct_rev (Treatment), color = fct_rev (Treatment)), binwidth = 1, 
  alpha = 0.5)
g <- g + theme_classic() + 
  theme (aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 6), 
  axis.text.y = element_text (size = 6),
  legend.position = "none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) + 
  scale_y_continuous (labels=scales::percent) +
  facet_grid ( . ~fct_rev(Treatment)) + 
  labs (x = '# of peaks / 2hr', y = "% cells")
print (g)
ggsave ("SuppFig5A.tiff", width = 2)