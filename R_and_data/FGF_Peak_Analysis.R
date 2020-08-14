#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                  Peak call in Control vs FGF Treatments                      #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          FGF_data.R
#          smoothpks_FGF.csv output from Matlab script PeakCall_smooth_FGF.m
#          smoothdat_FGF.csv output from Matlab script PeakCall_smooth_FGF.m
# Outputs: Supplemental Figure 5B
#          Supplemental Figure 5D
#          Supplemental Figure 5E
#          Supplemental Figure 5F
#          Supplemental Figure 5G
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('FGF_data.R')

palette <- c(Control = "#F8766D", FGF4 = "#C77CFF")

write_csv(good_spots, "FGF_peaks.csv")

#------------------------------------------------------------------------------# 
# Run Matlab FindPeaks Script ---- PeakCall_smooth_FGF.m
#------------------------------------------------------------------------------# 

#------------------------------------------------------------------------------# 
# Read Matlab FindPeaks Output ---- 
#------------------------------------------------------------------------------# 

pksmatlab <- read.csv("smoothpks_FGF.csv")
embryos_CN_FGF <- read.csv("smoothdat_FGF.csv") 

# Note Matlab CSV which removes . so TE.ICM now TEICM

#------------------------------------------------------------------------------# 
# Run Matlab Threshold  Script ----
#------------------------------------------------------------------------------# 

pksmatlab <- pksmatlab %>% 
  filter(pks > 1.2142) # threshold peaks as being > 1

colnames(pksmatlab)[colnames(pksmatlab) == "ID"] <- "MtUniqueID"

Treatment <- embryos_CN_FGF %>% 
  group_by (MtUniqueID, TEICM, Treatment, EmbryoId) %>% 
  summarise (freq = n()) 

min <- embryos_CN_FGF %>% 
  group_by (MtUniqueID) %>% 
  summarise(min = min (smoothed), med = median (smoothed))

pks <- dplyr::full_join(pksmatlab, Treatment, by = "MtUniqueID")

pks <- dplyr::left_join(pks, min, by = "MtUniqueID")

pks$h <- pks$pk - pks$min

pks$hmed <- pks$pk - pks$med

pks$wmin <- pks$w * 5

pks_y <- pks %>% 
  filter(!is.na(pks)) %>% 
  group_by (EmbryoId, TEICM,  MtUniqueID, Treatment, freq) %>% 
  summarise(numPk = n())

pks_n <- pks %>% 
  filter(is.na(pks)) %>% 
  group_by(EmbryoId, TEICM, MtUniqueID, Treatment, freq) %>% 
  summarise(numPk = 0)

summary <- rbind(pks_y, pks_n)

summary$freqPk <- summary$numPk / (summary$freq*5/60) 

counts <- summary %>% group_by(Treatment, TEICM, numPk) %>% summarise(freq=n())

pks$bsl <- pks$pks - pks$p


# Grouping for stats test
cell_group <- pks %>% 
  group_by (MtUniqueID, TEICM, EmbryoId, Treatment) %>% 
  summarise( p = mean(p), w = mean(wmin), max = mean(pks), bsl = mean(bsl))
cell_group <- subset (cell_group, !is.na(p))
emb_group <- cell_group %>% 
  group_by (EmbryoId, TEICM, Treatment) %>% 
  summarise (mean_p = mean(p), mean_w = mean(w), mean_max = mean(max), 
            mean_bsl = mean(bsl))
ICM_group <- subset (emb_group, TEICM == "ICM")

#------------------------------------------------------------------------------# 
# Supplemental Figure 5B ----
g <- summary %>% 
  filter (TEICM == "ICM") %>% 
  ggplot (aes (x = numPk)) 
g <- g + geom_histogram (aes (y = stat(width * density), fill = Treatment, 
  color = Treatment), binwidth = 1, alpha = 0.5)
g <- g + theme_classic ()
g <- g + theme(aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 6), 
  axis.text.y = element_text (size = 6),
  legend.position="none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) + 
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) + 
  scale_y_continuous(labels=scales::percent) +
  facet_grid(.~Treatment) + 
  labs (x = '# of peaks / 2hr', y = "% cells")
print (g)
ggsave ("SuppFig5B.tiff", width = 2)


#------------------------------------------------------------------------------# 
# Supplemental Figure 5D ----

g <- pks %>% 
  filter(TEICM == "ICM") %>%
  ggplot (aes (x = Treatment, y = wmin, color = Treatment)) 
g <- g + geom_boxplot (aes (fill = Treatment), alpha = 0.5, outlier.size = 0, 
  color = "black", width = 0.5)
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.5, 
  size = 1, stroke = 0) 
g <- g + theme_classic() + 
  theme (aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 7, angle = 90), 
  axis.text.y = element_text (size = 6),
  legend.position = "none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) + labs (y = 'Peak duration (min)', x = "") +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) 
print(g)
ggsave("SuppFig5D.tiff", width = 1.5)

wilcox.test(ICM_group$mean_w ~ ICM_group$Treatment) # p = 0.743 n.s.


#------------------------------------------------------------------------------# 
# Supplemental Figure 5E ----

g <- pks %>% 
  filter(TEICM == "ICM") %>%
  ggplot(aes(x = Treatment, y = log(p), color = Treatment))
g <- g + geom_boxplot (aes (fill = Treatment), alpha = 0.5, outlier.size = 0, 
  color="black", width = 0.5)
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.5, 
  size = 1, stroke = 0) 
g <- g + theme_classic() +
  theme(aspect.ratio = 1,
  axis.line = element_line(size = 0.5),
  axis.text.x = element_text (size = 7, angle = 90), 
  axis.text.y = element_text (size = 6),
  legend.position = "none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) + 
  labs(y = "Log(Peak prominence)", x = "") +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette)
print(g)
ggsave("SuppFig5E.tiff", width = 1.5)

wilcox.test(ICM_group$mean_p ~ ICM_group$Treatment) # p = 0.0274 * Peak prominence


#------------------------------------------------------------------------------# 
# Supplemental Figure 5F ----

g <- pks %>% 
  filter (TEICM == "ICM") %>%
  ggplot (aes (x = Treatment, y = pks, color = Treatment)) 
g <- g + geom_boxplot (aes (fill = Treatment), alpha = 0.5, outlier.size = 0, 
  color = "black", width = 0.5)
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.5, 
  size = 1, stroke = 0) 
g <- g + theme_classic() + 
  theme (aspect.ratio = 1,
  axis.line = element_line(size = 0.5),
  axis.text.x = element_text(size = 7, angle = 90), 
  axis.text.y = element_text(size = 6),
  legend.text = element_text(size = 7), 
  legend.title = element_text(size = 7),  
  legend.position = "none", 
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) + 
  labs(y = "Peak maximum", x = "") +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) 
print(g)
ggsave("SuppFig5F.tiff", width = 1.5)

wilcox.test(ICM_group$mean_max ~ ICM_group$Treatment) # p = 0.03595 * Peak Maxima


#------------------------------------------------------------------------------# 
# Supplemental Figure 5G ----

g <- pks %>%
  filter (TEICM == "ICM") %>%
  ggplot (aes (x = Treatment, y = pks, color = Treatment))
g <- g + geom_boxplot (aes (fill = Treatment), alpha = 0.5, outlier.size = 0, 
                       color = "black", width = 0.5)
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.3, 
  size = 1, stroke = 0) 
g <- g + theme_classic() +
  theme(aspect.ratio = 1, 
        axis.line = element_line (size = 0.5),
        axis.text.x = element_text (size = 7, angle = 90), 
        axis.text.y = element_text (size = 6),
        legend.position = "none", 
        legend.text = element_text (size = 7), 
        legend.title = element_text (size = 7),  
        strip.background = element_rect (size = 0.5),
        text = element_text (size = 7)) + 
  labs(y = "Peak nearest minimum", x = "") +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette)
print(g)
ggsave("SuppFig5G.tiff", width = 1.5)

wilcox.test(ICM_group$mean_bsl ~ ICM_group$Treatment) # p = 0.0009872 * Peak bsl
