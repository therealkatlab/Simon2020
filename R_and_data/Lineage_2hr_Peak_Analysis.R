#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                      Peak call in 2hr time-lapse movies                      #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          Lineage_2hr_Data.R
#          smoothpks.csv output from Matlab script PeakCall_smooth.m
#          smoothdat.csv output from Matlab script PeakCall_smooth.m
# Outputs: Figure 5B
#          Figure 5E
#          Figure 5F
#          Figure 5G
#          Supplemental Figure 5C
#          Figure 5H
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

## Load packages, themes and data from other file ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('Lineage_2hr_Data.R')

palette <- c (DP = "purple", EPI = "red", PrE = "blue", "Polar TE" = "#80FF80", 
  "Mural TE" = "#00A500", ICM = "purple3")

write_csv(good_spots, "embryos_CN.csv")

#------------------------------------------------------------------------------# 
# Run Matlab FindPeaks Script ---- PeakCall_smooth.m
#------------------------------------------------------------------------------# 

#------------------------------------------------------------------------------# 
# Read Matlab FindPeaks Output ----
#------------------------------------------------------------------------------# 

pksmatlab <- read.csv("smoothpks.csv")
embryos_CN <- read.csv("smoothdat.csv") 
# Note Matlab CSV which removes . so TE.ICM now TEICM

#------------------------------------------------------------------------------# 
# Run Matlab Threshold  Script ----
#------------------------------------------------------------------------------# 

pksmatlab <- pksmatlab %>% filter (pks > 1.2142) # threshold peaks fron script

# Data wrangling
 
colnames (pksmatlab)[colnames (pksmatlab) == "ID"] <- "MtUniqueID"

Identity <- embryos_CN %>% group_by(MtUniqueID, Cell_number, Stage, Identity, 
  EmbryoId) %>% summarise(freq=n())

min <- embryos_CN %>% group_by(MtUniqueID) %>% summarise(min=min(smoothed), 
  med = median (smoothed))

pks <- dplyr::full_join (pksmatlab, Identity, by = "MtUniqueID")

pks <- dplyr::left_join (pks, min, by = "MtUniqueID")

pks$h <- pks$pk - pks$min

pks$hmed <- pks$pk - pks$med

pks$wmin <- pks$w * 5

pks_y <- pks %>% 
  filter (!is.na (pks)) %>% 
  group_by (EmbryoId, Cell_number, Stage, MtUniqueID, Identity, freq) %>% 
  summarise (numPk = n ())

pks_n <- pks %>% 
  filter (is.na (pks)) %>% 
  group_by (EmbryoId, Cell_number, Stage, MtUniqueID, Identity, freq) %>% 
  summarise(numPk = 0)

summary <- rbind (pks_y, pks_n)

summary$freqPk <- summary$numPk / (summary$freq * 5 / 60) 

## Counts for numPk ####

counts <- summary %>% 
  group_by (Identity, numPk) %>% 
  summarise (freq = n ())

counts2 <- summary %>% 
  group_by (Identity, numPk, Stage) %>% 
  summarise (freq = n ())

pks$bsl <- pks$pks - pks$p

# grouping for stats test

cell_group <- pks %>% 
  group_by (Stage, MtUniqueID, EmbryoId, Identity) %>% 
  summarise (p = mean (p), w = mean (wmin), max = mean (pks), bsl = mean (bsl))

cell_group <- subset (cell_group, !is.na (p))

emb_group <- cell_group %>% 
  filter (Stage != "<32"  & 
            (Identity == "PrE" | Identity == "DP" | Identity == "EPI")) %>% 
  group_by (Stage, EmbryoId, Identity) %>% 
  summarise (mean_p = mean (p), mean_w = mean (w), mean_max = mean (max), 
             mean_bsl = mean (bsl))




#------------------------------------------------------------------------------# 
# Figure 5B ----

sample <- c ("082119_p11000000008Parent",  "082119_p11000000002Parent", 
            "082119_p11000000000Parent")

g <- embryos_CN %>% 
  filter (MtUniqueID %in% sample) %>% 
  ggplot (aes (x = TimeM, y = smoothed))
g <- g + geom_line (aes (color = Identity), size = 0.5)
g <- g + geom_point (data = subset (pks, MtUniqueID %in% sample), 
  aes(x = TimeM, y = pks), size = 0.5, color = "black")
g <- g + facet_wrap (Identity ~ MtUniqueID)
g <- g + theme_classic () + 
  theme (aspect.ratio = 0.75,
  axis.text.x = element_text (size = 6, angle = 90),
  axis.text.y = element_text (size = 6),
  legend.position = "top", 
  legend.text = element_text (size = 7), 
  legend.title = element_blank (), 
  panel.border = element_blank (),
  panel.grid.major = element_blank (), 
  panel.grid.minor = element_blank (), 
  panel.spacing.x = unit (.5, "lines"), 
  panel.spacing.y = unit (.5, "lines"), 
  strip.background = element_blank (), 
  strip.text = element_blank (), 
  text = element_text (size = 7)) +
  labs (y = 'ERK-KTR C:N', x = "Time (min)") + 
  scale_color_manual (values = palette) + 
  scale_y_continuous (breaks = c (1,2,3))
print (g)
ggsave ("Fig5B.tiff", width = 3)

#------------------------------------------------------------------------------# 
# Figure 5E ----
# Peak Maxima

g <- pks %>%
  filter (Stage != "<32" & (Identity == "DP" | Identity == "EPI" | Identity == "PrE")) %>%
  ggplot (aes (x = Identity, y = pks, color = Identity)) 
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.3, 
                      size = 1, stroke = 0) 
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
                       color = "black", width = 0.5)
g <- g + theme_classic ()
g <- g + labs (y = "Peak maximum", x = "") +
  theme (aspect.ratio = 1, 
         axis.line = element_line (size = 0.5),
         axis.text.x = element_text (size = 7, angle = 90), 
         axis.text.y = element_text (size = 6),
         legend.position = "none", 
         legend.text = element_text (size = 7), 
         legend.title = element_text (size = 7),  
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 7)) + 
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) 
print(g)
ggsave("Fig5E.tiff", width = 1.5)

pairwise.wilcox.test(emb_group$mean_max, emb_group$Identity,
                     p.adjust.method = "bonferroni") # n.s


#------------------------------------------------------------------------------# 
# Figure 5F ----
# Peak duration (width)

g <- pks %>% 
  filter ( (Identity == "DP" | Identity == "EPI" | Identity == "PrE") 
  & Stage != "<32") %>% 
  ggplot (aes (x = Identity, y = wmin, color = Identity))
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.3, 
  size = 1, stroke = 0)
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", width = 0.5)
g <- g + theme_classic ()
g <- g + labs (y = "Peak duration (min)", x="") +
  theme(axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 7, angle = 90), 
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7),  
        legend.position="none", 
        axis.line = element_line(size = 0.5),
        strip.background = element_rect(size = 0.5),
        text = element_text(size=7),
        aspect.ratio = 1) +  scale_color_manual(values = palette) + 
  scale_fill_manual(values = palette) 
print(g)
ggsave("Fig5F.tiff", width = 1.5)

pairwise.wilcox.test(emb_group$mean_w, emb_group$Identity,
                     p.adjust.method = "bonferroni") # n.s.


#------------------------------------------------------------------------------# 
# Figure 5G ----
# Peak prominence

g <- pks %>% 
  filter (Stage != "<32" & 
  (Identity == "DP" | Identity == "EPI" | Identity == "PrE")) %>% 
  ggplot (aes (x = Identity, y = log(p), color = Identity))
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.3, 
  size = 1, stroke = 0) 
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", width = 0.5)
g <- g + theme_classic ()
g <- g + labs (y = "Log(Peak prominence)", x = "") +
  theme (aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 7, angle = 90), 
  axis.text.y = element_text (size = 6),
  legend.position = "none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette)
print(g)
ggsave("Fig5G.tiff", width = 1.5)

pairwise.wilcox.test (emb_group$mean_p, emb_group$Identity,
                      p.adjust.method = "bonferroni") # EPI PrE p = 0.016


#------------------------------------------------------------------------------# 
# Supplemental Figure 5C ----
# Nearest minimum

g <- pks %>% 
  filter (Stage != "<32" & 
  (Identity == "DP" | Identity == "EPI" | Identity == "PrE")) %>% 
  ggplot (aes (x = Identity, y = bsl, color = Identity))
g <- g + geom_jitter (position = position_jitter (width = 0.2), alpha = 0.3, 
  size = 1, stroke = 0)
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", width = 0.5)
g <- g + theme_classic ()
g <- g + labs (y = "Peak nearest minimum", x = "") +
  theme(aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 7, angle = 90), 
  axis.text.y = element_text (size = 6),
  legend.position = "none", 
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette)
print (g)
ggsave ("SuppFig5C.tiff", width = 1.5)

pairwise.wilcox.test(emb_group$mean_bsl, emb_group$Identity,
                     p.adjust.method = "bonferroni") # EPI PrE n.s p = 0.088


#------------------------------------------------------------------------------# 
# Figure 5H ----
# Histogram of # peaks / cell

g <- summary %>% 
  filter ( (Identity == "DP" | Identity == "EPI" | Identity == "PrE") & 
  Stage != "<32") %>%
  ggplot (aes (x = numPk))  
g <- g + geom_histogram (aes (y = stat (width * density), fill = Identity, 
  color = Identity), binwidth = 1, alpha = 0.5)
g <- g + theme_classic ()
g <- g + theme(aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 6), 
  axis.text.y = element_text (size = 6),
  legend.text = element_text (size = 7), 
  legend.title = element_text (size = 7),  
  legend.position = "none", 
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 7)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  scale_y_continuous (labels = scales::percent) +
  facet_grid (. ~ Identity) + 
  labs (x = "# of peaks / 2hr", y = "% cells")
print (g)
ggsave ("Fig5H.tiff", width = 2)