#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                           Control vs FGF Treatments                          #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          FGF_data.R
# Outputs: Figure 2E
#          Figure 2G
#          Supplemental Figure 2D
#          Supplemental Figure 2F
#          Supplemental Figure 2G
#          Figure 3H
#          Supplemental Figure 3G
#..............................................................................#


#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('FGF_data.R')


#------------------------------------------------------------------------------# 
# Plots ----

# Colour palette
palette <- c ( Control = "#F8766D", FGF4 = "#C77CFF", ICM = "purple", 
  DP = "purple", EPI = "red", PrE = "blue", "Polar TE" = "#80FF80", 
  "Mural TE" = "#00A500", "TE" = "green")

# Appearance for boxplots
box_theme <- theme (aspect.ratio = 1,  
  axis.line = element_line (size = 0.5), 
  axis.text.x = element_text (size = 7), 
  axis.text.y = element_text (size = 7),
  legend.position = "none", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8))

#------------------------------------------------------------------------------# 
# Supplemental Figure 2D ----
# Heatmap of ERK-KTR C:N values over time, by embryo and treamtment

# Duplicate parent tracks so they are double height of daughters
# This allows parent track to "split" into two daughters in the heatmap
# and keeping the alignments in-line

parent <- filter (embryos_CN_FGF, Mitosis == "Parent")
parent$MtUniqueID <- paste0 (parent$UniqueID, 1)

parent2 <- filter (embryos_CN_FGF, Mitosis == "Parent")
parent2$MtUniqueID <- paste0 (parent2$UniqueID, 2)

d1 <- filter (embryos_CN_FGF, Mitosis == "Daughter1")
d1$MtUniqueID <- paste0 (d1$UniqueID, 1)

d2 <- filter (embryos_CN_FGF, Mitosis == "Daughter2")
d2$MtUniqueID <- paste0 (d2$UniqueID, 2)

duplicate <- rbind(parent, parent2, d1, d2)

# Plot heatmap
g <- duplicate %>% 
  filter (TE.ICM == "ICM" & !grepl ("Mitosis", Label)) %>%
  ggplot (aes (x = TimeM, y = reorder (MtUniqueID, trackAvg)))
  # Plot individual ICM cells, during interphase, over time, 
  # ordered by mean ERK-KTR C:N of each cell
g <- g + geom_raster (aes (fill = Value), hjust = 0.5, vjust = 0.5, 
  interpolate = FALSE) + scale_x_continuous (expand = c (0, 0))
  # Heatmap with ERK-KTR C:N values
g <- g + viridis::scale_fill_viridis (option="A", 
  rescaler = function (x, to = c (0, 1), from = NULL) {
  ifelse(x < 3, scales::rescale (x, to = to, from = c (min (x, na.rm = TRUE), 
  3)), 1)
  })
  # Set scale for heatmap
g <- g + facet_wrap (Treatment ~ Cell_number ~ Embryo.Id, scales = "free", 
  ncol = 5)
  # Divide and order plots by embryo, cell number and treatment group
g <- g + theme_classic () + 
  theme(aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7, angle = 90), 
  axis.text.y = element_blank (), 
  axis.ticks.y = element_blank (), 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  legend.position = "bottom", 
  panel.background = element_rect (fill = "grey20", color  =  NA),
  panel.spacing.x = unit (.05, "lines"),
  panel.spacing.y = unit (.05, "lines"),
  strip.background = element_blank (),
  strip.text = element_text (size = 6),
  text = element_text (size = 8)) +
  labs (x = 'Time (min)', y = 'ERK-KTR C:N')
  # Customise plot appearance
print (g)

ggsave ("SuppFig2B.tiff", width = 3) # Save TIFF
ggsave ("SuppFig2B.pdf", width = 3) # Save PDF

#------------------------------------------------------------------------------# 
# Figure 2E ----
# Line plot of ERK-KTR C:N values over time per cell, by treatment

# Plot line plot
g <- embryos_CN_FGF %>% 
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track &
  !grepl ("Mitosis", Label)) %>% 
  ggplot (aes (x = TimeM, y = Value))
  # Plot ERK-KTR C:N values for ICM cells tracked for > 5 frames, 
  # during interphase, over time
g <- g + geom_line (aes (group = MtUniqueID), color = "black", alpha = 0.1, 
  size = 0.5) 
  # Plot line plot for individual cells
g <- g + geom_smooth (aes (color = factor (Treatment, 
  levels = c ("Control", "FGF4"))), 
  size=1, fill = "white", se = TRUE, method = loess) + 
  guides (color = guide_legend (title = "Population Mean"))
  # Overlay population mean for each treatment
g <- g + facet_wrap (Treatment ~ ., nrow = 1)
  # Divide plots by treatment
g <- g + theme_classic() + 
  theme(axis.line = element_line(size = 0.5),
  axis.text = element_text(size = 7), 
  legend.position="none", 
  legend.text = element_text(size = 6), 
  legend.title = element_text(size = 6),  
  strip.background = element_rect(size = 0.5),
  text = element_text(size=8)) + 
  coord_fixed(25) +
  scale_color_manual (values = palette) +
  scale_x_continuous (expand = c (0, 0)) + 
  labs (x = 'Time (min)', y = 'ERK-KTR C:N')
  # Customise plot appearance
print (g)

ggsave ("Fig2E.tiff", width = 2.3) # Save TIFF


#------------------------------------------------------------------------------# 
# Figure 2G ----
# Boxplot of mean ERK-KTR C:N per cell, by treatment

# Plot boxplot
g <- avg %>%
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = factor (Treatment, levels = c ("Control", "FGF4")), 
  y = trackAvg))
  # Plot mean ERK-KTR C:N values for ICM cells tracked for > 5 frames
g <- g + geom_boxplot (aes (fill = factor (Treatment, 
  levels = c ("Control", "FGF4"))), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
  # Plot boxplots by treatment
g <- g + geom_jitter (aes (color = factor (Treatment, 
  levels = c ("Control", "FGF4"))), position = position_jitter (width = 0.2), 
  alpha = 0.5,  size = 1, stroke = 0) 
  # Overlay with points for mean ERK-KTR C:N values in individual cells
g <- g + theme_classic () + 
  box_theme + 
  coord_cartesian (ylim = c (0.5, 3)) + 
  labs (y = "Mean per cell", x = "") + 
  guides (fill = guide_legend (title = ""), color = guide_legend (title = "")) +
  scale_color_manual (values = palette) +  
  scale_fill_manual (values = palette)
  # Customise plot appearance
print (g)

ggsave("Fig2F.tiff", width = 1.25) # Save TIFF


# Wilcoxon test pairwise comparison

icm <- avg %>%
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track)
# For ICM cells tracked > 5 frames

emb_avg <- icm %>%
  group_by (Embryo.Id, Treatment) %>% 
  summarise (emb_mean = mean (trackAvg))
# Calculate mean ERK-KTR C:N per embryo

wilcox.test (emb_mean ~ Treatment, data = emb_avg) 
# p-value = 0.0003199


#------------------------------------------------------------------------------# 
# Supplemental Figure 2F ----
# Boxplot of coefficient of variance in ERK-KTR C:N per cell, by treatment

# Plot boxplot
g <- avg %>%
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = factor (Treatment, levels = c ("Control", "FGF4")), 
  y = trackCV))
  # Plot CV in ERK-KTR C:N values for ICM cells tracked for > 5 frames
g <- g + geom_boxplot (aes (fill = factor (Treatment, 
  levels = c("Control", "FGF4"))), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
  # Plot boxplots by treatment
g <- g + geom_jitter (aes (color = factor (Treatment, 
  levels = c("Control", "FGF4"))), position = position_jitter (width = 0.2), 
  alpha = 0.5, size = 1, stroke = 0) 
  # Overlay with points for CV ERK-KTR C:N values in individual cells
g <- g + theme_classic () +
  box_theme + 
  labs (y = "CV per cell", x = "") + 
  guides (fill = guide_legend (title = ""), color = guide_legend (title = "")) +
  scale_color_manual (values = palette) +  
  scale_fill_manual (values = palette)
  # Customise plot appearance
print (g)

ggsave ("SuppFig2F.tiff", width = 1.25) # Save TIFF


#------------------------------------------------------------------------------# 
# Figure 3H ----

# Plot boxplot
g <- avg %>%
  filter (TE.ICM == "TE" & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = Identity, y = log (trackAvg)))
  # Plot log mean ERK-KTR C:N values for ICM cells tracked for > 5 frames
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
  # Plot boxplots by treatment
g <- g + geom_jitter (aes (color = Identity), 
  position = position_jitter (width = 0.2), alpha = 0.5,  size = 1, stroke = 0) 
  # Overlay with points for mean ERK-KTR C:N values in individual cells
g <- g + theme_classic () + 
  labs (y = "Log(Mean ERK activity)", x = "") +
  box_theme + 
  facet_grid ( ~ Treatment) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette)
  # Customise plot appearance
print (g)

ggsave ("Fig3H.tiff", width = 2) # Save as TIFF

# Wilcoxon non-parametric test pairwise comparison

subset_TE <- subset (avg, TE.ICM == "TE" & MtUniqueID %in% long_track)
# For TE cells tracked for > 5 frames

TE_avg <- subset_TE %>%
  group_by (Embryo.Id, Treatment, Identity) %>%
  summarise (emb_mean = mean (trackAvg))
# Calculate mean ERK-KTR C:N per embryo

# Test within lineage, effect of treatment
wilcox.test (emb_mean ~ Treatment, data = subset (TE_avg, Identity == "Polar TE")) 
# n.s.
wilcox.test (emb_mean ~ Treatment, data = subset (TE_avg, Identity == "Mural TE")) 
# p = 4.571e-05

# Test within treatment, effect of Identity
wilcox.test (emb_mean ~ Identity, data = subset (TE_avg, Treatment == "FGF4")) 
# n.s.
wilcox.test (emb_mean ~ Identity, data = subset (TE_avg, Treatment == "Control")) 
# p = 0.0004871


#------------------------------------------------------------------------------# 
# Supplemental Figure 3G ----
# Scatter plot of mean distance of TE to nearest ICM neighbor against mean 
# ERK-KTR C:N value per cell

# Calculate TE nearest ICM neighbor for each frame

# Group data into TE and ICM
TE <- subset (embryos_CN_FGF, TE.ICM == "TE")
ICM <- subset (embryos_CN_FGF, TE.ICM == "ICM")

# Get ID info for individul TE cells at each time point
SpotID <- unique (TE$SpotID)

# Create empty table for distance measurements
dist <- TE
dist$d <- 0
dist <- dist[NULL, ]

# Store X , Y and Z coordinates for subsetting data
coord <- c ("X", "Y", "Z")

# For each TE cell get all ICM coordinates for that embryo at each time point
for (i in 1: length (SpotID)) {
  spoti <- TE[TE$SpotID == SpotID[i], ]
  Embryo.Idi <- spoti$Embryo.Id
  Timei <- spoti$Time
  ICMi <- subset (ICM, Embryo.Id == Embryo.Idi & Time == Timei)
  ICM.loc <- ICMi[coord] 
  TE.loc <- spoti[coord]

# Caclulate the distance between the TE cell and ICM cells
  a <- NULL
  for (j in 1: length (ICM.loc)) {
    sub <- ICM.loc[j, ]
    d.X <- sub$X - TE.loc$X # Caclulate distance in each dimension
    d.Y <- sub$Y - TE.loc$Y
    d.Z <- sub$Z - TE.loc$Z
    d <- sqrt (d.X ^ 2 + d.Y ^ 2 + d.Z ^ 2) # Calculate distance in 3D
    a <- rbind(a, d) 
  }
  
# Find the distance between the TE cell and its nearest neighbor ICM
  nn <- min (a) 
  
# Aggregate data on distance betweeen TE cell and its nearest ICM neighbor
  spoti$d <- nn
  dist <- rbind (dist, spoti)
}

# Summary stats on TE-ICM nearest neighbor distances
summary_TE <- dist %>% 
  group_by (MtUniqueID, Identity, Embryo.Id, trackAvg, Cell_number, Stage, 
  Treatment) %>% 
  summarise (d = mean (d))
summary_TE <- subset (summary_TE, MtUniqueID %in% long_track)

# Plot scatter plot
g <- summary_TE %>% 
  ggplot(aes (x = d, y = log (trackAvg)))
  # Plot individual TE cells mean distance from nearest ICM neighbor and
  # log mean ERK-KTR C:N values
g <- g + geom_point (aes (color = Treatment), size = 0.5, alpha = 0.5)
  # Colour code by treatment
g <- g + theme_classic()  +
  theme (aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "none",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) + 
  facet_wrap ( ~ Treatment) +
  labs (y = "Log[mean ERK activity]", x = "Mean distance from ICM (μM)") +  
  scale_color_manual (values = palette)
  # Customise plot appearance
print (g)

ggsave ("SuppFig3G.tiff", width = 2) # Save TIFF


# Correlation of TE-ICM distance and ERK activity in TE cells

# Pearson coefficient, control embryos
control <- subset (summary_TE, Treatment == "Control")
cor (control$d, log (control$trackAvg), method = c ("pearson")) # r = -0.6470163

# Pearson coefficient, FGF treated embryos
fgf <- subset (summary_TE, Treatment == "FGF4")
cor (fgf$d, log (fgf$trackAvg), method = c ("pearson")) # r = -0.4405173