#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                        Control vs MEKi Treatments                           #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          MEKi_Data.R
# Outputs: Figure 2D
#          Figure 2E
#          Supplemental Figure 2C
#          Supplemental Figure 2E
#..............................................................................#


#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('MEKi_Data.R')

#------------------------------------------------------------------------------# 
# Plots ----

#------------------------------------------------------------------------------# 
# Figure 2D ----
# Line plot of ERK-KTR C:N values over time per cell, by treatment

# Plot line plot
g <- embryos_CN_PD03 %>% 
  filter (TE.ICM == "ICM" & !grepl ("Mitosis", Label)) %>%
  ggplot (aes (x = TimeM, y = Value))
  # Plot ERK-KTR C:N values for ICM cells tracked for > 5 frames, 
  # during interphase, over time
g <- g + geom_line (aes (group = MtUniqueID), color = "black", alpha = 0.1, 
  size = 0.5) 
  # Plot line plot for individual cells
g <- g + geom_smooth (aes (color = factor (Treatment, 
  levels = c ("Control", "1uM PD03"))), size = 1, fill = "white", 
  method = loess, se=TRUE) + guides (color = guide_legend (title = "Population Mean"))
  # Overlay with population mean for each treatment
g <- g + facet_wrap (fct_rev (Treatment) ~ ., nrow = 1)
  # Divide plots by Treatment
g <- g + theme_classic() + scale_x_continuous (expand = c (0, 0)) + 
  coord_fixed (25) +
  theme ( axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7), 
  legend.position = "none", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +
  labs (x = "Time (min)", y = "ERK-KTR C:N") 
  # Modify plot appearance
print (g)
ggsave ("Fig2D.tiff", width = 2.3) # Save as TIFF


#------------------------------------------------------------------------------# 
# Figure 2E ----
# Boxplot of mean ERK-KTR C:N values per cell by treatement: Control vs MEKi

g <- avg %>% 
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = factor(Treatment, levels = c ("Control", "1uM PD03")), 
  y = trackAvg))
  # Plot mean ERK-KTR C:N values for each cell, by treatment
g <- g + geom_boxplot (aes (fill = factor (Treatment, 
  levels = c ("Control", "1uM PD03"))), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
  # Plot boxplot of mean ERK-KTR C:N values for each cell color coded by treatment
g <- g + geom_jitter (aes (color = factor (Treatment, 
  levels = c ("Control", "1uM PD03"))), 
  position = position_jitter (width = 0.2), alpha=0.5,  size = 1, stroke = 0) 
  # Overlay plot of individual data points of mean ERK-KTR C:N values for each cell
g <- g + theme_classic() + coord_cartesian (ylim = c (0.5, 3))
  # Set axis limits to directly compare with FGF treatment plots
  theme(aspect.ratio = 1,
  axis.line = element_line(size = 0.5),
  axis.text.x = element_text (size = 7), 
  axis.text.y = element_text (size = 7),
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  legend.position = "none", 
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +
  labs(y = 'Mean per cell', x="") + 
  scale_fill_discrete(name = "Treatment", labels = c("Control", "MEKi"))
  # Modfiy Plot appearance
print(g)
ggsave("Fig2C.tiff", width = 1.25) # Save as TIFF

# Wilcoxon test
subset <- subset (avg, TE.ICM == "ICM" & MtUniqueID %in% long_track)
  # Plot only ICM cells 
  # Include cells tracked for > 5 frames
emb_avg <- subset %>% 
  group_by (Embryo.Id, Treatment) %>% 
  summarise (emb_mean = mean (trackAvg))
  # Calculate mean for each embryo
wilcox.test (emb_mean ~ Treatment, data = emb_avg)
  # Wilcoxon non-parametric test p-value = 0.01587

#------------------------------------------------------------------------------# 
# Supplemental Figure 2C ----
# Heatmap of ERK-KTR C:N values over time, by embryo and treamtment

# Duplicate parent tracks so they are double height of daughters
# This allows parent track to "split" into two daughters in the heatmap
# and keeping the alignments in-line
parent <- filter (embryos_CN_PD03, Mitosis == "Parent") 
parent$MtUniqueID <- paste0 (parent$UniqueID, 1)
  # Parent cell duplicate 1
parent2 <- filter (embryos_CN_PD03, Mitosis == "Parent")
parent2$MtUniqueID <- paste0 (parent2$UniqueID, 2)
  # Parent cell duplicate 2
d1 <- filter (embryos_CN_PD03, Mitosis == "Daughter1") 
d1$MtUniqueID <- paste0 (d1$UniqueID, 1)
  # Daughter cell 1
d2 <- filter (embryos_CN_PD03, Mitosis == "Daughter2") 
d2$MtUniqueID <- paste0 (d2$UniqueID, 2)
  # Daughter cell 2
duplicate <- rbind (parent, parent2, d1, d2) 
  # Bind together for plotting heatmap

# Plot heatmap
g <- duplicate %>% 
  filter (TE.ICM == "ICM" & !grepl ("Mitosis", Label)) %>%
  ggplot (aes (x = TimeM, y = reorder (MtUniqueID, trackAvg)))
  # Plot individual ICM cells, during interphase, over time, 
  # ordered by mean ERK-KTR C:N of each cell
g <- g + geom_raster (aes (fill = Value), hjust = 0.5, vjust = 0.5, 
  interpolate = FALSE) + scale_x_continuous (expand = c (0, 0))
  # Heatmap of ERK-KTR C:N values
g <- g + viridis::scale_fill_viridis (option = "A", 
  rescaler = function(x, to = c (0, 1), from = NULL) {
  ifelse(x < 3, scales::rescale (x, to = to, 
  from = c (min (x, na.rm = TRUE), 3)), 1)
  })
  # Colour scale Viridis, set ERK-KTR C:N limits between 0 - 3
g <- g + facet_wrap (fct_rev (Treatment) ~ Cell_number ~ Embryo.Id, 
  scales = "free", ncol = 5)
  # Plot each heatmap by embryo, treatment, and cell number
g <- g + theme_classic() + 
  theme(aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7, angle = 90), 
  axis.ticks.y = element_blank (), 
  axis.text.y = element_blank (),
  legend.position = "bottom", 
  legend.text = element_text (size = 6), 
  legend.title = element_text(size = 6),  
  panel.background = element_rect (fill = "grey20", color  =  NA),
  panel.spacing.x = unit (.05, "lines"),
  panel.spacing.y = unit (.05, "lines"),
  strip.background = element_blank (),  
  strip.text = element_text (size = 6), 
  text = element_text (size = 8)) + 
  labs (x = "Time (min)", y = "ERK-KTR C:N") 
  # Modify plot appearance
print(g)

ggsave("SuppFig2C.tiff", width = 3) # Save as TIFF
ggsave("SuppFig2C.pdf", width = 3) # Save as PDF


#------------------------------------------------------------------------------# 
# Supplemental Figure 2E ----
# Boxplot of coefficient of variance per cell by treatment


g <- ggplot(data = subset(avg, TE.ICM =="ICM" & MtUniqueID %in% long_track & Cell_number >= 32), 
            aes(x = factor(Treatment, levels = c("Control", "1uM PD03")), y = trackCV))
g <- g + geom_boxplot(aes(fill= factor(Treatment, levels = c("Control", "1uM PD03"))), 
                      alpha = 0.5, outlier.size = 0, color = "black", size = 0.5)
g <- g + geom_jitter(aes(color = factor(Treatment, levels = c("Control", "1uM PD03"))), 
                     position=position_jitter(width=0.2), alpha=0.5,  size = 1, stroke = 0) 
g <- g + theme_classic() +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size=7), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 6),  
        legend.position="none", 
        axis.line = element_line(size = 0.5),
        strip.background = element_rect(size = 0.5),
        text = element_text(size=8),
        aspect.ratio = 1) + 
  labs(y = "CV per cell", x="") +
  scale_fill_discrete(name = "Treatment", labels = c("Control", "MEKi"))
print(g)

ggsave("SuppFig2E.tiff", width = 1.25) # Save as TIFF