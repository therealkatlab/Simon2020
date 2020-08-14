#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                     Lineages in 2hr time-lapse movies                       #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          Lineage_2hr_Data.R
#          totals1.csv
#          totals2.csv
# Outputs: Figure 3C
#          Figure 3D
#          Figure 3E
#          Figure 3G
#          Supplemental Figure 3B
#          Figure 4D
#          Figure 4E
#          Supplemental Figure 3B
#          Supplemental Figure 4A
#          Supplemental Figure 4C
#          Supplemental Figure 4D
#          Supplemental Figure 4E
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('Lineage_2hr_Data.R')

###


## Aesthetics for plots


palette <- c (ICM = "purple", DP = "purple", EPI = "red", PrE = "blue", 
  "Polar TE" = "#80FF80", "Mural TE" = "#00A500", "TE" = "green")


#------------------------------------------------------------------------------# 
# Figure 3C ----
# Line plot of TE ERK-KTR C:N values over time per cell, by stage and identity

# Plot line plot
g <- embryos_CN %>%  
  filter (TE.ICM == "TE" & MtUniqueID %in% long_track & Cell_number >= 32 & 
            !grepl ("Mitosis", Label)) %>%
  ggplot (aes (x = TimeM, y = Value))
# Plot ERK-KTR C:N values for ICM cells in blastocysts (>32c) tracked for 
# > 5 frames, during interphase, over time, excluding badly segmented cells
g <- g + geom_line (aes (group = MtUniqueID), color = "black", alpha = 0.1, 
                    size = 0.5)
# Plot line plot for individual cells
g <- g + facet_grid (Stage ~ Identity)
# Divide plots by stage and cell identity
g <- g + theme_classic () + scale_x_continuous (expand = c (0, 0)) +
  theme (axis.line = element_line (size = 0.5),
         axis.text = element_text (size = 7),
         legend.text = element_text (size = 6),
         legend.title = element_text (size = 6),
         legend.position = "top",
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 8)) + coord_fixed (25) +
  labs (x = 'Time (min)', y = 'ERK-KTR C:N') +  
  scale_color_manual (values = palette)
# Customise plot appearance
print (g)

ggsave ("Fig3C.tiff", width = 2.5) # Save TIFF


#------------------------------------------------------------------------------# 
# Figure 3D ----
# Violinplot of mean ERK-KTR C:N per TE cell, by identity

# Plot violin plot
g <- avg %>%
  filter (TE.ICM == "TE" & MtUniqueID %in% long_track & Cell_number >= 32) %>%
  ggplot (aes (x = Identity, y = trackAvg))
# Plot mean ERK-KTR C:N values for TE cells in blastocysts, 
# tracked for > 5 frames
g <- g + geom_jitter (aes (color = Identity), 
                      position = position_jitter (width = 0.2), alpha = 0.5,  size = 1, stroke = 0)
# Plot points for mean ERK-KTR C:N values in individual TE cells
g <- g + geom_violin (aes (fill = Identity), alpha = 0.5, color = "black", 
                      size = 0.5)
# Plot violin plots for mean ERK-KTR C:N values by TE cell identity
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
                       color = "black", size = 0.5, width = 0.1)
# Plot violin plots for mean ERK-KTR C:N values by TE cell identity
g <- g + theme_classic () + labs (y = "Mean per cell", x = "") +
  guides (fill = guide_legend (title = "Identity"), 
          color = guide_legend (title = "Identity")) +
  theme (aspect.ratio = 1, 
         axis.line = element_line (size = 0.5),
         axis.text.x = element_blank (), 
         axis.text.y = element_text (size = 7),
         legend.position = "right", 
         legend.text = element_text (size = 6), 
         legend.title = element_text (size = 6),  
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 8)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  facet_wrap ( ~ Stage)
# Customise plot appearance, and divide plots by stage
print (g)

ggsave ("Fig3D.tiff", width = 3) # save TIFF


#------------------------------------------------------------------------------# 
# Figure 3E, Figure 3G & Supplemental Figure 3B ----

# Calculate TE nearest ICM neighbor for each frame

# Group data into TE and ICM
TE <- subset (embryos_CN, TE.ICM == "TE")
ICM <- subset (embryos_CN, TE.ICM == "ICM")

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
    sub <- ICM.loc[j,]
    d.X <- sub$X - TE.loc$X # Caclulate distance in each dimension
    d.Y <- sub$Y - TE.loc$Y
    d.Z <- sub$Z - TE.loc$Z
    d <- sqrt (d.X ^ 2 + d.Y ^ 2 + d.Z ^ 2) # Calculate distance in 3D
    a <- rbind(a, d) 
  }
  
  # Find the distance between the TE cell and its nearest neighbor ICM
  nn <- min(a)
  
  # Aggregate data on distance betweeen TE cell and its nearest ICM neighbor
  spoti$d <- nn
  dist <- rbind(dist, spoti)
}


#------------------------------------------------------------------------------# 
# Figure 3E ----

delta <- dist %>% 
  filter (Value != 0 & (!grepl("Mitosis", Label)) & 
            (!grepl ("Apoptosis", Label))) %>%
  group_by (MtUniqueID) %>% 
  mutate (d.intercept = t (coef (lm (d ~ TimeM)))[ , 1], 
          d.slope = t (coef (lm (d ~ TimeM)))[ , 2], 
          v.intercept = t (coef (lm (Value ~ TimeM)))[ , 1], 
          v.slope = t (coef (lm (Value ~ TimeM)))[ , 2])

delta_summary <- delta %>% 
  group_by(MtUniqueID, Cell_number, Identity, d.slope, Stage) %>% 
  summarise()

g <- delta_summary %>%
  filter (Cell_number >= 32 & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = d.slope*60))
g <- g + geom_density (aes (fill = Identity), color = "black", alpha = 0.75)
g <- g + geom_vline (xintercept = 0, linetype = 3) 
g <- g + theme_classic ()  +
  theme(aspect.ratio = 0.75,
        axis.line = element_line (size = 0.5),
        axis.text = element_text (size = 6),
        legend.position = "top",
        legend.text = element_text (size = 6),
        legend.title = element_text (size = 6),
        strip.background = element_rect (size = 0.5),
        text = element_text (size = 8))  +  
  facet_grid ( . ~ Stage) +
  labs (x = expression (Delta*"distance from nearest ICM cell/"*Delta*"time (μm/hr)")) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  scale_y_continuous (breaks = seq (0, 0.2, by = 0.1))
print(g)

ggsave("Fig3E.tiff", width = 2.5)


#------------------------------------------------------------------------------# 
# Figure 3G ----
# Scatter plot of mean distance of TE to nearest ICM neighbor against mean 
# ERK-KTR C:N value per cell

# Summary stats on TE-ICM nearest neighbor distances
summary_TE <- dist %>% 
  group_by (MtUniqueID, Identity, Embryo.Id, trackAvg, Cell_number, Stage) %>% 
  summarise (d = mean (d))
summary_TE <- subset (summary_TE, MtUniqueID %in% long_track & 
                        Cell_number >= 32)
# summary_emb <- dist %>% 
# group_by (Identity, Embryo.Id, Cell_number, Stage) %>% 
# summarise( d = mean (d), mean = mean (trackAvg))

# Plot scatter plot
g <- summary_TE %>% 
  ggplot (aes (x = d, y = log (trackAvg)))
# Plot individual TE cells mean distance from nearest ICM neighbor and 
# log mean ERK-KTR C:N values
g <- g + geom_point (size = 0.5, alpha = 0.5) 
# Colour code by treatment
g <- g + geom_smooth (method = "lm", color = "red") 
# Plot linear model
g <- g + theme_classic ()  +
  theme (aspect.ratio = 1, 
         axis.line = element_line (size = 0.5),
         axis.text = element_text (size = 7),
         legend.position = "top",
         legend.text = element_text (size = 6),
         legend.title = element_text (size = 6),
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 8)) + 
  labs (y = "Log[mean ERK activity]", x = "Mean distance from ICM (μM)") +  
  scale_color_manual (values = palette)
# Customise plot appearance
print(g)

ggsave("Fig3G.tiff", width = 3) # Save TIFF


# Linear model to look at the effect of distance on log mean ERK activity
# (i.e. exponential decay in ERK activity with increasing distance from ICM)
model.lm <- lm (log (trackAvg) ~ d, data = summary_TE)
summary (model.lm) # slope -0.0175105 r2 0.4521 p-value: < 2.2e-16


#------------------------------------------------------------------------------# 
# Supplemental Figure 3B ----
# Line plot of mural and polar TE population mean distance from ICM over time

# Plot line plot
g <- dist %>% 
  filter (Cell_number >= 32 & MtUniqueID %in% long_track) %>%
  ggplot (aes (x = TimeM, y = d))
# Plot TE cells  distance from nearest ICM neighbor in TE cells tracked for 
# > 5 frames in blastocysts
g <- g + stat_smooth (aes (color = Identity)) 
# Plot population mean for Mural and Polar TE identities
g <- g + theme_classic ()  +
  theme(aspect.ratio = 1,
        axis.line = element_line(size = 0.5),
        axis.text = element_text (size = 7),
        legend.position = "top",
        legend.text = element_text (size = 6),
        legend.title = element_text (size = 6),
        strip.background = element_rect(size = 0.5),
        text = element_text (size = 8) ) + 
  facet_wrap( . ~ Stage) +
  labs(x = "Time (min)", y = "Distance from ICM (μm)") +  
  scale_color_manual (values = palette)
# Customise plot appearance and divide by stage
print(g)

ggsave("FigS3B.tiff", width = 3) # Save TIFF


#------------------------------------------------------------------------------# 
# Supplemental Figure 4A ----
# Heatmap of ICM ERK-KTR C:N values over time, by embryo

# Duplicate parent tracks so they are double height of daughters
# This allows parent track to "split" into two daughters in the heatmap
# and keeping the alignments in-line

parent <- filter (embryos_CN, Mitosis == "Parent")
parent$dMtUniqueID <- paste0 (parent$UniqueID, 1)

parent2 <- filter (embryos_CN, Mitosis == "Parent")
parent2$dMtUniqueID <- paste0 (parent2$UniqueID, 2)

d1 <- filter (embryos_CN, Mitosis == "Daughter1")
d1$dMtUniqueID <- paste0 (d1$UniqueID, 1)

d2 <- filter(embryos_CN, Mitosis == "Daughter2")
d2$dMtUniqueID <- paste0 (d2$UniqueID, 2)

duplicate <- rbind (parent, parent2, d1, d2)

# Plot heatmap
g <- duplicate %>% 
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track & Cell_number >= 32 & 
  !grepl("Mitosis", Label)) %>% 
  ggplot (aes (x = TimeM, y = reorder (MtUniqueID, trackAvg)))
  # Plot individual ICM cells in blastocysts (>32c), during interphase, 
  # over time, ordered by mean ERK-KTR C:N of each cell
g <- g + geom_raster (aes (fill = Value), hjust=0.5, vjust=0.5, 
  interpolate=FALSE) + scale_x_continuous (expand = c (0, 0))
  # Heatmap with ERK-KTR C:N values
g <- g + viridis::scale_fill_viridis (option = "A", 
  rescaler = function(x, to = c (0, 1), from = NULL) {
  ifelse (x < 3, scales::rescale (x, to = to, from = c (min (x, na.rm = TRUE), 
  3)), 1)
  })
  # Set scale for heatmap
g <- g + facet_wrap (Cell_number ~ ., scales = "free", ncol = 6)
  # Divide and order plots by embryo and cell number
g <- g + theme_classic () +
  theme (aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 6),
  axis.text.y = element_blank (), 
  axis.ticks.y = element_blank (), 
  panel.background = element_rect (fill = "grey20", color  =  NA),
  panel.border = element_blank (),
  panel.grid.major = element_blank (),
  panel.grid.minor = element_blank (),
  panel.spacing.x = unit (.1, "lines"),
  panel.spacing.y = unit (.1, "lines"),
  legend.position = "bottom", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_blank (), 
  strip.text = element_blank (), 
  text = element_text (size = 8)) +
  labs(x = "Time (min)", y = "ERK-KTR C:N")
  # Customise plot appearance
print(g)

ggsave("FigS4A.tiff", width = 3) # save TIFF


#------------------------------------------------------------------------------# 
# Figure 4D ----
# Line plot of ICM ERK-KTR C:N values over time per cell, by stage and identity

# Plot line plot
g <- embryos_CN %>%  
  filter (TE.ICM == "ICM" & Identity != "ICM" & MtUniqueID %in% long_track & 
  Cell_number >= 32 & !grepl ("Mitosis", Label) & 
  !grepl ("Bad Segmentation", Label)) %>%
  ggplot (aes (x = TimeM, y = Value))
  # Plot ERK-KTR C:N values for ICM cells in blastocysts (>32c), with known cell
  # identities, tracked for > 5 frames, during interphase, over time,
  # excluding badly segmented cells
g <- g + geom_line (aes (group = MtUniqueID), color = "black", alpha = 0.1, 
  size = 0.5)
  # Plot line plot for individual cells
g <- g + geom_smooth (aes (x = TimeM, y = Value, color = Identity), size = 1, 
  fill = "white", se = TRUE) + 
  guides (color = guide_legend (title = "Population Mean"))
  # Overlay population mean for each cell identity
g <- g + facet_grid (Stage ~ Identity)
  # Divide plots by stage and cell identity
g <- g + theme_classic () + scale_x_continuous (expand = c (0, 0)) +
  theme (axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  legend.position = "top",
  strip.background = element_rect(size = 0.5),
  text = element_text (size = 8)) + coord_fixed (40) +
  labs(x = 'Time (min)', y = 'ERK-KTR C:N') +  
  scale_color_manual (values = palette)
  # Customise plot appearance
print (g)

ggsave ("Fig4D.tiff", width = 3)


# Calculate percentage change in ERK-KTR C:N over movie from first to last frame

# Subset first frame and last frame for ICM cells with known cell identities at 
# blastocyst stage, tracked for > 5 frames,  excluding cells with poor 
# segmentation or mitosis at these time points

clean_t1 <- embryos_CN %>% 
  filter (TE.ICM == "ICM" & Identity != "ICM" & MtUniqueID %in% long_track & 
  Cell_number >= 32 & Value > 0 & !grepl("Mitosis", Label) & 
  !grepl("Bad Segmentation", Label) & Time == 1)

clean_t24 <- embryos_CN %>% 
  filter (TE.ICM == "ICM" & Identity != "ICM" & MtUniqueID %in% long_track & 
  Cell_number >= 32 & Value > 0 & !grepl("Mitosis", Label) & 
  !grepl("Bad Segmentation", Label) & Time == 24)

# Identify those cells tracked at both the beginning and end of the time-lapse

frame1_ID <- clean_t1$MtUniqueID
frame24_ID <- clean_t24$MtUniqueID
long_ID <- frame1_ID[frame24_ID %in% frame1_ID]

# Keep data for those cells present at begininng and end of time-lapse

clean_t1 <- clean_t1 %>% 
  filter (MtUniqueID %in% long_ID)

clean_t24 <- clean_t24 %>% 
  filter (MtUniqueID %in% long_ID)

# Keep key variables and bind first and last time-points

vars <- c("Value", "Identity", "Stage", "MtUniqueID", "Time")
clean_t1 <- clean_t1[vars]
clean_t24 <- clean_t24[c("MtUniqueID", "Value")]

together <- left_join(clean_t1, clean_t24, by="MtUniqueID")

# Calculate percentage change for each individual cell

together <- together %>% 
  mutate(pc = ( (Value.y - Value.x) / Value.x) * 100)

together <- together %>% 
  filter(!is.na(pc))

# Calculate mean percentage stage for each stage

together %>% 
  group_by (Stage) %>% summarise (mean = mean(pc)) 
  # 32-64 = 17.9% 
  # 64-128 = 0.642%


#------------------------------------------------------------------------------# 
# Figure 4E ----
# Boxplot of mean ERK-KTR C:N per ICM cell, by identity

# Plot boxplot
g <- avg %>%
  filter (TE.ICM =="ICM" & Cell_number >= 32 & Identity != "ICM" & 
  MtUniqueID %in% long_track) %>%
  ggplot (aes (x = Identity, y = trackAvg))
  # Plot mean ERK-KTR C:N values for ICM cells in blastocysts, 
  # tracked for > 5 frames, with known cell identities
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
  # Plot boxplots by cell identity
g <- g + geom_jitter (aes (color = Identity), 
  position = position_jitter (width = 0.2), alpha = 0.5,  size = 1, stroke = 0)
  # Overlay with points for mean ERK-KTR C:N values in individual cells
g <- g + theme_classic () + labs (y = "Mean per cell", x = "") + 
  guides (fill = guide_legend (title = "Identity"), 
  color = guide_legend (title = "Identity")) +
  theme (aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_blank (), 
  axis.text.y = element_text (size = 7),
  legend.position = "right", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  facet_wrap ( ~ Stage)
  # Customise plot appearance, and divide plots by stage
print(g)

ggsave("Fig4E.tiff", width = 3) # Save TIFF


# Linear mixed model to test identity

long <- subset (avg, MtUniqueID %in% long_track & 
                  (Identity == "PrE" | Identity == "DP" | Identity == "EPI"))
early <- subset(long, Stage == "32_64")
late <- subset(long, Stage == "64_128")

model.lmer <- lmer (trackAvg ~ Identity + (1 | Embryo.Id / Identity), 
                    data = early)
summary(model.lmer)
anova(model.lmer) # Identity p = 0.02017
rand(model.lmer) # Random effects n.s.
marginal = lsmeans (model.lmer, ~ Identity)
pairs(marginal, adjust = "tukey") # EPI-PrE 0.0271

model.lmer <- lmer (trackAvg ~ Identity + (1 | Embryo.Id / Identity), 
                    data = late) 
summary (model.lmer)
anova (model.lmer) # Identity p = 0.04617
rand (model.lmer) # Random effects n.s.
marginal = lsmeans (model.lmer, ~ Identity)
pairs (marginal, adjust = "tukey") # Pairwaise post-hoc n.s.


#------------------------------------------------------------------------------# 
# Supplemental Figure 4C ----
# Crossbar and dotplot of mean ERK-KTR C:N per ICM cell, by identity, 
# for each embryo

# Plot crossbar and dotplot
g <- avg %>%
  filter (TE.ICM =="ICM" & Cell_number >= 32 & Identity != "ICM" & 
  MtUniqueID %in% long_track) %>%
  ggplot (aes (x = Identity, y = trackAvg))
 # Plot mean ERK-KTR C:N values for ICM cells in blastocysts, 
 # tracked for > 5 frames, with known cell identities
g <- g + stat_summary (aes (color = Identity), fun.y = mean, fun.ymin = mean, 
  fun.ymax = mean, geom = "crossbar", width = 0.8, alpha = 0.5, size = 0.2)
  # Plot crossbar for mean ERK-KTR C:N values in each cell identity
g <- g + geom_jitter (aes (color = Identity), 
  position = position_jitter (width = 0.2), alpha = 0.4,  size = 1, 
  stroke = 0.5) 
  # Overlay with points for mean ERK-KTR C:N values in individual cells of 
  # each cell identity
g <- g + theme_classic () + labs (y = "Mean per cell", x = "") + 
  guides (fill = guide_legend (title = "Identity"), 
  color = guide_legend (title = "Identity")) +
  theme(aspect.ratio = 1.5,
  axis.line = element_line (size = 0.5),
  axis.line.x = element_blank (),
  axis.text.x = element_blank (), 
  axis.text.y = element_text (size = 7),
  axis.ticks.x = element_blank (),
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  legend.position = "bottom", 
  strip.background = element_blank (),
  text = element_text (size = 8)) +  
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  facet_wrap (Cell_number ~ Embryo.Id, ncol = 14)
  # Customise plot appearance, and divide plots by embryo ordered by cell number
print (g)

ggsave ("Fig4SC.tiff", width = 4) # save TIFF
ggsave ("Fig4SC.pdf", width = 7) # save PDF


#------------------------------------------------------------------------------# 
# Supplemental Figure 4D & Suppemental Figure 4E ----

# Get lineage info from MINS
totals <- read_csv("totals1.csv")
totals2 <- read_csv("totals2.csv")

# Calculate total cells in each lineage
total_all <- rbind (totals, totals2)

id <- embryos_CN %>% 
  group_by (Embryo.Id) %>% 
  summarise()
ls_id <- id$Embryo.Id
sub_tot <- subset (total_all, Embryo.Id %in% ls_id)
icm <- sub_tot %>% 
  filter (Identity.km == "TE") %>%
  mutate (ICMcount = Cell_number - value)
icm2 <- icm[c ("Embryo.Id", "ICMcount")]
tot <- left_join (sub_tot, icm2, by = "Embryo.Id")
epi <- tot %>%
  filter (Identity.km == "EPI") %>%
  mutate (pc.EPI = value / ICMcount)
pc <- epi[c ("Embryo.Id", "pc.EPI")]

avg_icm <- left_join(avg, pc, by="Embryo.Id")

emb <- avg_icm %>% 
  group_by(Embryo.Id, Identity, Cell_number, Stage, pc.EPI) %>%
  summarise(mean=mean(trackAvg))

#------------------------------------------------------------------------------# 
# Supplemental Figure 4D ----

g <- emb %>% 
  filter (Cell_number >= 32 & 
  (Identity == "EPI" | Identity == "PrE" | Identity == "DP")) %>%
  ggplot (aes (x = Identity, y = mean))
g <- g + geom_boxplot (aes (fill = Identity), alpha = 0.5, outlier.size = 0, 
  color = "black", size = 0.5)
g <- g + geom_jitter (aes (color = Identity), 
  position = position_jitter (width = 0.2), alpha = 0.5, size = 1, stroke = 0) 
g <- g + theme_classic () + 
  theme(aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text.x = element_blank (), 
  axis.text.y = element_text (size = 7),
  legend.position = "right", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) + 
  labs (y = "Mean per embryo", x = "") + 
  guides (fill = guide_legend (title = "Identity"), 
  color = guide_legend (title = "Identity")) +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) +
  facet_wrap ( ~ Stage)
print (g)

ggsave ("FigS4D.tiff", width = 3)



#------------------------------------------------------------------------------# 
# Supplemental Figure 4E ----

g <- emb %>% 
  filter (Cell_number >= 32 & 
  (Identity == "EPI" | Identity == "PrE" | Identity == "DP")) %>%
  ggplot (aes (x = (pc.EPI * 100), y = mean))
g <- g + geom_smooth (method = "lm", color = "grey30") 
g <- g + geom_point (aes (color = Identity), size = 1, alpha = 0.5) 
g <- g + theme_classic () + facet_grid (  . ~ Identity)
g <- g + labs(y = "Mean per embryo", x = "% EPI in ICM") +
  theme (aspect.ratio = 1,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 7), 
  axis.text.y = element_text (size = 7),
  legend.position = "none", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) + 
  scale_color_manual (values = palette)
print(g)

ggsave("FigS4E.tiff", width = 4)


model.lm <- lm (mean ~ pc.EPI, data = subset (emb, Identity == "EPI"))
summary(model.lm) # p-value: 0.698 adj r2 -0.05956 slope 0.09041

model.lm <- lm (mean ~ pc.EPI, data = subset (emb, Identity == "PrE"))
summary(model.lm) #p-value: 0.002382 adj r2 0.4135 slope =  0.70335

model.lm <- lm(mean ~ pc.EPI, data = subset(emb, Identity == "DP"))
summary(model.lm) #p-value: 0.7231 adj r2  -0.05748 slope = 0.07888