#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                      Lineages in 12hr time-lapse movies                      #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          Lineage_2hr_Data.R
# Outputs: Supplemental Figure 3C
#          Supplemental Figure 3D
#          Supplemental Figure 3E
#          Supplemental Figure 3F
#          Supplemental Figure 4F
#          Supplemental Figure 4G
#          Supplemental Figure 4H
#          Supplemental Figure 4I
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('Lineage_12hr_Data.R')

## Aesthetics for plots

palette <- c(DP = "purple", EPI = "red", PrE = "blue", "Polar TE" = "#80FF80", 
  "Mural TE" = "#00A500", ICM = "grey", "DP&EPI" = "#FF00FF", 
  "DP&PrE" = "#7700FF", "EPI&PrE" = "purple", "EPI&EPI" = "red", 
  "PrE&PrE" = "blue", "DP&DP" = "purple", "NA" = "white")


#------------------------------------------------------------------------------# 
# Supplemental Figure 3C ----

g <- embryos_CN %>% 
  filter (TE.ICM == "TE" & MtUniqueID %in% long_track & Cell_number >= 32 & 
            Cell_number < 64 & Value > 0 & !grepl("Mitosis", Label)) %>%
  ggplot (aes (x = TimeM, y = Value))
g <- g + geom_smooth (aes (color = Identity), method = loess, size = 1, 
                      fill = "grey", se = TRUE) + 
  guides (color = guide_legend (title = "Population Mean"))
g <- g + theme_classic () +
  theme (aspect.ratio = 1, 
         axis.line = element_line (size = 0.5),
         axis.text = element_text (size = 7),
         legend.position = "top",
         legend.text = element_text (size = 6),
         legend.title = element_text (size = 6),
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 8)) +
  scale_color_manual (values = palette) + 
  scale_x_continuous (expand = c (0, 0)) +
  labs(x = "Time (min)", y = "ERK-KTR C:N") 
print(g)
ggsave("SuppFig3C.tiff", width = 3)



#------------------------------------------------------------------------------# 
# Supplemental Figure 3D ----


# Spatial position of TE
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


g <- dist %>% filter (Cell_number >= 32 & Cell_number < 64 & 
  MtUniqueID %in% long_track) %>% 
  ggplot (aes (x = TimeM, y = d))
g <- g + stat_smooth (aes (color = Identity)) 
g <- g + theme_classic ()  +
  theme (aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "top",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +
  labs (x = 'Time (min)', y = 'Distance from ICM (μm)') +  
  scale_color_manual (values = palette)
print (g)
ggsave ("SuppFig3D.tiff", width = 1.5)

#------------------------------------------------------------------------------# 
# Supplemental Figure 3E ----

delta <- dist %>% filter (Value != 0 & (!grepl ("Mitosis", dist$Label)) & 
  (!grepl ("Apoptosis", dist$Label))) %>%
  group_by (MtUniqueID) %>% 
  mutate (d.intercept = t (coef (lm (d ~ TimeM)))[,1], 
  d.slope = t (coef (lm (d ~ TimeM)))[,2], 
  v.intercept = t (coef (lm (Value ~ TimeM)))[,1], 
  v.slope = t (coef (lm (Value ~ TimeM)))[,2])

delta_summary <- delta %>% 
  group_by (MtUniqueID, Cell_number, Identity, d.slope, Stage) %>% summarise ()

g <- delta_summary %>% filter (Cell_number >= 32 & Cell_number < 64 & 
  (MtUniqueID %in% long_track)) %>%
  ggplot (aes (x = d.slope * 60))
g <- g + geom_density (aes (fill = Identity), color = "black", alpha = 0.75)
g <- g + geom_vline (xintercept = 0, linetype=3) 
g <- g + theme_classic ()  +
  theme(aspect.ratio = 0.75, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "top",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8))  +
  labs(x = expression (Delta * "distance from nearest ICM cell/" * Delta * "time (μm/hr)")) +
  scale_color_manual (values = palette) +  
  scale_fill_manual (values = palette)
print(g)
ggsave("SuppFig3E.tiff", width = 2)

#------------------------------------------------------------------------------# 
# Supplemental Figure 3F ----

summary_TE <- dist %>% 
  group_by (MtUniqueID, Identity, Embryo.Id, trackAvg, Cell_number, Stage) %>% 
  summarise (d = mean(d))
summary_TE <- subset (summary_TE, MtUniqueID %in% long_track & 
                       Cell_number >= 32 & Cell_number < 64)

g <- summary_TE %>% 
  filter (Cell_number >= 32 & Cell_number < 64) %>%
  ggplot (aes (x = d, y = log(trackAvg)))
g <- g + geom_point (size = 0.5, alpha = 0.5) 
g <- g + geom_smooth (method = "lm", color = "red") 
g <- g + theme_classic ()  +
  theme (aspect.ratio = 1, 
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "top",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +
  labs (y = "Log[mean ERK activity per cell]", 
  x = "Mean distance from ICM over time (μM)") +  
  scale_color_manual (values = palette)
print(g)
ggsave("SuppFig3F.tiff", width = 1.5)


# statistics
model.lm <- lm (log (trackAvg) ~ d, 
  data = subset (summary_TE, Cell_number >= 32 & Cell_number < 64))
summary (model.lm) # -0.0084013 r2 = 0.4892 p-value: < 2.2e-16

#------------------------------------------------------------------------------# 
# Supplemental Figure 4F ----

g <- embryos_CN %>%
  filter(MtUniqueID %in% long_track & 
  (End_Identity == "DP" | End_Identity == "EPI" | End_Identity == "PrE") & 
  Cell_number >= 32 & Cell_number < 64 & Value > 0 & 
  !grepl("Mitosis", Label) & 
  !grepl("Bad Segmentation", Label) & 
  !grepl("Apoptosis", Label)) %>% 
  ggplot (aes( x = (TimeM / 60), y = Value))
g <- g + geom_line (aes (group = MtUniqueID), color = "black", 
  alpha = 0.1, size = 0.5)
g <- g + geom_smooth (aes (color = End_Identity), method = loess, size = 1, 
  fill = "white", se = TRUE) + 
  guides (color = guide_legend (title = "Population Mean"))
g <- g + facet_grid ( ~ End_Identity)
g <- g + theme_classic() + 
  theme (aspect.ratio = 0.6,
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "top",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect(size = 0.5),
  text = element_text (size = 8)) + 
  labs (x = "Time (hr)", y = "ERK-KTR C:N") +  
  scale_color_manual (values = palette) + 
  scale_x_continuous (expand = c (0, 0))
print (g)
ggsave ("SuppFig4F.tiff", width = 3)


#------------------------------------------------------------------------------# 
# Supplemental Fig 4G ----

g <- avg %>%
  filter (MtUniqueID %in% long_track & 
            (End_Identity == "DP" | End_Identity == "EPI" | End_Identity == "PrE") & 
            Cell_number >= 32 & Cell_number < 64) %>%
  ggplot(aes(x = End_Identity, y = trackAvg))
g <- g + geom_boxplot (aes (fill = End_Identity), alpha = 0.5, outlier.size = 0, 
                       color = "black", size = 0.5)
g <- g + geom_jitter (aes (color = End_Identity), 
                      position = position_jitter (width = 0.2), alpha = 0.5, size = 1, stroke = 0) 
g <- g + theme_classic () + 
  theme (aspect.ratio = 1,
         axis.line = element_line (size = 0.5),
         axis.text.x = element_text (size = 7, angle = 90), 
         axis.text.y = element_text (size = 7),
         legend.position = "top", 
         legend.text = element_text (size = 6), 
         legend.title = element_text (size = 6),  
         strip.background = element_rect (size = 0.5),
         text = element_text (size = 8))  + 
  labs (x = "", y = "Mean per cell") + 
  guides (fill = guide_legend (title = "Identity"), 
          color = guide_legend (title = "Identity")) +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) 
print (g)
ggsave ("SuppFig4G.tiff", width = 1.5)

# stats test

endid <- subset (avg, MtUniqueID %in% long_track & TE.ICM == "ICM" & 
                   (End_Identity == "DP" | End_Identity == "EPI" | End_Identity == "PrE") & 
                   Cell_number >= 32 & Cell_number < 64)

model.lmer <- lmer (trackAvg ~ End_Identity + (1 | Embryo.Id / End_Identity), 
                    data = endid)
summary (model.lmer)
anova (model.lmer)
rand (model.lmer)

marginal = lsmeans (model.lmer, ~ End_Identity)
pairs (marginal, adjust = "tukey") # EPI PRE P = 0.0191


#------------------------------------------------------------------------------# 
# Supplemental Figure 4H ----

g <- CN_mitosis_t %>%
  filter (TE.ICM == "ICM" & MtUniqueID %in% long_track & 
  Cell_number >= 32 & Cell_number < 64 & !grepl ("Mitosis", Label) & 
  !grepl ("Bad Segmentation", Label) & !grepl ("Apoptosis", Label)) %>%
  ggplot (aes (x = MitosisRT / 60, y = Value))
g <- g + geom_line (aes (group = MtUniqueID), color = "grey40", alpha = 0.1, 
  size = 0.5)
g <- g + geom_vline (aes (xintercept = 0), linetype = "dotted")
g <- g + geom_smooth(data=.%>%filter(MitosisRT < 0), color = "black", method = loess, 
  size = 1, fill = "grey", se = TRUE) + 
  guides (color = guide_legend (title = "Population Mean"))
g <- g + geom_smooth(data=.%>%filter(MitosisRT >= 0), color = "black", method = loess, 
  size = 1, fill = "grey", se = TRUE) + 
  guides (color = guide_legend (title = "Population Mean"))
g <- g + theme_classic () + 
  theme (aspect.ratio = 0.75,
  axis.line = element_line (size = 0.5),
  axis.text = element_text (size = 7),
  legend.position = "top",
  legend.text = element_text (size = 6),
  legend.title = element_text (size = 6),
  strip.background = element_rect (size = 0.5),
  text = element_text (size = 8)) +  scale_x_continuous (expand = c (0, 0)) +
  labs (x = "Relative Mitotic Time (hr)", y = "ERK-KTR C:N")
print (g)
ggsave ("SuppFig4H.tiff", width = 2)


#------------------------------------------------------------------------------# 
# Supplemental Figure 4I ----
# bin over time of mitosis, 2 hrs

similar <- c("EPI","EPI&EPI","PrE","PrE&PrE","DP","DP&DP")

CN_mitosis_t$same_ID <- ifelse (CN_mitosis_t$End_Identity == "EPI&EPI" | 
                                 CN_mitosis_t$End_Identity == "EPI", "EPI",
                         ifelse (CN_mitosis_t$End_Identity == "PrE&PrE" |
                                 CN_mitosis_t$End_Identity == "PrE", "PrE",
                          ifelse (CN_mitosis_t$End_Identity == "DP&DP" |
                                CN_mitosis_t$End_Identity == "DP", "DP", "ICM")
                          ))

# specify interval/bin labels
tags <- c("[-12_-10)","[-10_-8)","[-8_-6)", "[-6_-4)", "[-4_-2)", "[-2_0)", 
  "[0_2)","[2_4)", "[4_6)", "[6_8)", "[8_10)","[10_12)")

# bucketing values into bins
bin <- CN_mitosis_t %>% 
  mutate(tag = case_when(
    MitosisRT < -600                     ~ tags[1],
    MitosisRT >= -600 & MitosisRT < -480 ~ tags[2],
    MitosisRT >= -480 & MitosisRT < -360 ~ tags[3],
    MitosisRT >= -360 & MitosisRT < -240 ~ tags[4],
    MitosisRT >= -240 & MitosisRT < -120 ~ tags[5],
    MitosisRT >= -120 & MitosisRT < 0    ~ tags[6],
    MitosisRT >= 0    & MitosisRT < 120  ~ tags[7],
    MitosisRT >= 120  & MitosisRT < 240  ~ tags[8],
    MitosisRT >= 240  & MitosisRT < 360  ~ tags[9],
    MitosisRT >= 360  & MitosisRT < 480  ~ tags[10],
    MitosisRT >= 480  & MitosisRT < 600  ~ tags[11],
    MitosisRT >= 600                     ~ tags[12]
  ))

bin.Avg <- bin %>% 
  group_by (MtUniqueID, Cell_number, Stage, Identity, End_Identity, same_ID, 
  tag) %>% 
  filter (!grepl("Mitosis", Label) & !grepl("Bad Segmentation", Label) & 
  MtUniqueID %in% long_track) %>% 
  summarise (n = n (), mean = mean (Value))
bin.Avg <- bin.Avg %>% 
  filter( n >= 3)
bin.Avg$tag <- factor (bin.Avg$tag, levels = tags)

# Boxplot binned by mitotic time

end_tags <- c ("[0_2)", "[2_4)", "[4_6)", "[6_8)", "[8_10)", "[10_12)")

g <- bin.Avg %>%
  filter (MtUniqueID %in% long_track & 
  (End_Identity == "EPI" | End_Identity == "PrE"| End_Identity == "DP") &
  Cell_number >= 32 & Cell_number < 64 & (tag %in% end_tags)) %>%
  ggplot (aes (x = tag, y = mean)) 
g <- g + geom_point (aes (group = End_Identity, color = End_Identity), 
  position = position_dodge (width = 0.75), alpha = 0.5, size = 0.5)
g <- g + geom_boxplot (aes (x = tag, fill = End_Identity), alpha = 0.5, 
  outlier.size = 0, color = "black", size = 0.5)
g <- g + theme_classic () + 
  theme (aspect.ratio = 0.6,
  axis.line = element_line (size = 0.5),
  axis.text.x = element_text (size = 7, angle = 90), 
  axis.text.y = element_text (size = 7),
  legend.position = "right", 
  legend.text = element_text (size = 6), 
  legend.title = element_text (size = 6),  
  strip.background = element_rect (size = 0.5),
  text = element_text(size=8)) + 
  labs (y = "Mean per cell", x = "Time since mitosis (hr)") + 
  guides (fill = guide_legend (title = "Identity"), 
  color = guide_legend (title = "Identity")) +
  scale_color_manual (values = palette) + 
  scale_fill_manual (values = palette) 
print(g)
ggsave("SuppFig4I.tiff", width = 3.5)