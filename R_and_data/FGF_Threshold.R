#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#        Cell tracks thresholded in Control vs FGF treated embryos             #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          FGF_data.R
# Outputs: Supplemental Figure 2I
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('FGF_data.R')


## Data wrangling ####

palette <- c(DP = "purple", EPI = "red", PrE = "blue", "Polar TE" = "#80FF80", 
  "Mural TE" = "#00A500", ICM = "purple3")


#------------------------------------------------------------------------------# 
# ICM ----

d <- filter (embryos_CN_FGF, TE.ICM == "ICM")

good_spots <- filter (d, Value != 0 )
good_spots <- good_spots %>% 
  filter (!grepl ("Mitosis", Label))
good_spots <- good_spots %>% 
  filter (!grepl ("Apoptosis", Label))
good_spots <- good_spots %>% 
  filter (!grepl("Bad Segmentation", Label))

dat <- filter (good_spots, MtUniqueID %in% long_track)

FGF4 <- dat %>% 
  filter (Treatment =="FGF4")

dat$val_thresh <- ifelse (dat$Value > 1.2142, "HIGH", "LOW")
# threshold from Matlab Threshold_finder.m

summary <- dat %>% 
  group_by (MtUniqueID, Embryo.Id, Treatment, val_thresh) %>% 
  summarise (n = n()) %>% 
  mutate (freq = n / sum(n))

summary <- summary %>% 
  arrange (val_thresh, desc(freq)) 
first <- subset (summary, val_thresh == "LOW")
ls <- first$MtUniqueID
second <- subset (summary, val_thresh == "HIGH" & freq == 1)
ls2 <- second$MtUniqueID
sort <- c(ls, ls2)

summary$MtUniqueID <- factor (summary$MtUniqueID, levels = sort)


#------------------------------------------------------------------------------# 
# Supplemental Figure 2I ----

g <- summary %>% 
  ggplot (aes (x = MtUniqueID, y = freq, fill = factor(val_thresh)))
g <- g + geom_bar(stat = "identity", position = "fill")
g <- g + theme_classic () + 
  theme(aspect.ratio = 0.4, 
        axis.line = element_line (size = 0.5),
        axis.text.x = element_blank (),
        axis.text.y = element_text (size = 6), 
        axis.ticks = element_blank (), 
        legend.position = "top", 
        legend.text = element_text (size = 6), 
        legend.title = element_text(size = 6),  
        strip.background = element_rect(size = 0.5),
        text = element_text (size = 8)) +
  facet_wrap(. ~ Treatment, nrow = 1, scales = "free") + 
  labs (y = "% of cell track", x = "Individual cells") + 
  scale_fill_viridis_d (option = "A", begin = 0.3, end = 0.8, direction = -1) + 
  scale_y_continuous (labels = scales::percent)
print(g)

ggsave("SuppFig2I.tiff", width = 2)
ggsave("SuppFig2I.pdf", width = 4)

