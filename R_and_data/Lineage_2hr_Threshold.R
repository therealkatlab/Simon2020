#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#             Cell tracks thresholded in 2hr time-lapse movies                 #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Packages.R
#          Lineage_2hr_Data.R
# Outputs: Supplemental Figure 4B
#          Supplemental Figure 3A
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

## Load packages, themes and data from other file ----
setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('Lineage_2hr_Data.R')


#------------------------------------------------------------------------------# 
# TE ----

d <- filter (embryos_CN, TE.ICM == "TE")

good_spots <- filter (d, Value != 0 )
good_spots <- good_spots %>% 
  filter (!grepl("Mitosis", good_spots$Label))
good_spots <- good_spots %>% 
  filter (!grepl("Apoptosis", good_spots$Label))
good_spots <- good_spots %>% 
  filter (!grepl("Bad Segmentation", good_spots$Label))
dat <- filter (good_spots, MtUniqueID %in% long_track)

dat$val_thresh <- ifelse(dat$Value > 1.2142, "HIGH", "LOW")

summary <- dat %>% 
  group_by (MtUniqueID, Embryo.Id, Identity, Cell_number, Stage, val_thresh) %>% 
  summarise (n = n()) %>%   mutate (freq = n / sum(n))

summary <- summary %>% 
  arrange(val_thresh, desc(freq)) 
first <- subset (summary, val_thresh == "LOW")
ls <- first$MtUniqueID
second <- subset (summary, val_thresh == "HIGH" & freq == 1)
ls2 <- second$MtUniqueID
sort <- c(ls, ls2)

summary$MtUniqueID <- factor (summary$MtUniqueID, levels = sort)


#------------------------------------------------------------------------------# 
# Supplemental Figure 3A ----

g <- summary %>% 
  # filter (Stage == "32_64") %>% # early blastocyst stage
  filter (Stage == "64_128") %>% # mid blastocyst stage
  ggplot(aes(x = MtUniqueID, y = freq, fill = factor(val_thresh)))
g <- g + geom_bar(stat = 'identity', position = 'fill')
g <- g + theme_classic() + 
  theme(aspect.ratio = 0.4,
        axis.text.x = element_blank (),
        axis.ticks = element_blank (), 
        legend.position = "top", 
        legend.title = element_blank(),  
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.spacing.x = unit (.1, "lines"),
        panel.spacing.y = unit (.1, "lines"),
        strip.background = element_blank (), 
        strip.text = element_blank ()) + 
  scale_y_continuous (labels = scales::percent) + 
  facet_wrap (Identity ~ ., nrow = 3, scales = "free") + 
  scale_fill_viridis_d (option = "A", begin = 0.3, end = 0.8, direction = -1) + 
  labs(y = "% of cell track", x = "Individual cells")
print (g)

save("SuppFig3A.tiff", width = 2)
ggsave("SuppFig3A.pdf", width = 2)


#------------------------------------------------------------------------------# 
# ICM ----

# clearn up data 

fate <- c("DP", "EPI", "PrE")

d <- filter (embryos_CN, Identity %in% fate)

good_spots <- filter (d, Value != 0 )
good_spots <- good_spots %>%
  filter (!grepl("Mitosis", good_spots$Label))
good_spots <- good_spots %>% 
  filter(!grepl("Apoptosis", good_spots$Label))
good_spots <- good_spots %>% 
  filter(!grepl("Bad Segmentation", good_spots$Label))
dat <- filter(good_spots, MtUniqueID %in% long_track)

dat$val_thresh <- ifelse(dat$Value > 1.2142, "HIGH", "LOW") 
# threshold from Matlab

summary <- dat %>% 
  group_by (MtUniqueID, Embryo.Id, Identity, Cell_number, Stage, val_thresh) %>% 
  summarise (n = n()) %>%   mutate (freq = n / sum(n))

summary <- summary %>% 
  arrange (val_thresh, desc(freq)) 
first <- subset(summary, val_thresh == "LOW")
ls <- first$MtUniqueID
second <- subset (summary, val_thresh == "HIGH" & freq == 1)
ls2 <- second$MtUniqueID
sort <- c(ls, ls2)

summary$MtUniqueID <- factor (summary$MtUniqueID, levels = sort)
summary %>% 
  group_by (Stage, Identity) %>% 
  summarise(n = n_distinct(MtUniqueID))


#------------------------------------------------------------------------------# 
# Supplemental Figure 4B ----

g <- summary %>% 
# filter (Stage == "32_64") %>% # early blastocyst stage
  filter (Stage == "64_128") %>% # mid blastocyst stage
  ggplot (aes (x = MtUniqueID, y = freq, fill = factor (val_thresh)))
g <- g + geom_bar (stat = 'identity', position = 'fill')
g <- g + theme_classic() + 
  theme(aspect.ratio = 0.4,
  axis.text.x = element_blank (),
  axis.ticks = element_blank (), 
  legend.position = "top", 
  legend.title = element_blank(),  
  panel.border = element_blank(), 
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), 
  panel.spacing.x = unit (.1, "lines"),
  panel.spacing.y = unit (.1, "lines"),
  strip.background = element_blank (), 
  strip.text = element_blank ()) + 
  scale_y_continuous (labels = scales::percent) + 
  facet_wrap (Identity ~ ., nrow = 3, scales = "free") + 
  scale_fill_viridis_d (option = "A", begin = 0.3, end = 0.8, direction = -1) + 
  labs(y = "% of cell track", x = "Individual cells")
print (g)

ggsave("SuppFig4B.tiff", width = 2)
ggsave("SuppFig4B.pdf", width = 2)
