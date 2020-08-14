#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#        Cell tracks thresholded in Control vs MEKi treated embryos            #
#------------------------------------------------------------------------------#



#..............................................................................#
# Input:   Packages.R
#          MEKi_data.R
# Outputs: Supplemental Figure 2G
#          Supplemental Figure 2H
#..............................................................................#



#------------------------------------------------------------------------------# 
# Load packages, themes and data ----

setwd("~/Documents/Papers/ERKKTR/Github")
source('Packages.R')
source('MEKi_data.R')

d <- filter (embryos_CN_PD03, TE.ICM == "ICM")

good_spots <- filter (d, Value != 0 )
good_spots <- good_spots %>% 
  filter (!grepl("Mitosis", Label))
good_spots <- good_spots %>% 
  filter(!grepl("Apoptosis", Label))
good_spots <- good_spots %>% 
  filter(!grepl("Bad Segmentation", Label))

dat <- filter (good_spots, MtUniqueID %in% long_track)

dat$val_thresh <- ifelse(dat$Value > 1.2142, "HIGH", "LOW") 
# threshold from Matlab script

summary <- dat %>% 
  group_by (MtUniqueID, Embryo.Id, Treatment, val_thresh) %>% 
  summarise (n = n()) %>%   
  mutate(freq = n / sum(n))

summary <- summary %>% 
  arrange (val_thresh, desc(freq)) 
first <- subset (summary, val_thresh == "LOW")
ls <- first$MtUniqueID
second <- subset (summary, val_thresh == "HIGH" & freq == 1)
ls2 <- second$MtUniqueID
sort <- c(ls, ls2)
treat <- c("Control", "1uM PD03")

summary$MtUniqueID <- factor (summary$MtUniqueID, levels = sort)
summary$Treatment <- factor (summary$Treatment, levels = treat)

#------------------------------------------------------------------------------# 
# Supplemental Figure 2G ----

thresh <- c (LOW = "#641980", HIGH = "#FE9F6D")

g <- ggplot (data = subset (dat, MtUniqueID == "022019_p21000000019Parent"))
g <- g + geom_line (aes (group = MtUniqueID, x = TimeM, y = Value), 
                    color = "grey", size = 1) 
g <- g + geom_point (aes (x = TimeM, y = Value, color = val_thresh), alpha = 1, 
                     size = 1.5) 
g <- g + geom_hline (yintercept =  1.2142)
g <- g + theme_classic() + 
  theme(axis.line = element_line (size = 0.5),
        axis.text = element_text (size = 7), 
        legend.position="top", 
        legend.text = element_text (size = 6), 
        legend.title = element_text (size = 6),  
        strip.background = element_rect(size = 0.5),
        text = element_text (size = 8)) + 
  coord_fixed (25) +
  labs(x = "Time (min)", y = "ERK-KTR C:N")  + 
  scale_x_continuous (expand = c (0, 0)) + 
  scale_color_manual (values = thresh) 
print(g)

ggsave("SuppFig2G.tiff", width = 3)

#------------------------------------------------------------------------------# 
# Supplemental Figure 2H ----

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

ggsave("SuppFig2H.tiff", width = 2)
ggsave("SuppFig2H.pdf", width = 4)
