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
# Input:   FGF.csv
#..............................................................................#

embryos_CN_FGF <- read.csv("FGF.csv")

# Calculate track average values for each unique track ----

# Ignore when cells are undergoing mitosis or apoptosis, and poor segmentation

good_spots <- filter(embryos_CN_FGF, Value != 0 )
good_spots <- good_spots %>% 
  filter(!grepl("Mitosis", Label))
good_spots <- good_spots %>% 
  filter(!grepl("Apoptosis", Label))

# Treats parent and daughter cells the same

x <- good_spots
tracks <- unique(x$UniqueID)
average <- data.frame(UniqueID = "Empty", trackAvg = 0)
average <- average[NULL, ]

for (i in 1:length(tracks)) {
  xi <- x[x$UniqueID==tracks[i],]
  UniqueIDi <- xi$UniqueID[1]
  trackAvgi <- mean(xi$Value)
  trackVari <- var(xi$Value)
  yi <- data.frame(UniqueID = UniqueIDi, 
                   trackAvg = trackAvgi,
                   trackVar = trackVari
  )
  average <- rbind(average, yi) 
  
}

average$UniqueID<- as.character(average$UniqueID)

embryos_CN_FGF <- dplyr::left_join(embryos_CN_FGF, average, by = "UniqueID")

# Calculate track average values for each unique track ----

# Ignore when cells are undergoing mitosis or apoptosis, and poor segmentation

good_spots <- filter (embryos_CN_FGF, Value != 0 )
good_spots <- good_spots %>% 
  filter (!grepl ("Mitosis", Label))
good_spots <- good_spots %>% 
  filter (!grepl ("Apoptosis", Label))

# Calculate track avaerage treating every parent and daughter cell as unique track

x <- good_spots
Mttracks <- unique(x$MtUniqueID)
avg <- data.frame(MtUniqueID = "Empty", trackAvg = 0)
avg <- avg[NULL, ]

for (i in 1:length(Mttracks)) {
  xi <- x[x$MtUniqueID == Mttracks[i],]
  MtUniqueIDi <- xi$MtUniqueID[1]
  Identityi <- xi$Identity[1]
  TE.ICMi <- xi$TE.ICM[1]
  Treatmenti <- xi$Treatment[1]
  Embryo.Idi <- xi$Embryo.Id[1]
  trackAvgi <- mean(xi$Value)
  trackCVi <- ((sd(xi$Value)/mean(xi$Value))*100)
  yi <- data.frame(MtUniqueID = MtUniqueIDi, 
                   Identity = Identityi, 
                   TE.ICM = TE.ICMi, 
                   Treatment = Treatmenti, 
                   Embryo.Id = Embryo.Idi, 
                   trackAvg = trackAvgi, 
                   trackCV = trackCVi)
  avg <- rbind(avg, yi) 
  
}

### Long tracks for analysis ----

long_track <- count(good_spots, c("MtUniqueID"))
long_track <- filter(long_track, freq >= 5)
long_track <- dplyr::pull(long_track, MtUniqueID)