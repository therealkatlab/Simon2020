#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                         Control vs MEKi Treatments                           #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   PD03.csv
#..............................................................................#


embryos_CN_PD03 <- read.csv("PD03.csv")


# Calculate track average values for each unique track ----

# Ignore when cells are undergoing mitosis or apoptosis, and poor segmentation

good_spots <- filter(embryos_CN_PD03, Value != 0 )
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

embryos_CN_PD03 <- dplyr::left_join(embryos_CN_PD03, average, by = "UniqueID")
# Calculate track avaerage treating every parent and daughter cell as unique track

x <- good_spots
Mttracks <- unique(x$MtUniqueID)
avg <- data.frame(MtUniqueID = "Empty", trackAvg = 0)
avg <- avg[NULL, ]

for (i in 1:length(Mttracks)) {
  xi <- x[x$MtUniqueID==Mttracks[i],]
  MtUniqueIDi <- xi$MtUniqueID[1]
  Identityi <- xi$Identity[1]
  Cell_numberi <- xi$Cell_number[1]
  Stagei <- xi$Stage[1]
  TE.ICMi <- xi$TE.ICM[1]
  Treatmenti <- xi$Treatment[1]
  Embryo.Idi <- xi$Embryo.Id[1]
  trackAvgi <- mean(xi$Value)
  trackCVi <- ((sd(xi$Value)/mean(xi$Value))*100)
  yi <- data.frame(MtUniqueID = MtUniqueIDi, 
                   Identity = Identityi, 
                   Cell_number = Cell_numberi,
                   Stage = Stagei, 
                   TE.ICM = TE.ICMi, 
                   Treatment = Treatmenti, 
                   Embryo.Id = Embryo.Idi, 
                   trackAvg = trackAvgi,
                   trackCV = trackCVi)
  avg <- rbind(avg, yi) 
  
}

# Compute relative time to mitotic and apoptotic event ----

## For each mitotic track find the point of nuclear breakdown -> t = 0
## Relative mitotic time is before (-ve) or after (+ve) this event
## Input for this is manual annotation of Imaris spot

mitosis <- filter(embryos_CN_PD03, grepl("Mitosis", embryos_CN_PD03$Label) )
CN_mitosis <- dplyr::semi_join(embryos_CN_PD03, mitosis, by = "UniqueID")
x <- CN_mitosis
CN_mitosis_t <- x[NULL, ]
tracks <- unique(x$UniqueID)

for (i in 1:length(tracks)) {
  xi <- x[x$UniqueID==tracks[i],]
  yi <- filter(xi, grepl("Mitosis", xi$Label) )
  zi <- min(yi$TimeM)
  xi$MitosisT <- zi
  xi$MitosisRT <- (xi$TimeM - xi$MitosisT)
  CN_mitosis_t <- rbind(CN_mitosis_t, xi) 
  
}

### Long tracks for analysis ----

long_track <- count(good_spots, c("MtUniqueID"))
long_track <- filter(long_track, freq >= 5)
long_track <- dplyr::pull(long_track, MtUniqueID)