#------------------------------------------------------------------------------#
#  This script accompanies the manuscript                                      #
#  Simon et al., (2020) Developmental Cell                                     #
#  Repository available on https://github.com/therealkatlab                    #
#  Please consult READ_ME for more information                                 #
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#
#                           Lineages in 12hr movies                           #
#------------------------------------------------------------------------------#


#..............................................................................#
# Input:   Lineage_12hr.csv.csv
#..............................................................................#

## Read in combined data from IMARIS and MINS


embryos_CN <- read.csv("Lineage_12hr.csv")


# Get identities based on final time point ----

final <- filter (embryos_CN, (Mitosis == "Daughter1" | Mitosis == "Daughter2") 
  & Time == 48)
variables <- c("UniqueID", "MtUniqueID", "Identity")
final_ids <- final[variables]

x <- final_ids %>% 
  group_by (UniqueID) %>% 
  distinct (Identity)
tracks <- unique (x$UniqueID)
endid <- data.frame (UniqueID = "Empty", End_Identity = "Empty")
endid <- endid[NULL, ]

for (i in 1:length(tracks)) {
  xi <- x[x$UniqueID == tracks[i],]
  UniqueIDi <- xi$UniqueID[1]
  y <- as.factor (xi$Identity)
  if (length(y) == 1) {
    z <- paste0(y, sep=";", y) }
  else {
    z <- paste0(y[[1]], sep=";", y[[2]]) }
  a <- data.frame(UniqueID = UniqueIDi, End_Identity = z)
  endid <- rbind(endid, a)
}

# Generate the final identities of cells and their progeny

parents <- filter(embryos_CN, Mitosis == "Parent")
parents <- as.factor(parents$UniqueID)
endid_parent <- subset(endid, UniqueID %in% parents)
endid_parent$MtUniqueID <- paste0(endid_parent$UniqueID, sep = "_", "Parent")
endid_parent <- endid_parent[-1]

cell_divide <- embryos_CN %>% 
  group_by (UniqueID) %>% 
  filter ( (Mitosis == "Daughter1" | Mitosis == "Daughter2") & Time == 48) %>% 
  summarise(n = n())
both <- cell_divide %>% 
  filter(n == 2)
both.Id <- both$UniqueID
one <- cell_divide %>% 
  filter(n == 1)
one.Id <- paste0(one$UniqueID, sep = "_", "Parent")

endid_merge <- endid_parent

endid_merge$End_Identity <- as.character(endid_merge$End_Identity)
endid_merge1 <- endid_merge %>% 
  filter(MtUniqueID %in% one.Id) %>% 
  mutate(End_Identity = paste0(End_Identity, sep = ";", "ICM"))
endid_merge2 <- endid_merge %>% 
  filter(!(MtUniqueID %in% one.Id)) %>% 
  mutate(End_Identity = End_Identity)
endid_merge <- rbind(endid_merge1, endid_merge2)

endid_merge$End_Identity[endid_merge$End_Identity == "DP;DP"] <- "DP&DP"
endid_merge$End_Identity[endid_merge$End_Identity == "EPI;EPI"] <- "EPI&EPI"
endid_merge$End_Identity[endid_merge$End_Identity == "PrE;PrE"] <- "PrE&PrE"

endid_merge$End_Identity[endid_merge$End_Identity == "DP;DP;ICM"] <- "DP;ICM"
endid_merge$End_Identity[endid_merge$End_Identity == "EPI;EPI;ICM"] <- "EPI;ICM"
endid_merge$End_Identity[endid_merge$End_Identity == "PrE;PrE;ICM"] <- "PrE;ICM"

endid_merge$End_Identity[endid_merge$End_Identity == "DP;EPI"] <- "DP&EPI"
endid_merge$End_Identity[endid_merge$End_Identity == "EPI;DP"] <- "DP&EPI"
endid_merge$End_Identity[endid_merge$End_Identity == "EPI;PrE"] <- "EPI&PrE"
endid_merge$End_Identity[endid_merge$End_Identity == "PrE;EPI"] <- "EPI&PrE"
endid_merge$End_Identity[endid_merge$End_Identity == "PrE;DP"] <- "DP&PrE"
endid_merge$End_Identity[endid_merge$End_Identity == "DP;PrE"] <- "DP&PrE"

embryos_CN <- dplyr::left_join(embryos_CN, endid_merge, by = "MtUniqueID")
embryos_CN$Identity2 <- as.character(embryos_CN$Identity)
embryos_CN$End_Identity[is.na(embryos_CN$End_Identity)] <- as.character(embryos_CN$Identity2[is.na(embryos_CN$End_Identity)])
embryos_CN$Identity2 <- NULL
embryos_CN <- arrange(embryos_CN, MtUniqueID, Time)


# Calculate track average values for each unique track ----

# Ignore when cells are undergoing mitosis or apoptosis, and poor segmentation

good_spots <- filter(embryos_CN, Value != 0 )
good_spots <- good_spots %>% 
  filter(!grepl("Mitosis", Label))
good_spots <- good_spots %>% 
  filter(!grepl("Apoptosis", Label))
good_spots <- good_spots %>% 
  filter(!grepl("Bad Segmentation", Label))

# Calculate track averages when treating parent and daughter cells the same unit
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

embryos_CN <- dplyr::left_join(embryos_CN, average, by = "UniqueID")

# Calculate track avaerages treating every parent and daughter cell as a unique unit

x <- good_spots
Mttracks <- unique(x$MtUniqueID)
avg <- data.frame(MtUniqueID = "Empty", trackAvg = 0)
avg <- avg[NULL, ]

for (i in 1:length(Mttracks)) {
  xi <- x[x$MtUniqueID==Mttracks[i],]
  MtUniqueIDi <- xi$MtUniqueID[1]
  Identityi <- xi$Identity[1]
  End_Identityi <- xi$End_Identity[1]
  Cell_numberi <- xi$Cell_number[1]
  Stagei <- xi$Stage[1]
  TE.ICMi <- xi$TE.ICM[1]
  Embryo.Idi <- xi$Embryo.Id[1]
  trackAvgi <- mean(xi$Value)
  trackCVi <- ((sd(xi$Value)/mean(xi$Value))*100)
  yi <- data.frame(MtUniqueID = MtUniqueIDi, 
                   Identity = Identityi, 
                   End_Identity = End_Identityi, 
                  Cell_number = Cell_numberi,
                  Stage = Stagei, 
                   TE.ICM = TE.ICMi, 
                   Embryo.Id = Embryo.Idi, 
                   trackAvg = trackAvgi, 
                   trackCV = trackCVi)
  avg <- rbind(avg, yi) 
  
}

# Compute relative time to mitotic and apoptotic event ----

## For each mitotic track find the point of nuclear breakdown -> t = 0
## Relative mitotic time is before (-ve) or after (+ve) this event
## Input for this is manual annotation of Imaris spot

mitosis <- filter(embryos_CN, grepl("Mitosis", embryos_CN$Label) )
CN_mitosis <- dplyr::semi_join(embryos_CN, mitosis, by = "UniqueID")
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