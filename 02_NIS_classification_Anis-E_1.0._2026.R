

## Script used to classify species into Native (NATS), NIS and Uncertain categories.
## Classification performed using ANIS-E database v.1.1 form Violet et al., 2026


# Table needed: 
# Species_occurrence_table : 
## table with the taxonomic lineage of each species and their occurrence across sites. This table must include:
##  1) one "Species_Worms_Accepted" column with all species names validated by worms  
##  2) AphiaID column
##  3) OTU_seed column
# Metadata table: 
## this table must include:
## 1) one "Ecoregion" column, reporting the correct ecoregion where samples were collected.
## 2) one Sample_Name column

## For this study the ecoregion of interest is Western Mediterraenan


# Libraries needed
library(anise)
library(dplyr)

## PREPARE DATASETS FOR NIS CLASSIFICATION ####

# Import dataset

## Label NIS from the ps object

# Import completed dataset 
combined_taxa_ps_summer.sp<-readRDS("ps_538_species_summer2021.RData")

## Extract taxa table

df<- as.data.frame(tax_table(combined_taxa_ps_summer.sp ))

colnames(df)
df$OTU_seed<-rownames(df)

# select taxonomy columns, FARTA, OTU_Seed and AphiaID
tax<-df[,c(1:7,43,45,8)]

# select OTU_seed and sample columns with counts
counts<-as.data.frame(otu_table(combined_taxa_ps_summer.sp))
occurrences<-counts
occurrences$OTU_seed<-rownames(occurrences)

# Import the meta data 

meta<- sample_data(combined_taxa_ps_summer.sp)%>% data.frame() 

# metadata with one Ecoregion column 
meta$Ecoregion<-"Western Mediterranean"

# make sure there is a column named "Sample_Name" in the metadata
colnames(meta)[colnames(meta)%in% "sample_name"]<-"Sample_Name" # replace "sample_name" with any other identifier you used for this column

# Upload ANISE datasets: intorduction, meow, origin, taxonomy, taxonomy_identifiers

data_anise() 

# combine introduciton and taxonomy tables 
intro.tax<-left_join(introduction,taxonomy, by="taxonID")

meow.df<-meow %>% select (-Lat_Zone, -WKT) # remove the last column of meow because too heavy to handle otherwise (we do not need geo polygons at this stage)

# finally combine the three tables to have a full table with the taxonID, the taxonomy and the actual ecoregion and province in which the species is introduced
intro.tax.meow<-left_join(intro.tax,as.data.frame(meow.df),by="ECO_CODE_X") 

# In the introduction table, the taxonID is a different ID from the AphiaID,
# but we need the AphiaID to make the joining of our species list and the nis table 
# so we have to use the taxonomy_identitfier table to recover the AphiaID label for each nis taxonID

taxNISAphiaID<-taxonomy_identifiers[taxonomy_identifiers$title%in%"AphiaID",]

taxNISAphiaID<- taxNISAphiaID %>% dplyr::select (taxonID,identifier)

colnames(taxNISAphiaID)[2]<-"AphiaID" 

# combine the AphiaID with the taxon intro table

intro.tax.meow<- left_join(intro.tax.meow,taxNISAphiaID, by="taxonID")

# Filter the ANISE table by your species

#AphiaID must be numeric
tax$AphiaID<-as.numeric(tax$AphiaID) # The AphiaID must be numeric column

class(tax$AphiaID)

###  EXTRACT SPECIES FROM ANISE TO CHECK THE NIS FOUND in the WHOLE DATASET

nis<-subset(intro.tax.meow, AphiaID%in% tax$AphiaID)

# how many unique nis did we found: 
length(unique(nis$AphiaID)) 

# Note that this is not the real number of NIS, 
# becasue some of those may be native of your ecoregion

# CHECK HOW MANY NIS WHERE OBSERVED IN OUR ECOREGIONS 

#Assign the syntax equal to the one of ANISe to the Ecoregions reported in the datsaet 
# reduce the function to just your ecoregion in case you only have one.

meta <- meta %>% mutate(Ecoregion=case_when(
  # Ecoregion %in%  "North_Sea" ~ "North Sea",
                                                    # Ecoregion %in%  "Aegean Sea" ~ "Aegean Sea",
                                                    # Ecoregion %in%  "South_European_Atlantic_Shelf" ~ "South European Atlantic Shelf",
                                                    # Ecoregion %in%  "Alboran Sea" ~ "Alboran Sea",
                                                    # Ecoregion %in%  "Adriatic Sea" ~ "Adriatic Sea",
                                                    # Ecoregion %in%  "Southern_Norway" ~ "Southern Norway",
                                                    # Ecoregion %in%  "Caeltic_Seas" ~ "Celtic Seas",
                                                    Ecoregion %in%  "Western Mediterranean" ~ "Western Mediterranean"))



ER.b5d<-unique(meta$Ecoregion)
 
# Check how many NIS from the ANISE db found in my dataset have been already recorded in my Ecoregion of interest
length(unique(nis$acceptedNameUsage[nis$ECOREGION %in% ER.b5d])) 

# Extract the NIS from the taxonomy table and their OTU_seed identifier

ids<-tax[tax$AphiaID %in% nis$AphiaID,c("OTU_seed","AphiaID","Species_Worms_Accepted")]

occurrence.nis<-subset(occurrences, OTU_seed %in% ids$OTU_seed)

occurrence.nis<-occurrence.nis[,c(66,1:65)]
# Check if there are samples with no nis at all 
colSums(occurrence.nis[,-1])[colSums(occurrence.nis[,-1])==0] # only bastia Mar0621_BAS_E3 

# Make long table with as many rows as many OTUs found for each sample
occ.nis.lg<-occurrence.nis %>%
  tidyr::pivot_longer(
    cols = -OTU_seed,
    names_to = "Sample_Name",
    values_to="val"
  ) %>%
  filter(val != 0)

# combine with metadata 

occ.nis.lg<- left_join(occ.nis.lg,meta, by="Sample_Name") %>% 
  dplyr::select(OTU_seed,Sample_Name,Ecoregion, Locality)

# Now combine with the tax table to obtain the AphiaID
occ.nis.lg<-left_join(occ.nis.lg, ids, by="OTU_seed") 

# Now highlight those species that are found in one ECOREGION but they are actually 
# native of that ECOREGION

# Filter the "origin" ANISe dataset to obtain the origin of each NIS 

# Obtaine the taxonID for each NIS found in the dataset:

taxID.nis<-taxNISAphiaID[taxNISAphiaID$AphiaID%in% unique(occ.nis.lg$AphiaID),]

# Filter the origin table by these identifiers
origin.sub<-origin[origin$taxonID %in% taxID.nis$taxonID,]

# combine with info from meow
origin.sub<-left_join(origin.sub,meow.df,by="ECO_CODE_X") 

# Add info about species taxon
origin.sub<-left_join(origin.sub,taxonomy,by="taxonID")

# Add AphiaID to the origin table
origin.sub<-left_join(origin.sub,taxID.nis ,by="taxonID")

# ASSIGN "NIS" - "CRYPTOGENIC" - "NATIVE" STATUS #### 

orig.check<-data.frame()

for (i in unique(occ.nis.lg$AphiaID)){ # for each AphiaID of our "putative NIS"
  
  t<-subset(occ.nis.lg, AphiaID%in% i) # obtain the rows in the introduciton (long) dataset corresponding to this APhiaID
  ER.intro<-unique(t$Ecoregion) # extract the unique Ecoregions of introduciton recorded for this NIS
  ER.orig<-unique(origin.sub$ECOREGION[origin.sub$AphiaID %in% i]) # extract the unique Ecorigions of origin for this NIS
  
  #  handle NA origin as Cryptogenic 
  
  if (all(is.na(ER.orig))) {
    
    t <- t %>% mutate(status = "Cryptogenic")  # or "Unknown origin"
    
  } else {
    
    inters<-intersect(ER.intro, ER.orig) # intersect the Ecoregions f introduciton and those of origin for this NIS
    
    
    if(all(ER.intro %in% inters)){ # if all the Ecoregions of introduction are among those intersected
      
      t <-t %>% mutate (status=case_when(t$Ecoregion %in% inters ~ "Native")) # then we can consider this NIS native in all the Ecoregions
      
    } else { # if some or none of the Ecoregions of introduction are among those intersected 
      
      t <-t %>% mutate (status=case_when(t$Ecoregion %in% inters ~ "Native", # we flag as "Native" Ecoregion, those found in the intersect
                                         !t$Ecoregion %in% inters ~ "NIS")) # we flag as NIS or Cryptogenic those that are not in the intersect 
      # - meaning that those are true NIS becasue their ECOREGION of orgining 
      # and that of introduciton does not overlap.
    }
  }
  
  orig.check<-rbind(orig.check,t) 
}

# How many Species are really NIS or Cryptogenic ?

length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "NIS"])) # 37
length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "Cryptogenic"])) #16
length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "Native"])) #21

# However this script does not assign cryptogenic label to those species that are cryptogenic in some locaitons but not in others,
# these species do not have NA in the origin table. But they have the label "cryptogenic" in the introduciton table.

putative.crypto<-intro.tax.meow[intro.tax.meow$AphiaID%in%orig.check$AphiaID[orig.check$status%in%"NIS"],]

# filter only those mentions in the West Med Sea
putative.crypto.Med <- putative.crypto[putative.crypto$ECOREGION%in% "Western Mediterranean",]

# select the AphiaID that have the label cryptogenic in the establishmentmeans column at least one time
# if in the Country column there is France --> assign the establishment mean of France to that NIS 
# if in the Country column there is no France, then assign Cryptogenic

crypto.Med.real<-putative.crypto.Med %>% group_by(AphiaID) %>% filter(any(establishmentMeans=="Cryptogenic"))

for (i in unique(crypto.Med.real$AphiaID)){
  
  if ( any(crypto.Med.real$country[crypto.Med.real$AphiaID%in% i] =="France")){
  
      orig.check$status[orig.check$AphiaID%in% i]<-crypto.Med.real$establishmentMeans[crypto.Med.real$AphiaID %in% i & crypto.Med.real$country%in%"France"]
  } else {
    
    orig.check$status[orig.check$AphiaID%in% i]<-unique(crypto.Med.real$establishmentMeans[crypto.Med.real$AphiaID %in% i])
  }
}



orig.check$status[orig.check$status%in%"Non-indigenous"]<-"NIS"


# How many Species are really NIS or Cryptogenic ?

length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "NIS"])) # 35
length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "Cryptogenic"])) #18
length(unique(orig.check$Species_Worms_Accepted[orig.check$status%in% "Native"])) #21


# how many NIS found in each Ecoregion (if you have more than one)
orig.check %>% filter(status %in% "NIS") %>%
  group_by(Ecoregion) %>% 
  summarise(unique_sp=length(unique(Species_Worms_Accepted)),
            unique_otu=length(unique(OTU_seed)))

orig.check %>% filter(status %in% "Cryptogenic") %>%
  group_by(Ecoregion) %>% 
  summarise(unique_sp=length(unique(Species_Worms_Accepted)),
            unique_otu=length(unique(OTU_seed)))


# We are working on samples from the Western Mediterranean Sea, however, any species native from the larger Mediterraenan province should be considered native to Mediterraenan
# Check if any of the species classified as NIS in the list are actually native from other Ecoregions of the Mediterraena sea but not from "Western Mediterraenan"

other.Med.ER<-unique(meow.df$ECOREGION[meow.df$PROVINCE%in%"Mediterranean Sea"])

aphia.nis<-orig.check$AphiaID[orig.check$status%in% "NIS"]

#see if they intersect
intersect( other.Med.ER,unique(origin.sub$ECOREGION[origin.sub$AphiaID %in% aphia.nis]))
# None 

#### IMPORT THE NEW STATUS IN THE ORIGINAL TABLE ####

# select only nis and cryotogenic

orig.check.nis<-orig.check %>% subset(status %in% c("NIS" ))
orig.check.cryp<-orig.check %>% subset(status %in% c("Cryptogenic" ))

ll.nis<-orig.check.nis %>% group_by(Ecoregion) %>% summarise(ids.nis=list(unique(OTU_seed)))
ll.cryp<-orig.check.cryp %>% group_by(Ecoregion) %>% summarise(ids.nis=list(unique(OTU_seed)))

## Rename columns of NIS assignation made with old ANIS-E dataset

colnames(df)


# assign the status of NIS or Cryptogenic to the OTUs taxa table according the Aphia ID 

df.NIS <- df %>%
  mutate(
    NIS_Western_Mediterranean = case_when(
      OTU_seed %in% ll.nis$ids.nis[ll.nis$Ecoregion == "Western Mediterranean"][[1]] ~ "NIS",
      OTU_seed %in% ll.cryp$ids.nis[ll.cryp$Ecoregion == "Western Mediterranean"][[1]] ~ "Cryptogenic",
      TRUE ~ NA
       )
    )

# put the NIS status before the sample columns
colnames(df.NIS)
# df.NIS<-df.NIS[,c(1:24,92,27:91)]


#### ADD EUROPEAN ORIGIN COLUMN ####
### Create a column where for each NIS with known origin in Europe we report a list of the European Ecoregions in which this NIS is native
# This will be useful if we want to analyse the genetic diversity of a NIS in its native range in Europe and its invaded range in Europe

# extact the nis fro the df.NIS datase: 
df.NIS$AphiaID<-as.numeric(df.NIS$AphiaID)

n<-df.NIS$AphiaID[df.NIS$NIS_Western_Mediterranean %in% "NIS"]

# obtain their orgin from the origin dataset 
n.orig<-origin.sub[origin.sub$AphiaID %in% n,]

# extract the NIS with EU origin
list.true.european.nis<-n.orig[n.orig$nativeEurope=="TRUE",c("AphiaID","ECOREGION")]

# create a variable with the list of the origin European ecoregions for each NIS
colnames(list.true.european.nis)[2]<-"ECOREGION_EU_ORIGIN"

list.true.european.nis<-list.true.european.nis %>% group_by(AphiaID) %>% 
  summarise(NIS_EU_ORIGIN=paste0(unique(ECOREGION_EU_ORIGIN), collapse=", "))

# Add this column to the dataset 

df.NIS<-left_join(df.NIS,list.true.european.nis,by="AphiaID")

rownames(df.NIS)<-df.NIS$OTU_seed

## Save a table where for each nis and cryptogeneic species found you report the list of introdrocution sites and the origin

nis.crypt.final<-subset(intro.tax.meow, AphiaID%in% df.NIS$AphiaID[df.NIS$NIS_Western_Mediterranean%in%c("NIS","Cryptogenic")])

length(unique(nis.crypt.final$AphiaID))

colnames(nis.crypt.final)

nis.crypt.final<-nis.crypt.final[,c("AphiaID","acceptedNameUsage","year","country","ECOREGION","PROVINCE","REALM",
                                    "establishmentMeans",
                                    "degreeOfEstablishment",
                                   "occurrenceRemarks","associatedReferences","references")]

# make list of NIS that have not been recorded yet in the Mediterraenan sea
nis.crypt.final<-as.data.frame(nis.crypt.final )

alreadyObservedMedSea<-unique(nis.crypt.final$acceptedNameUsage[nis.crypt.final$PROVINCE%in% "Mediterranean Sea"])

newEU_nis<-unique(nis.crypt.final$acceptedNameUsage[!nis.crypt.final$acceptedNameUsage %in% alreadyObservedMedSea])

# same for the orgin table

nis.crypt.origin<-subset(origin.sub, AphiaID%in% df.NIS$AphiaID[df.NIS$NIS_Western_Mediterranean%in%c("NIS","Cryptogenic")])
colnames(nis.crypt.origin)

nis.crypt.origin<-nis.crypt.origin[,c("AphiaID","acceptedNameUsage", "nativeEurope",
                                     "ECOREGION","PROVINCE","REALM", "references" )]

# write.csv(nis.crypt.final,file="Final_files_pubblication/NIS_Crypto_introduction_ranges.csv")
# write.csv(nis.crypt.origin,file="Final_files_pubblication/NIS_Crypto_native_ranges.csv")

# MANUAL SPECIES CHECK ####

### To the automatic classification performed using ANIS-E a manual one, followed 
### to check for those species that are not easily distinguishable form native one based 
### on the markers considered in the study

# Ciona intestinalis and Mytilus edulis 
# The former is only identified by 18S marker and this is not enough to distinguish it from C. robusta. Excluded from the analyses 
# The latter is only identified by 18S marker and this is not enough to distinguish it from M. galloprovincialis and M. trossolus. Excluded from the analyses.

df.NIS<-subset(df.NIS, !FARTA %in% c("Ciona intestinalis",
                                     "Mytilus edulis"))


# According with the discussion had with Argyro Zenetos we should keep Phallusia mammillata as native and Jassa slatteryi as NIS-MED 

df.NIS<-df.NIS %>% mutate( STATUS = case_when(FARTA %in% newEU_nis & NIS_Western_Mediterranean %in% "NIS" ~ "NIS-EU",
                                                     !FARTA %in% newEU_nis & NIS_Western_Mediterranean%in% "NIS"~"NIS-MED"))

df.NIS$STATUS[df.NIS$NIS_Western_Mediterranean%in% "Cryptogenic"] <- "Cryptogenic"


df.NIS$STATUS[df.NIS$FARTA %in% "Phallusia mammillata"]<-"Native from Med Sea (and not introduced elsewhere in Eur)"
df.NIS$STATUS[df.NIS$FARTA %in% "Jassa slatteryi"]<-"NIS-MED"

# Four species are classified as NIS-EU by ANIS-E classification, however this is a mistake because in the old version of ANISE they were classified as introdocued in the Mediterraenanen Sea

# Eriocheir sinensis
# Watersipora subatra
# Aurelia coerulea
# Pseudodiaptomus marinus

df.NIS$STATUS[df.NIS$FARTA %in% "Eriocheir sinensis"]<-"NIS-MED"
df.NIS$STATUS[df.NIS$FARTA %in% "Watersipora subatra"]<-"NIS-MED"
df.NIS$STATUS[df.NIS$FARTA %in% "Aurelia coerulea"]<-"NIS-MED"
df.NIS$STATUS[df.NIS$FARTA %in% "Pseudodiaptomus marinus"]<-"NIS-MED"

setdiff(sort(df.NIS$FARTA[df.NIS$STATUS%in% "NIS-EU"]),
        sort(fred.table$Species[fred.table$`Updated_Status_TO USE BY ETIENNE`%in% c("NIS-EU (ANIS-E - version 1.0)" ,
                                                                                    "NIS-EU")]))
length(df.NIS$Species[df.NIS$STATUS%in% c("NIS-MED","NIS-EU","Cryptogenic")]) # 50 in total
length(df.NIS$Species[df.NIS$STATUS%in% c("NIS-MED")]) # 50 in total
length(df.NIS$Species[df.NIS$STATUS%in% c("Cryptogenic")]) # 50 in total



write.csv(df.NIS,file="Final_files_pubblication/538_species_NIS_classified_ANIS-E_version_1.0.csv")

# Save new phyloseq object with classification


## PROCEDE WITH MANUAL ASSIGNATION OF NATIVE SPECIES TO MEDITERRAENAN SEA #######
# Reload the taxonomy table with annotated with labels for species non classified by ANIS-E as NIS or Cryptogenic
# Native: species that are found in litterature with a distribution rangeincludin galso Mediterranenan Sea, or cosmopolitan species;
# Uncertain: species that in litterature do not show "Mediterraenan Sea" in their distribution range; so potetial new NIS, but not native from Med Sea.

new.class.tax<-read.csv("Final_files_pubblication/536_species_NIS_classified_ANIS-E_version_1.0.csv") #file with new classificaiton of Native and Uncertain
rownames(new.class.tax)<-new.class.tax$OTU_seed

# Remove Ciona intestinalis and Mytilus edulis from the whole phyloseq object 

combined_taxa_ps_summer.sp.nociona.mytilus<-subset_taxa(combined_taxa_ps_summer.sp, !FARTA %in% c("Ciona intestinalis",
                                                                                                  "Mytilus edulis"))

### reload and save in the 

# Import in the phyloseq object
tax_table(combined_taxa_ps_summer.sp)<-as.matrix(new.class.tax)

saveRDS(combined_taxa_ps_summer.sp,"~ps_536_species_summer2021_new_NIS_classification_ANISE1.1.RData")
