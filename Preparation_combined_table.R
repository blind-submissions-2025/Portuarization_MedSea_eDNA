# Processing of the OTUs tables generated for each marker (12S, 18S, COI) to:
# 1) Rarefying each OTU table to the minimum numeber of reads
# 2) Selecting only OTU assigned to Animalia
# 3) Remove all non-marine organisms
# 4) Create the non-redundant combined taxonomic dataset: joining the three filtered and rarefied OTU table into a final taxonomic table.
# 5) Subset the dataset, keeping only taxa classified down to species level


# Load libraries 
pkgs <- c("phyloseq","dplyr","tidyr","broom","worrms","tibble","readxl","stringr")
lapply(pkgs, require,character.only=TRUE)

# Load OTU tables
data<-load("med_int_ext_allBLAST.RData")
names(otu_phyloseq_objects) <- c("12S", "18S", "COI")

# Function to rename taxa in a phyloseq object
rename_phyloseq <- function(physeq_obj, prefix_replace = "ASV_", prefix_new = "OTU_") {
  # Get current taxa names
  taxa_names(physeq_obj) <- gsub(prefix_replace, prefix_new, taxa_names(physeq_obj))
  
  # Update OTU table
  if (!is.null(otu_table(physeq_obj))) {
    taxa_names(otu_table(physeq_obj)) <- taxa_names(physeq_obj)
  }
  
  # Update taxonomy table
  if (!is.null(tax_table(physeq_obj))) {
    taxa_names(tax_table(physeq_obj)) <- taxa_names(physeq_obj)
  }
  
  return(physeq_obj)
}

# Apply renaming function to all phyloseq objects for 12S, 18S, and COI
otu_phyloseq_objects <- lapply(otu_phyloseq_objects, rename_phyloseq)

modify.meta<-function(phy){
  meta<-sample_data(phy)%>% data.frame()
  meta$collection_date<-as.Date(meta$collection_date, "%d/%m/%Y")
  meta<-meta %>% mutate(Season=case_when(collection_date >= "2021-03-17" & collection_date < "2021-06-15" ~"Spring",
                                         collection_date >= "2021-06-15" & collection_date <= "2021-07-08" ~ "Summer"))
  
  meta<-meta %>% mutate(Type_habitat_broad = case_when(Type_Habitat == "Natural"~"Natural",
                                                       Type_Habitat == "Natural_OutsidePort"~"Natural",
                                                       Type_Habitat == "Marina"~"Artificial",
                                                       Type_Habitat == "Harbour"~"Artificial"))
  
  
  meta<-meta %>% mutate(Region = case_when(Locality == "Sete" ~ "Continent",
                                           Locality == "Seyne-sur-mer"  ~ "Continent",
                                           Locality == "Marseille" ~ "Continent",
                                           Locality == "Nice" ~ "Continent",
                                           Locality == "Port-la-Nouvelle" ~ "Continent",
                                           Locality == "Toulon" ~ "Continent", 
                                           Locality == "Ajaccio" ~ "Corsica",
                                           Locality == "Saint-Tropez"~ "Continent",
                                           Locality == "Bastia"~ "Corsica",
                                           Locality == "Bonifaccio" ~ "Corsica",
                                           Locality == "Calvi"~ "Corsica",
                                           Locality == "Porto-Vecchio"~ "Corsica"))
  sample_data(phy)<-sample_data(meta)
  return(phy)
}

otu_phyloseq_objects <- lapply(otu_phyloseq_objects, modify.meta)

# Select only summer samples

summer.select<-function(ps){
  subset<-subset_samples(ps, Season %in% "Summer")
  return(subset)
}

otu_phyloseq_objects_summer<-lapply(otu_phyloseq_objects, summer.select)

# Function to replace Chordata with Vertebrata or Tunicata in the Phylum column, based on Class assignments

replace_chordata_phylum <- function(tax_table_df) {
  # Copy the original Phylum column to the col "initial_phylum_assignment"
  tax_table_df <- tax_table_df %>%
    mutate(initial_phylum_assignment = Phylum) %>%
    
    # Modify the Phylum column based on Class values
    mutate(Phylum = case_when(
      Class %in% c("Actinopteri", "Mammalia", "Aves", "Amphibia", "Chondrichthyes", "Lepidosauria") ~ "Vertebrata",
      Class %in% c("Ascidiacea", "Thaliacea", "Appendicularia") ~ "Tunicata",
      Class %in% c("Leptocardii") ~ "Cephalochordata",
      Phylum == "Chordata" & is.na(Class) ~ "Unclassified-Chordata",
      Phylum == "Chordata"~ "Vertebrata",# For any remaining Chordata, set to Vertebrata
      TRUE ~ Phylum  # Keep original Phylum value if no condition matches
    ))
  
  return(tax_table_df)
}

# Apply the replacement to each ASV and OTU phyloseq object

for (i in seq_along(otu_phyloseq_objects_summer)) {
  tax_table_df <- data.frame(tax_table(otu_phyloseq_objects_summer[[i]]))  # Extract tax table as data frame
  tax_table_df <- replace_chordata_phylum(tax_table_df)  # Replace Chordata in Phylum column
  tax_table(otu_phyloseq_objects_summer[[i]]) <- as.matrix(tax_table_df)  # Update tax table in phyloseq object
}


# Check results for one phyloseq object to confirm replacements

phyla_list_combined <- unique(
  unlist(lapply(c(otu_phyloseq_objects_summer), function(physeq) {
    # Extract the tax_table, convert it to a data frame, and select the Phylum column
    as.data.frame(tax_table(physeq))$Phylum
  }))
)

# CHUNK 1: Rarefaction ----
# Function to rarefy phyloseq objects to 100000 reads unless the minimum number of reads is higher 
rarefy_per_marker <- function(phyloseq_list) {
  lapply(phyloseq_list, function(ps) {
    
    # Determine the smallest sample size
    min_reads <- min(sample_sums(ps))
    if (min_reads>= 1e+05){
      cat(paste("rarefied to minimum number of reads"))
      
      # Rarefy the phyloseq object to the minimum if this is >= 100000
      rarefied_ps <- rarefy_even_depth(ps, 
                                       sample.size = min_reads, 
                                       rngseed = 123,  # Set seed for reproducibility
                                       verbose = TRUE)}
    else{
      cat(paste("rarefied to 1e+05 reads"))
      rarefied_ps <- rarefy_even_depth(ps, 
                                       sample.size = 1e+05, 
                                       rngseed = 123,  # Set seed for reproducibility
                                       verbose = TRUE)
    }
    return(rarefied_ps)
  })
}

# Rarefy OTU phyloseq objects
set.seed(1234)
rarefied_otu_summer_objects <- rarefy_per_marker(otu_phyloseq_objects_summer)

ntaxa(rarefied_otu_summer_objects$`12S`)*100/ntaxa(prune_taxa(taxa_sums(otu_phyloseq_objects_summer$`12S`)!=0,otu_phyloseq_objects_summer$`12S`))
ntaxa(rarefied_otu_summer_objects$`18S`)*100/ntaxa(prune_taxa(taxa_sums(otu_phyloseq_objects_summer$`18S`)!=0,otu_phyloseq_objects_summer$`18S`))
ntaxa(rarefied_otu_summer_objects$`COI`)*100/ntaxa(prune_taxa(taxa_sums(otu_phyloseq_objects_summer$`COI`)!=0,otu_phyloseq_objects_summer$`COI`))

# Report the same samples in each dataset
# I do this in case the threshold to rarefaction was not the minimum number of reads
# Different samples were removed from each dataset but at the end the samples have to be the same in the three datasets.

# Function to list samples below a threshold
list_samples_below_threshold <- function(phyloseq_list, threshold = 1e+05) {
  lapply(phyloseq_list, function(ps) {
    reads_per_sample <- sample_sums(ps)
    samples_below <- names(reads_per_sample[reads_per_sample < threshold])
    samples_below
  })
}

# Function to list OTU and ASV lost due to rarefaction
list_otu_asv_lost <- function(ps,rare.ps) {
  rare.kept<-taxa_names(rare.ps)
  all<-taxa_names(ps)
  rare.left.out<-setdiff(all,rare.kept)
  rare.left.out
}

# List samples below 1E+05 for OTU phyloseq objects
otu_samples_below_threshold_summer <- list_samples_below_threshold(otu_phyloseq_objects_summer, threshold = 1e+05)

# List OTU left out after rarefaction
otu_rare_left_out_12S_summer<- list_otu_asv_lost(otu_phyloseq_objects_summer$`12S` ,rarefied_otu_summer_objects$`12S`)

otu_rare_left_out_18S_summer<- list_otu_asv_lost(otu_phyloseq_objects_summer$`18S` ,rarefied_otu_summer_objects$`18S`)

otu_rare_left_out_COI_summer<- list_otu_asv_lost(otu_phyloseq_objects_summer$COI ,rarefied_otu_summer_objects$COI)

# Print results

print("Samples Below 1E+05 Reads for OTU Objects:")
print(otu_samples_below_threshold_summer)

# Remove samples below threshold from all the rarefied objects

# On OTU objects
remove_below_threshold<-function(phyloseq_list) {
  
  lapply(phyloseq_list, function(rarefied.ps){
    
    rarefied.ps<-subset_samples(rarefied.ps,!sample_name %in% unname(unlist(otu_samples_below_threshold_summer)))
    
    rarefied.ps<-prune_taxa(taxa_sums(rarefied.ps)!=0, rarefied.ps)
    
    return(rarefied.ps)
  })
}


rarefied_otu_summer_objects<-remove_below_threshold(rarefied_otu_summer_objects)

summary.samples<-function(ps){
  meta<-sample_data(ps)%>% data.frame()
  
  meta.summer<-subset(meta, Season %in% "Summer") # Check only summer since that is the season you are interested in 
  
  summary<- meta.summer %>% group_by(Locality, Type_habitat_broad) %>% summarise(numb.samples=n())
  
  return(summary)
}

# Pre-rarefaction

summary.samples.pre.rarefaciton.otu<-lapply(otu_phyloseq_objects_summer, summary.samples)

# Post rarefaction

summary.samples.post.rarefaciton.otu<-lapply(rarefied_otu_summer_objects, summary.samples)

# Function to count OTUs before and after rarefaction
count_asvs <- function(initial_list, rarefied_list) {
  lapply(seq_along(initial_list), function(i) {
    initial_count <- ntaxa(initial_list[[i]])  # Number of ASVs/OTUs before rarefaction
    retained_count <- ntaxa(rarefied_list[[i]])  # Number of ASVs/OTUs after rarefaction
    list(
      marker = names(initial_list)[i],  # Marker name
      initial = initial_count,
      retained = retained_count
    )
  })
}


# Calculate counts for OTU phyloseq objects
otu_counts <- count_asvs(otu_phyloseq_objects_summer, rarefied_otu_summer_objects)

# Print results for OTU
cat("\nOTU Counts Before and After Rarefaction:\n")
print(otu_counts)

## CHUNK 2: Explore then prepare Animalia data ----

# Obtain all phyla in the taxonomy tables of the rarefied objects

library(worrms)
# Produce a list of all phyla and check those  

phyla_list_combined.rarefied.otu <- unique(
  unlist(lapply(c(rarefied_otu_summer_objects), function(physeq) {
    # Extract the tax_table, convert it to a data frame, and select the Phylum column
    as.data.frame(tax_table(physeq))$Phylum
  }))
)

animalia.list.worms<-wm_children_(name= "Animalia", marine_only = FALSE)

# create new column in the taxonomy table to assign the Kingdom "Animalia"

assign.animalia<-function(phyloseq_list){
  
  lapply(phyloseq_list, function (ps){
    
    tax.table<-as.data.frame(tax_table(ps))
    
    tax.table<-tax.table %>% mutate(Category =
                                      case_when(is.na(initial_phylum_assignment)~"Unassigned",
                                                initial_phylum_assignment %in% animalia.list.worms$scientificname ~"Animalia",
                                                !initial_phylum_assignment %in% animalia.list.worms$scientificname ~"Other"))      
    
    tax_table(ps)<-tax_table(as.matrix(tax.table))
    return(ps)
  })
}

rare_otu_phyloseq_objects.animalia<-assign.animalia(rarefied_otu_summer_objects)

print("12S")
length(taxa_names(subset_taxa(rare_otu_phyloseq_objects.animalia$`12S`, Category%in%"Animalia")))
print("18S")
length(taxa_names(subset_taxa(rare_otu_phyloseq_objects.animalia$`18S`, Category%in%"Animalia")))
print("COI")
length(taxa_names(subset_taxa(rare_otu_phyloseq_objects.animalia$COI, Category%in%"Animalia")))

# Check all the non animalia phyla to be sure you haven't skipped some

list.non.animalia<-function(phyloseq_list){
  unique(unlist(unname(lapply(phyloseq_list, function (ps){
    
    tax.table<-as.data.frame(tax_table(ps))%>% subset(Category == "Other")
    
    non.animalia.phyla<-unique(tax.table$initial_phylum_assignment)
    
  }))))
}


list.non.animalia.otu<-list.non.animalia(rare_otu_phyloseq_objects.animalia)

# Obtain the list of OTUs not included in the animalia kingdom
list.non.animalia.ASVs.OTUs<-function(phyloseq_list){
  
  lapply(phyloseq_list, function (ps){
    
    list.non.animalia<- subset_taxa(ps, !Category %in% "Animalia") %>% 
      
      taxa_names()
    
  })
}

non.Animalia.OTUs<-list.non.animalia.ASVs.OTUs(rare_otu_phyloseq_objects.animalia)

#Manually check if there are some Animalia Phyla that have been left

cat("Phyla not included in the animalia Kingdom in WORMS")

print(list.non.animalia.otu)

# Filter to keep only Animalia OTUs ----

filter_phyloseq_for_animalia<-function (phyloseq_list){
  
  lapply(phyloseq_list, function (ps){
    
    ps.animalia<-subset_taxa(ps, Category %in% "Animalia")
    
    tax.table<- as.data.frame(tax_table(ps.animalia))
    
    tax.table<- tax.table %>% mutate(initial_Kingdom_assignation=Kingdom) %>%
      mutate(Kingdom=Category)
    
    tax_table(ps.animalia) <- tax_table(as.matrix(tax.table))
    return(ps.animalia)
    
  }
  
  )
}

rare_otu_phyloseq_objects.animalia.bf.filtering<-rare_otu_phyloseq_objects.animalia

# Apply filtering to all otu_phyloseq objects
rare_otu_phyloseq_objects.animalia <- filter_phyloseq_for_animalia(rare_otu_phyloseq_objects.animalia)

# Print final number of OTUs assogned to Animalia
dim(tax_table(rare_otu_phyloseq_objects.animalia$`12S`))[1] # 562
dim(tax_table(rare_otu_phyloseq_objects.animalia$`18S`))[1] # 376
dim(tax_table(rare_otu_phyloseq_objects.animalia$`COI`))[1] # 422

# CHUNK 3: Remove non marine Insecta, Mammalia and Birds ----
# Function to recognize non marine taxa in general, with a special focus on Insecta and Mammalia for this dataset

marine.org <- function(ps, organism, taxlevel) {
  # Extract taxonomy and counts
  tax <- as.data.frame(tax_table(ps))
  counts <- as.data.frame(otu_table(ps))
  
  # Filter for specified organisms
  if (!is.null(organism)) {
    tax <- subset(tax, Class %in% organism)
  }
  
  # Query WoRMS for marine taxa
  valid_names <- unique(tax[, taxlevel])
  valid_names <- valid_names[!is.na(valid_names) & valid_names != ""]
  
  marine_list <- tryCatch({
    wm_records_names(name = valid_names, marine_only = TRUE)
  }, error = function(e) {
    warning("WoRMS query failed: ", conditionMessage(e))
    return(list())
  })
  
  marine_list <- purrr::compact(marine_list)  # Remove NULL or empty elements
  
  # Extract scientific names of marine taxa
  retained_names <- if (length(marine_list) > 0) {
    unlist(lapply(marine_list, function(x) x$scientificname))
  } else {
    character()
  }
  
  # Retain ASVs matching the retained names
  retained_asvs <- rownames(tax)[tax[, taxlevel] %in% retained_names]
  retained_table <- tax[retained_asvs, ]
  
  # Identify and save removed ASVs
  removed_asvs <- setdiff(rownames(tax), retained_asvs)
  removed_table <- cbind(tax[removed_asvs, ], counts[removed_asvs, ])
  
  return(list(
    retained_asvs = retained_asvs,
    retained_table = retained_table,
    removed_asvs = removed_asvs,
    removed_table = removed_table
  ))
}



# Initialize empty data frames for results
list_OTU_Mammalia_removed <- data.frame()
list_OTU_Insecta_removed <- data.frame()
list_OTU_Aves_removed <- data.frame()
list_OTU_Arachnida_removed <- data.frame()

# Loop through the otu_phyloseq_filtered_objects
for (i in seq_along(rare_otu_phyloseq_objects.animalia)) { #seq_along(rare_otu_phyloseq_objects.animalia)
  ps <- rare_otu_phyloseq_objects.animalia[[i]]
  marker <- paste0("Marker_", names(rare_otu_phyloseq_objects.animalia[i]))  # Replace with actual marker names if available
  
  # Apply marine.org to Mammalia
  result_mammalia <- marine.org(ps, organism = "Mammalia", taxlevel = "Species")
  if (length(result_mammalia$removed_asvs) > 0) {  # Note: renamed to "removed_otus" for OTUs
    to_remove_mammalia <- data.frame(marker = marker, OTU = result_mammalia$removed_asvs, result_mammalia$removed_table)
    list_OTU_Mammalia_removed <- rbind(list_OTU_Mammalia_removed, to_remove_mammalia)
  }
  
  # Apply marine.org to Insecta
  result_insecta <- marine.org(ps, organism = "Insecta", taxlevel = "Species")
  if (length(result_insecta$removed_asvs) > 0) {
    to_remove_insecta <- data.frame(marker = marker, OTU = result_insecta$removed_asvs, result_insecta$removed_table)
    list_OTU_Insecta_removed <- rbind(list_OTU_Insecta_removed, to_remove_insecta)
  }
  
  # Apply marine.org to Arachnida
  result_arachnida <- marine.org(ps, organism = "Arachnida", taxlevel = "Species")
  if (length(result_arachnida$removed_asvs) > 0) {
    to_remove_arachnida <- data.frame(marker = marker, OTU = result_arachnida$removed_asvs, result_arachnida$removed_table)
    list_OTU_Arachnida_removed <- rbind(list_OTU_Arachnida_removed, to_remove_arachnida)
  }
  
  # Apply marine.org to Birds
  result_aves <- marine.org(ps, organism = "Aves", taxlevel = "Species")
  if (length(result_aves$removed_asvs) > 0) {
    to_remove_aves <- data.frame(marker = marker, OTU = result_aves$removed_asvs, result_aves$removed_table)
    list_OTU_Aves_removed <- rbind(list_OTU_Aves_removed, to_remove_aves)
  }
  
}


#remove non marine insecta, arachnida, mammalia and birds
# Functions to prune ASVs and OTUs from a phyloseq object

remove_otus_from_list <- function(ps, otus_to_remove) {
  pruned_ps <- prune_taxa(!rownames(tax_table(ps)) %in% otus_to_remove, ps)
  return(pruned_ps)
}

# Get the list of OTUs to remove from list_ASV_XXX_removed

otus_to_remove <- unique(c(list_OTU_Mammalia_removed$OTU, 
                           list_OTU_Insecta_removed$OTU,
                           list_OTU_Arachnida_removed$OTU,
                           list_OTU_Aves_removed$OTU))


OTU_to_discard=bind_rows(list_OTU_Arachnida_removed, list_OTU_Aves_removed, list_OTU_Insecta_removed, list_OTU_Mammalia_removed)


OTU_to_discard$marker<-gsub("Marker_","",OTU_to_discard$marker)
# Check if the sam eOTU is in the different marker dataset

otu.discard.summ.marker<-OTU_to_discard %>% group_by(marker)%>% summarise(OTUs=list(OTU))
# check if any 12S or 18S OTU is in COI ps object
intersect(unlist(otu.discard.summ.marker$OTUs[1:2]), taxa_names(rare_otu_phyloseq_objects.animalia$COI)) # 1 OTU shared between 12S/18S and COI 
intersect(unlist(otu.discard.summ.marker$OTUs[2:3]), taxa_names(rare_otu_phyloseq_objects.animalia$`12S`)) # 0 OTUs shared between 12S and the other markers
intersect(unlist(otu.discard.summ.marker$OTUs[1]), taxa_names(rare_otu_phyloseq_objects.animalia$`18S`)) # 3 OTUs sahred between 18S and 12S
intersect(unlist(otu.discard.summ.marker$OTUs[3]), taxa_names(rare_otu_phyloseq_objects.animalia$`18S`)) # 0 OTUs sahred between 18S and COI


rare_otu_phyloseq_Marine.animalia<-rare_otu_phyloseq_objects.animalia

for (marker in c("COI","18S", "12S")) {
  
  #marker=marker[1]
  
  ###extract only marker part of the tables to discard:
  
  
  OTU_to_discard_marker=filter(OTU_to_discard, marker =={{marker}})
  
  
  
  if (nrow(OTU_to_discard_marker)>0) {
    
    rare_otu_phyloseq_Marine.animalia[[marker]]=remove_otus_from_list(rare_otu_phyloseq_Marine.animalia[[marker]], OTU_to_discard_marker$OTU)
    
  }
}

# Print final number of OTUs assogned to marine organisms
dim(tax_table(rare_otu_phyloseq_Marine.animalia$`12S`))[1] # 537
dim(tax_table(rare_otu_phyloseq_Marine.animalia$`18S`))[1] # 371
dim(tax_table(rare_otu_phyloseq_Marine.animalia$`COI`))[1] # 410


#### Check which mammalia and birds are left in the datasets

list.kept.otu<-list()

for (i in seq_along(rare_otu_phyloseq_Marine.animalia)){
  
  marker <- paste0("Marker_", names(rare_otu_phyloseq_Marine.animalia[i]))  # Replace with actual marker names if available
  
  tax.table<-as.data.frame(tax_table(rare_otu_phyloseq_Marine.animalia[[i]]))
  
  mammalia.list<-unique(tax.table[tax.table$Class == "Mammalia","Order"])
  
  aves.list<-unique(tax.table[tax.table$Class == "Aves","Order"])
  
  insecta.list<-unique(tax.table[tax.table$Class == "Insecta","Order"])
  
  arhacn.list<-unique(tax.table[tax.table$Class == "Arhacnida","Order"])
  
  list<-c(mammalia.list,aves.list,insecta.list)
  
  names(list)<-marker
  
  list.kept.otu<-list(list.kept.otu,list)
  
}

list.kept.otu

cat("Mammalia, INsect, arachnida and Birds left in the Dataset")

print(list.kept.otu)

# CHUNK 4: Creating combined markers Table ----
# 1- fusion of identical taxonomy  2-discarding duplicates PER SITE

#################### OBJECTIVES
# This part has the objective of reconciling the three markers data in order to obtain a non-redundant taxonomy table that is our main object to analyze data in a biologically meaningful and comprehensive way. In order to do so we need to:
#   Step 1 - Merge phyloseq objects of the three markers as OTUs tables.
#   Step 2 - Then merge OTUs (defined by their ASVseed) with the same taxonomic assignment (FARTA), PER marker, and keep track of the merging in a specific column.
#   Step 3 - Screen for those OTUs that are not assigned down to the Species level with a particular marker, but may be redundant with other OTUs assigned with a better resolution with another marker.
#   Step 4 - When they occur in the same sample, we hypothesize they are the same and keep only the higher-resolution data, setting the count to 0 for the sample in the sample(s) where they co-occur, and only there. The implicit hypothesis is when they occur alone, they are not the same as the one identified by the other marker in case of co-occurrence; else it would have detected them (i.e., we discard the possibility of false negatives in this particular case).
#   Step 5 - Merge all identical FARTA, whatever the marker, and keep track of the merging. Remove all rows with sum = 0 as a result of previous operations.
#   Step 6 - Create a new phyloseq object that will be used for downstream analysis.

#################### TABLES PRODUCED AND SAVED:
# Step 1: `Combined_OTU_three_markers.tsv` - A table combining taxonomy and presence-absence (psce-abs) data of all three markers, adding a column with the marker name for traceback, and renaming OTUs with the marker name to avoid duplicated names.
# Step 2: `PER_Markers_merged_Tax_OTU_Table.tsv` - The table is transformed by adding for each OTU the FAR (Final Assignment Rank) and FARTA (Final Assignment Rank Taxonomy).
# Step 3: `PER_Markers_merged_Tax_OTU_Table.tsv` - The final table is obtained by grouping combined data by FARTA PER marker, summing counts, and setting all values ≥1 to 1 for presence-absence. Only the first OTU with a given FARTA is kept, and merged OTUs are tracked.
# Step 4: Step 4: all_matches (tibble) and all_matches_df (transformed as datatable): table listing all OTU that do not reach species as FAR (target_asvs), and their 'higher resolution counterpart", that have the same FAR+ additional taxonomic rank => saved as 
# Matches_merged_table: modified at step 4 build a table reconciled for any taxonomic assignment, by screening for each OTU low resolution if they co-occur with their higher resolution counterpart (listed in all_matches) screening each sample one by one. When they do co-ccur in a sample X, set count of the target low resolution otu to 0 for the sample X and go to the next sample => Tables saved as TAX_merged_traceability_log.csv and TAX_merged_table.csv"
# Step 5: final_pruned_table (TAX_merged_traceability_pruned_log.csv" & TAX_merged_final_pruned_table.csv"), we merge all identical FARTA keeping the first OTU/ASVseed we keep trace of the merged identical taxa in the following 12S, 18S and COI columns, we save traceability and OTU table as TAX_merged_traceability_log.csv"and TAX_merged_final_table.csv", row.names = FALSE)
# Step 6: `combined_taxa_phyloseq.rds` - Final phyloseq object reconstructed from the unique table.
# Step 7: compute a density distribution of the number of OTU per taxa (using the traceability inserted in the final_pruned_table), and etxract the 5% taxa outliers on the right tail. Comment: to do tis work properly, we cannot start from asv table (though tempting), but we should insert the traceability table from swarm to include also the size of each OTU before the whole merging. In fact, Tunicata and Vertebrata for ex, that I checked all along because of drastic reduction, have violent redundancy within and among markers

################### A-FUSION of Phyloseq objects and Identical Taxonomy PER MARKER

### Step 1 - Here, we create a presence absence Table with all markers merged

# Combine all marker data with taxonomy and counts, adding "Marker" and retaining OTU names

combined_data <- bind_rows(
  lapply(names(rare_otu_phyloseq_Marine.animalia), function(marker) {
    phyloseq_obj <- rare_otu_phyloseq_Marine.animalia[[marker]]
    
    # Extract taxonomy and counts tables
    tax_table <- as.data.frame(tax_table(phyloseq_obj))
    otu_table <- as.data.frame(otu_table(phyloseq_obj))
    
    # Add OTU names as a column
    otu_table <- otu_table %>% rownames_to_column(var = "OTU")
    
    # Remove "initial_phylum_assignment" and add Marker column
    tax_table <- tax_table %>%
      select(-initial_phylum_assignment) %>%
      mutate(Marker = marker)
    
    # Combine taxonomy and OTU count data
    combined_table <- bind_cols(tax_table, otu_table)
    
    return(combined_table)
  })
)
# # Write the combined data to a TSV file
# write.table(combined_data, "Preparation_combined_dataset/01_Combined_OTU_three_markers.tsv", sep = "\t", row.names = FALSE, quote = FALSE)



## Step 2- We merge OTU with the same FARTA ----
# (same assignation to the same rank) per marker (we don't merge markers)

#### a- Here we define FAR -final assignation rank- and FARTA -final assignation rank taxonomy- for each ASVseed/OTU and add columns with those information
#### Prepare variable needed for downstream analysis

# Define metadata columns typically present before sample columns
metadata_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Marker")

# Identify sample columns by excluding metadata columns and selecting only numeric columns
sample_cols <- names(combined_data)[!(names(combined_data) %in% metadata_cols) & sapply(combined_data, is.numeric)]

# Define taxonomy levels from species to phylum
taxonomy_levels <- c("Species", "Genus", "Family", "Order", "Class", "Phylum")


# Function to determine FAR and FARTA, skipping over NA values in ascending order
get_far_farta <- function(tax_row) {
  for (rank in taxonomy_levels) {
    if (!is.na(tax_row[[rank]])) {
      # Return the first non-NA level found in ascending order
      return(list(FAR = rank, FARTA = tax_row[[rank]]))
    }
  }
  
  # If all levels are NA, return NA for both FAR and FARTA
  return(list(FAR = NA, FARTA = NA))
}

###### b- Apply the function to each row in combined_data and add FAR/FARTA columns
far_farta_results <- map_dfr(1:nrow(combined_data), function(i) {
  get_far_farta(combined_data[i, ])
})

# Add FAR and FARTA columns to combined_data
combined_data <- combined_data %>%
  mutate(
    FAR = far_farta_results$FAR,
    FARTA = far_farta_results$FARTA
  )

# Ensure that row names are added as the "OTU" column
combined_data <- combined_data %>%
  mutate(OTU = paste0(Marker, "_", OTU))


###### Validate the taxonomical assignation of FARTA with WormS  

# Check that the column exist
if (!"FARTA" %in% colnames(combined_data)) {
  stop("La colonne 'species' est absente du fichier.")
}

# Step 2 : Validate names with WORMS and provide the valid name if needed
validate_species_name <- function(species_name) {
  records <- tryCatch(
    wm_records_name(name = species_name),
    error = function(e) NULL
  )
  
  if (!is.null(records) && nrow(records) > 0) {
    # Vérifier si le nom est accepté
    if (!is.na(records$valid_name[1])) {
      valid_name <- records$valid_name[1]
    } else {
      valid_name <- species_name
    }
    return(data.frame(
      original_name = species_name,
      valid_name = valid_name,
      status = ifelse(species_name == valid_name, "Accepted", "Synonym")
    ))
  } else {
    return(data.frame(
      original_name = species_name,
      valid_name = NA,
      status = "Not found"
    ))
  }
}

# Extract unique FARTA to run the validation
combined_data.unique.farta<-combined_data %>% group_by(FARTA) %>%
  summarize(freq=n())

# Apply the function above to the whole list of species
validated_names <- do.call(rbind, lapply(combined_data.unique.farta$FARTA, validate_species_name))

# Step 3 : Look and add AphiaIDs for the lines with a valid_name (non NA)
get_aphia_id <- function(species_name) {
  if (is.na(species_name)) {
    return(NA)  # Retourner NA si le nom valide est absent
  }
  
  records <- tryCatch(
    wm_records_name(name = species_name),
    error = function(e) NULL
  )
  
  if (!is.null(records) && nrow(records) > 0) {
    return(records$AphiaID[1])  # Retourner le premier AphiaID trouvé
  } else {
    return(NA)  # Retourner NA si aucun résultat
  }
}

# Apply the AphiaID to the whole list of species with a valid name
cat("Get Aphia ID for valid name...\n")
validated_names$AphiaID <- sapply(validated_names$valid_name, get_aphia_id)

# saveRDS(validated_names,file="Preparation_combined_dataset/tmp_validated.RData")
# 
# write.csv(validated_names,file="Preparation_combined_dataset/02_Unique_FARTA_validated_names_summerTEST.csv")
# 
# ########### IMPORTANT: 
############## THIS FIRST CHUNK is followed by a MANUAL check (excell) of the output file to understand why some species were not found. Some were actually terrestrial or freshwater species. Other were not spelled correctly.
############## After the manual checking, the new file is saved: 02_Unique_FARTA_validated_names_manual_check.cs


# Reload the new file with manual check done

validated_names_manual<-read.csv2( "Preparation_combined_dataset/02_Unique_FARTA_validated_names_manual_check.csv") %>%
  select(-Column1)

colnames(validated_names_manual)[1]<-c("FARTA") # Change the "original_name" column into "FARTA"

# Add AphiaID to those taxa that have been manually corrected

no_aphia<-subset(validated_names_manual, is.na(AphiaID))

no_aphia$AphiaID<-sapply(no_aphia$valid_name, get_aphia_id)

validated_names_manual<-subset(validated_names_manual, !valid_name  %in% no_aphia$valid_name)

validated_names_manual<-rbind(validated_names_manual, no_aphia)


# Add the valid name to the file
combined_data <-  merge(validated_names_manual,combined_data, by.x = "FARTA", all.x = TRUE)

# Step 4 : Rename the "valid_name" column as FARTA and the "FARTA" column as "old_unvalidated_FARTA"

combined_data<-combined_data %>% mutate(old_unchecked_FARTA=FARTA)%>%
  select(-FARTA) %>%
  rename(FARTA=valid_name)

#Step 5 Re-organize the column to get first the Aphia ID then the valid species name, then the status, then the original name

combined_data <- combined_data %>%
  relocate(AphiaID, FARTA, status, old_unchecked_FARTA, .before = everything())  # Reorganize columns
head (combined_data)

# Report the rownames as OTU_names

rownames(combined_data)<- combined_data$OTU

#Save the combined data with the valid name

# write.csv(combined_data,file="Preparation_combined_dataset/03_Combined_OTU_three_markers_Valid_Names.csv")

# Group by FARTA and Marker and summarize
final_table <- combined_data %>%
  group_by(FARTA, Marker) %>%
  summarise(
    # Retain the first OTU as the seed
    OTU_seed_retained = first(OTU),
    
    # Collect merged OTUs (excluding the retained OTU seed)
    merged_OTUs = list(OTU[-1]),
    
    # Retain taxonomy information
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    # Retain the Final Assignment Rank (FAR)
    FAR = first(FAR),
    
    # Retain the Final AphiaID 
    AphiaID = first(AphiaID),
    
    # Reconvert sample counts to presence/absence
    # for semiquantitative data change shut off this part
    across(all_of(sample_cols), ~ ifelse(sum(.) >= 1, 1, 0)),
    
    .groups = "drop"
  )

# Convert merged_OTUs lists to comma-separated strings
final_table <- final_table %>%
  mutate(
    merged_OTUs = sapply(merged_OTUs, function(x) if (length(x) == 0) "" else paste(x, collapse = ", "))
  )

# Replace any NA values with empty strings
final_table[is.na(final_table)] <- ""


###Keep track of this Table
# write.table(final_table, "Preparation_combined_dataset/04_PER_Markers_merged_Tax_OTU_Table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

## Now, build a  table reconciled for any taxonomic assignment ----
# by screening for each OTU low resolution if they co-occur with their counterpart in each sample. When they do, set to 0 for the target otu the sample where it co-occur with one of its higher resolution counterpart

### a Prepare tax levels and traceability table for downstream analysis
# Define taxonomy levels in hierarchical order
taxonomy_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Create the traceability table to log changes
traceability_taxonomy_reconciliation <- data.frame(
  erased_OTU = character(),
  erased_marker = character(),
  kept_OTU = character(),
  kept_marker = character(),
  sample = character(),
  stringsAsFactors = FALSE
)

# b: function to find higher resolution matches in taxonomic levels FAR

find_higher_resolution_matches <- function(final_table, target_otu, taxonomy_levels) {
  # Extract target OTU row
  current_row <- final_table %>% filter(OTU_seed_retained == target_otu)
  
  # Define the current FAR level and next level
  far_index <- match(current_row$FAR, taxonomy_levels)
  current_levels <- taxonomy_levels[1:far_index]
  next_level <- taxonomy_levels[far_index + 1]
  
  # Get the marker of the target OTU
  target_marker <- current_row$Marker
  
  # Filter for higher resolution matches from different markers
  higher_resolution_matches <- final_table %>%
    filter(
      # Exclude the target OTU itself
      OTU_seed_retained != target_otu,
      
      # Match taxonomy up to the FAR level
      Reduce(`&`, lapply(current_levels, function(level) {
        final_table[[level]] == current_row[[level]]
      })),
      
      # Ensure the next level is populated
      !is.na(.data[[next_level]]) & .data[[next_level]] != "",
      
      # Exclude matches from the same marker as the target OTU
      Marker != target_marker
    )
  
  return(higher_resolution_matches)
}


# c: Filter OTUs with FAR level not equal to "Species" that will be target OTU for this 'reconcile and avoid redundancy' section
target_otus <- final_table %>%
  filter(FAR != "Species") %>%
  pull(OTU_seed_retained) %>%
  unique()

# Initialize a data structure to store matches
all_matches <- list()

# d: Loop over each target OTU and apply the function to find higher-resolution matches
for (target_otu in target_otus) {
  #cat("Processing:", target_otu, "\n")
  
  # Find higher resolution matches for the current target OTU
  matches <- find_higher_resolution_matches(final_table, target_otu, taxonomy_levels)
  
  # Store results in the list if matches were found
  if (nrow(matches) > 0) {
    all_matches[[target_otu]] <- matches
  }
}

# e: Combine all results into a single data frame if needed
all_matches_df <- bind_rows(all_matches, .id = "Target_OTU")

# View the resulting matches
print(all_matches_df)

##f-Now, we screen the listed OTU with low resolution and all their counterparts with higher level of assignment, and we assign their number to 0 WHEN THEY CO-OCCUR with an OTU with the same initial taxonomic assignment as the target OTU FARTA, but higher taxonomic resolution (FARTA+ or ++ or +++...), and we create a traceability table to keep track of all changes per OTU target, OTU kept and sample for which OTU target count has been set to 0

# Initialize traceability log
traceability_log <- data.frame(
  erased_OTU = character(),
  erased_marker = character(),
  kept_OTU = character(),
  kept_marker = character(),
  sample = character(),
  stringsAsFactors = FALSE
)

# Process each target OTU
unique_target_otus <- unique(all_matches_df$Target_OTU)

# Proceed each target OTU along the list of matches, then along the list of samples
for (target_otu in unique_target_otus) {
  # Extract target OTU row
  target_row <- final_table %>% filter(OTU_seed_retained == target_otu)
  
  # Extract FARTA for target OTU
  target_farta <- target_row$FARTA
  
  # Get all higher-resolution matches for this target OTU
  higher_matches <- all_matches_df %>% filter(Target_OTU == target_otu)
  
  # Loop through each sample
  for (sample in sample_cols) {
    for (i in seq_len(nrow(higher_matches))) {
      # Extract higher-resolution OTU info
      higher_otu <- higher_matches$OTU_seed_retained[i]
      higher_marker <- higher_matches$Marker[i]
      
      # Extract FARTA for the higher-resolution OTU
      higher_farta <- final_table %>%
        filter(OTU_seed_retained == higher_otu) %>%
        pull(FARTA)
      
      # Check if the higher-resolution OTU is present in the same sample
      higher_otu_present <- final_table %>%
        filter(OTU_seed_retained == higher_otu) %>%
        pull(!!sym(sample)) >= 1 #garde le nom de l'échantillon si l'OTU "higher' a un count =>1
      
      # Check if a same-marker OTU with the same FARTA exists in the same sample
      same_marker_farta_exists <- final_table %>%
        filter(
          Marker == target_row$Marker,          # Same marker as target OTU
          FARTA == higher_farta                # Same FARTA as target OTU
        ) %>%
        pull(!!sym(sample))
      
      # Simplify the condition: Replace only if higher OTU exists and no same-marker FARTA exists
      if (target_row[[sample]] >= 1 && higher_otu_present && 
          (is.null(same_marker_farta_exists) || !any(same_marker_farta_exists))) {
        # Set target OTU in this sample to 0
        final_table[final_table$OTU_seed_retained == target_otu, sample] <- 0
        
        # Log the change
        traceability_log <- rbind(traceability_log, data.frame(
          erased_OTU = target_otu,
          erased_marker = target_row$Marker,
          kept_OTU = higher_otu,
          kept_marker = higher_marker,
          sample = sample,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}


# Output results
cat("Number of OTUs modified:", length(unique(traceability_log$erased_OTU)), "\n")
#cat("Traceability Log:\n")
#print(traceability_log)
# write.table(traceability_log, "Preparation_combined_dataset/TAX_merged_traceability_log.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
# # 
# write.table(final_table, "Preparation_combined_dataset/TAX_merged_final_table.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

####### Step 5- merge all identical FARTA, whatever the marker, and keep track of the merging, remove all rows with sum= 0 as a result of previous operations
# Group by FARTA and merge rows with identical FARTA
final_pruned_table <- final_table %>%
  
  group_by(FARTA) %>%
  summarize(
    
    collapsed_OTU_seeds=paste(OTU_seed_retained, collapse = " ,"),
    
    collapsed_Markers=paste(Marker, collapse = " ,"),
    
    # Retain the first OTU seed for the group
    
    OTU_seed_retained = first(OTU_seed_retained),
    
    # Retain the first Marker for the group (assuming Marker is consistent per FARTA)
    Marker = first(Marker),
    
    # Combine merged OTUs and previously merged OTUs
    merged_OTUs = paste(
      unique(c(
        unlist(strsplit(paste(merged_OTUs[merged_OTUs != ""], collapse = ", "), ", ")), # Previously merged OTUs
        OTU_seed_retained[-1] # Other OTUs in the group
      )),
      collapse = ", "
    ),
    
    # Retain taxonomy information
    Kingdom = first(Kingdom),
    Phylum = first(Phylum),
    Class = first(Class),
    Order = first(Order),
    Family = first(Family),
    Genus = first(Genus),
    Species = first(Species),
    
    # Retain the Final Assignment Rank (FAR)
    FAR = first(FAR),
    
    # Retain the AphiaID
    AphiaID = first(as.character(AphiaID)),
    
    # Merge sample counts by summing and reconverting to presence/absence
    # for semiquantitative data change shut off this part or introduce the semiquantitative transformation here
    
    across(where(is.numeric), ~ ifelse(sum(.) >= 1, 1, 0)),
    
    .groups = "drop"
  )


# Replace any NA values with empty strings
final_pruned_table[is.na(final_pruned_table)] <- ""

# Save the resulting pruned table

# write.table(final_pruned_table, "Preparation_combined_dataset/05_Final_Pruned_Table_With_Traceability.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#Finally, recreate a phyloseq object for downstream analysis ----

colnames(final_pruned_table)
# Identify taxonomy and sample columns
taxonomy_cols <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "FARTA", "FAR", "AphiaID", "collapsed_Markers")

# Create taxonomy_table with OTU_seed_retained as rownames
rownames(final_pruned_table)<-NULL
taxonomy_table <- final_pruned_table %>%
  select(OTU_seed_retained, all_of(taxonomy_cols)) %>%
  column_to_rownames("OTU_seed_retained") %>%
  as.matrix()

taxa_names<-rownames(taxonomy_table)

# Extract count table with OTU_seed_retained as rownames
count_table <- final_pruned_table %>%
  select(OTU_seed_retained, all_of(sample_cols)) %>%
  column_to_rownames("OTU_seed_retained") %>%
  as.matrix()


sam<- sample_data(rare_otu_phyloseq_Marine.animalia[[1]]) %>% data.frame()

combined_taxa_phyloseq <- phyloseq(otu_table(count_table, taxa_are_rows = TRUE),
                                   tax_table(taxonomy_table), taxa_names(taxa_names),
                                   sample_data(sam))

### Introduce accepted full taxonomy through Worms ----

full.tax<-as.data.frame(tax_table(combined_taxa_phyloseq))%>% rownames_to_column(var="OTU_seed")

classification.data<-data.frame("Kingdom", "Subkingdom","Infrakingdom","Phylum", "Phylum (Division)","Subphylum","Subphylum (Subdivision)","Infraphylum",
                                "Parvphylum","Gigaclass", "Megaclass","Superclass","Class", "Subclass",      
                                "Order", "AphiaID")

colnames(classification.data)<-c("Kingdom", "Subkingdom","Infrakingdom","Phylum","Phylum (Division)","Subphylum","Subphylum (Subdivision)","Infraphylum",
                                 "Parvphylum","Gigaclass", "Megaclass","Superclass","Class", "Subclass",     
                                 "Order","AphiaID")

classification.data[1,]<-NA

for (i in as.numeric(full.tax[,"AphiaID"])) {
  
  fis<-wm_classification(id=i , marine_only=FALSE)
  fis<-fis[,c("rank","scientificname")]
  fis<-as.data.frame(t(fis))
  colnames(fis)<-fis[1,]
  fis<-fis[-1,]
  fis$AphiaID<-i
  classification.data<-plyr::rbind.fill(classification.data,fis)
  classification.data<-classification.data[,c("Kingdom","Phylum",
                                              "Class","Order",
                                              "Family","Genus",
                                              "Species", "AphiaID")]
}

# Substitute in the Phylum column the Chordata into Vertebrata and Tuicata. Use the Phylum column of the orginal taxonomic table as guide (which was modified at the step 2 of this script.

classification.data<-classification.data[-1,]
colnames(classification.data)<-c("Kingdom_Worms_Accepted","Phylum_Worms_Accepted",
                                 "Class_Worms_Accepted","Order_Worms_Accepted",
                                 "Family_Worms_Accepted","Genus_Worms_Accepted",
                                 "Species_Worms_Accepted","AphiaID")

tax_Worm_accepted<-left_join(classification.data,full.tax,by="AphiaID")
tax_Worm_accepted<-column_to_rownames(tax_Worm_accepted,var="OTU_seed")

#apply the substitution for Chordata into Vertebrata and Tunicata using the old "Phylum" column

tax_Worm_accepted$Phylum_Worms_Accepted_old<-tax_Worm_accepted$Phylum_Worms_Accepted
tax_Worm_accepted <- tax_Worm_accepted %>% mutate(Phylum_Worms_Accepted=case_when(Phylum %in% "Vertebrata"~"Vertebrata",
                                                                                  Phylum %in% "Tunicata"~"Tunicata",
                                                                                  Phylum %in% "Cephalochordata"~"Cephalochordata",
                                                                                  Phylum %in% "Unclassified-Chordata"~"Unclassified-Chordata",
                                                                                  TRUE ~ Phylum_Worms_Accepted))


# By checking the old and new classification I noticed that for 3 taxa the old and new classification was not correct
# These were among the taxa that had a synonym in Worm and their name was substituted by th e "validate_species_name" function. 
# The validate_species_name assigned the worng aphiaID so the full taxonomical classificaiton is not correct and need to be fixed

worms_status<-validated_names_manual %>% select (FARTA,status) 
worms_status$status<-as.character(worms_status$status)
tax_Worm_accepted$OTU_seed<-rownames(tax_Worm_accepted)

tax_Worm_accepted_status<-left_join(tax_Worm_accepted,worms_status,by="FARTA")
# write.csv(tax_Worm_accepted_status,file="Preparation_combined_dataset/Taxonomy_classification_mistake.csv")

# Only 3 taxa had an incorrect lineage 
taxa_to_correct<-read.csv("Preparation_combined_dataset/Taxonomy_classification_mistake.csv")


correction.data<-data.frame("Kingdom", "Subkingdom","Infrakingdom","Phylum", "Phylum (Division)","Subphylum","Subphylum (Subdivision)","Infraphylum",
                            "Parvphylum","Gigaclass", "Megaclass","Superclass","Class", "Subclass",      
                            "Order", "AphiaID")

colnames(correction.data)<-c("Kingdom", "Subkingdom","Infrakingdom","Phylum","Phylum (Division)","Subphylum","Subphylum (Subdivision)","Infraphylum",
                             "Parvphylum","Gigaclass", "Megaclass","Superclass","Class", "Subclass",     
                             "Order","AphiaID")

correction.data[1,]<-NA

correct.AphiaID<-taxa_to_correct$Correct.Aphia.ID

for (i in as.numeric(correct.AphiaID)){
  
  fis<-wm_classification(id= i, marine_only=FALSE)
  fis<-fis[,c("rank","scientificname")]
  fis<-as.data.frame(t(fis))
  colnames(fis)<-fis[1,]
  fis<-fis[-1,]
  fis$AphiaID<-i
  correction.data<-plyr::rbind.fill(correction.data,fis)
}

colnames(correction.data)
correction.data<-correction.data %>% select (Kingdom,Phylum,
                                             Class,Order,
                                             Family,Genus,
                                             AphiaID)
correction.data<-correction.data[-1,]

# Change the classification of the taxa from the tax_Worm_accepted 
# tax_Worm_accepted<-as.data.frame(combined_taxa_phyloseq_WormsAcc@tax_table)

tax_Worm_correction<-tax_Worm_accepted  %>% subset(OTU_seed %in% taxa_to_correct$OTU_seed)%>% 
  
  mutate(Kingdom_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Kingdom[1],
                                            OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Kingdom[2],
                                            OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Kingdom[3])) %>%
  
  mutate(Phylum_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Phylum[1],
                                           OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Phylum[2],
                                           OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Phylum[3])) %>%
  
  mutate(Class_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Class[1],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Class[2],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Class[3])) %>% 
  
  mutate(Order_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Order[1],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Order[2],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Order[3]))  %>% 
  
  mutate(Family_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Family[1],
                                           OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Family[2],
                                           OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Family[3])) %>% 
  
  mutate(Genus_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Genus[1],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Genus[2],
                                          OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Genus[3]))%>% 
  
  mutate(Species_Worms_Accepted = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Species[1],
                                            OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Species[2],
                                            OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Species[3]))%>% 
  
  mutate(AphiaID = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Correct.Aphia.ID[1],
                             OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Correct.Aphia.ID[2],
                             OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Correct.Aphia.ID[3]))%>% 
  
  mutate(FARTA = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Genus[1],
                           OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Genus[2],
                           OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Genus[3])) %>%
  
  mutate(Phylum_Worms_Accepted_old = case_when(OTU_seed %in% taxa_to_correct$OTU_seed[1]~taxa_to_correct$Phylum[1],
                                               OTU_seed %in% taxa_to_correct$OTU_seed[2]~taxa_to_correct$Phylum[2],
                                               OTU_seed %in% taxa_to_correct$OTU_seed[3]~taxa_to_correct$Phylum[3])) 




tax_Worm_accepted<-tax_Worm_accepted %>% 
  subset(!OTU_seed %in% taxa_to_correct$OTU_seed)

tax_Worm_accepted<-rbind(tax_Worm_accepted, tax_Worm_correction)
rownames(tax_Worm_accepted)<-tax_Worm_accepted$OTU_seed


combined_taxa_phyloseq_WormsAcc<-phyloseq(otu_table(otu_table(combined_taxa_phyloseq,taxa_are_rows = TRUE)),
                                          tax_table(as.matrix(tax_Worm_accepted)), taxa_names(rownames(tax_Worm_accepted)),
                                          sample_data(sample_data(combined_taxa_phyloseq)))


#  CHUNK 5: Select only summer samples

combined_taxa_ps_summer<-subset_samples(combined_taxa_phyloseq_WormsAcc, Season %in% "Summer")
# #
combined_taxa_ps_summer<-prune_taxa(taxa_sums(combined_taxa_ps_summer)!=0,combined_taxa_ps_summer)
# 
combined_taxa_ps_summer<-combined_taxa_phyloseq_WormsAcc

# # Select only the OTU-seeds classified at Species level

combined_taxa_ps_summer.sp<-prune_taxa(as.data.frame(tax_table(combined_taxa_ps_summer))[,"FAR"]=="Species",combined_taxa_ps_summer)

# Final Taxonomiy table at species level 
ps_538_species_summer2021<-combined_taxa_ps_summer.sp


