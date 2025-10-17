# Diversity analyses (Alpha and Beta) as reported in the article "Expanding on 
# the portuarization syndrome from an ecological perspective: eDNA reveals rich diversity, 
# non-indigenous hotspots, and biotic homogenization in ports"
# The phyloseq object called in this script has been obtained by the combination 
# of the OTU and TAXA generated for three different markers 12S,18S, 
# COI accoridng to the "Preparation table script"
# Additional filtering for organisms classified only at species level and for samples 
# collected only during summer 2021 has been done.


# The analyses performed are organized as follow:

# 1) General information about dataset and communities
# 2) section a. - Alpha diversity analyses on whole community (COMM) and resident species only (RESS) 
# 3) section b., c., e., f. - Beta diversity analyses on whole community (COMM) and resident species only (RESS)
# including analyses n beta diversity components (turnover and nestedness), biotic homogenization and partial dbRDA.
# 4) section d. Venn diagrams for unique and shared species in ports and in natural habitats
# 5) section g-p: diversity analyses on NIS community



# Author: Ginevra Lilli 

# Load libraries -----------------------------------------------------------------
pkgs <- c("phyloseq","ggplot2","dplyr","vegan","tidyr",
          "cowplot","gridExtra","broom","worrms","car","purrr",
        "tibble","ggpubr","MicEco","readxl","indicspecies",
        "adespatial","stringr","pheatmap","grid","ggforce","betapart")


lapply(pkgs, require,character.only=TRUE)

# color palettes ---------------------------------------------------------------
c25<-get(load("C:/Users/Ginevra Lilli/Dropbox/Postdoc_CEAB/MarEEE-Med/Analyses/colors_barplot25.Rdata")) #load 25 colors palette for barplot 
new.pal<-c(
  "#E69F00",  # orange
  "#56B4E9",  # sky blue
  "#009E73",  # bluish green
  "#F0E442",  # yellow
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#CC79A7",  # reddish purple
  "#999999",  # grey
  "#117733",  # dark green
  "#332288",  # dark blue
  "#88CCEE",  # light blue
  "#44AA99",  # turquoise
  "#DDCC77",  # sand
  "#AA4499",  # purple-pink
  "#882255",  # wine red
  "#661100",  # brown
  "#6699CC",  # soft blue
  "#888888",  # mid grey
  "#F1A340",  # tan-orange
  "#998EC3")   # muted purple)

# general estetics of plots
fig <-theme(strip.background = element_blank(), strip.text.y = element_text(size=16)) +
  theme(legend.background = element_blank()) +
  theme(legend.key = element_blank()) +
  theme(legend.text = element_text(size=14)) +
  theme(legend.title  = element_text(size=15)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme (axis.text.x = element_text(size=13,  vjust=0))+ #,angle=90,
  theme (axis.text.y = element_text(size=13))+
  theme (axis.title.x  = element_text(size=15))+
  theme (axis.title.y= element_text(size=15))+
  theme(title = element_text(size=15))


# START -------------------------------------------------------------------------
# Load R data object with 538 species obatined from collapsing 3 markers, filtering species at species level and only samples from summer 2021

get(load("ps_538_species_summer2021.RData"))

# Load the species occurrence table 
traits.ninis.df<-read.csv2("C:/Users/Ginevra Lilli/Dropbox/Postdoc_CEAB/MarEEE-Med/Writing/Article/Manuscrip_MarEEEMed/MS_Reviewed/Final_to_coauthors/Species_occurrence_table_Watersipora_check.csv")

# Generate full community (COMM=NIS + LOCAL), RESS species (RESS) and non-indigenous species (NIS) objects

COMM <- combined_taxa_ps_summer.sp

metaCOMM<-sample_data(COMM)%>% data.frame()

RESS<-subset_taxa(combined_taxa_ps_summer.sp, !NINIS_Native %in% c("NIS"))

RESS<-prune_samples(sample_sums(RESS)!=0,RESS) 

df<-sample_data(RESS)%>% data.frame()

NIS<-subset_taxa(combined_taxa_ps_summer.sp, NINIS_Native %in% c("NIS"))

NIS<-prune_samples(sample_sums(NIS)!=0,NIS) 

metaNIS<-sample_data(NIS)%>% data.frame()

# 0. General info about dataset and community ----

# Determine how many species were classified with 1, 2 or 3 markers and which taxa 

unique_phylum<-unique(as.data.frame(tax_table(COMM))$Phylum_Worms_Accepted)
phylum_colors <- setNames(c25[1:length(unique_phylum)], unique_phylum)

freq.collaps.markers<-traits.ninis.df %>% group_by(collapsed_Markers, Phylum_Worms_Accepted)%>%summarise(freq=n(),
                                                                                                         freq.perc=(n()*100)/537)

total_species_by_marker<- freq.collaps.markers %>% group_by(collapsed_Markers) %>% summarise(tot.spec.by.mark=sum(freq))
freq.collaps.markers <- left_join(freq.collaps.markers,total_species_by_marker, by="collapsed_Markers")
  
freq.collaps.markers$perc.by.taxa<- (freq.collaps.markers$freq*100)/freq.collaps.markers$tot.spec.by.mark

# Set graphic env
theme_set(theme_minimal())


# Generate an upset plot for the combined markers 

library(UpSetR)

 input.upset<-c(
#   `12S`=160,
#   `18S`=450,
#   COI=313,
  "12S&18S"=1,
  "12S&COI"=53,
  "18S&COI"=43,
  "12S&18S&COI"=1,
  "18S"=157,
  "12S"=90,
  "COI"=193
  
)

upset.markers.plot<-upset(fromExpression(input.upset), 
      nintersects = 7, 
      nsets = 3, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 2.1, 
      point.size = 2.8, 
      line.size = 1
)



# Define the desired order
desired_order <- c(  "18S","COI",  "12S", "12S ,COI","18S ,COI","12S ,18S", "12S ,18S ,COI")

# Apply that order to the factor
freq.collaps.markers$collapsed_Markers <- factor(freq.collaps.markers$collapsed_Markers, 
                                                 levels = desired_order)

markers.plot.by.taxa<-ggplot(freq.collaps.markers,aes(x=collapsed_Markers, y=perc.by.taxa, fill=Phylum_Worms_Accepted))+ 
  geom_bar(stat="identity", position = "stack")+
  ggtitle("Number of Species classified by combining the 3 markers","Classificatio at Phylum level")+
  scale_fill_manual(values = phylum_colors)+fig+
  theme(legend.position = "bottom")

# Average number of species per locality----- 

sample.per.loc<-metaCOMM %>% group_by(Locality)%>%summarise(samples=sample_name)
otu<-as.data.frame(otu_table(COMM))

otu.pooled<-as.data.frame(matrix(nrow=538,ncol=12))
colnames(otu.pooled)[1:12]<-unique(metaNIS$Locality)
rownames(otu.pooled)<-taxa_names(COMM)

loc<-unique(sample.per.loc$Locality)

for (i in loc){
  sub<-subset (sample.per.loc,Locality %in% i)
  try<- otu %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled[,i]<-try[,i]
}

otu.pooled[otu.pooled>1]<-1 # transform in presence absence 
mean(colSums(otu.pooled))# 125.9
sd(colSums(otu.pooled))#23.1



# Barplot taxa by localities pooled ----

theme_set(theme_minimal())
barplot_habitat.loc<-function(ps, taxrank){
  ps.ra<-transform_sample_counts(ps, function(x){x/sum(x)})
 
   plot.bar<-plot_bar(ps.ra, x= "Locality", fill={{taxrank}})+
    
     geom_bar(aes(color=NULL), position="stack", stat="identity")+
     
    facet_grid(~Type_habitat_broad,scale="free_x",
               labeller = labeller(Type_habitat_broad = c(
                 "Artificial" = "Port",
                 "Natural" = "Natural")))+
    scale_fill_manual(values = phylum_colors)+
     scale_y_continuous(limits = c(0,1.08))+
    scale_x_discrete(limits= c("Port-la-Nouvelle", "Sete_1", "Sete_2", "Seyne-sur-mer",
                               "Toulon",
                               "Saint-Tropez", "Nice", "Calvi", "Ajaccio", "Bonifaccio", 
                               "Porto-Vecchio", "Bastia"),
                     labels=c("Port-la-Nouvelle", "Sete_1", "Sete_2", "Seyne-sur-mer",
                              "Toulon",
                              "Saint-Tropez", "Nice", "Calvi", "Ajaccio", "Bonifacio", 
                              "Porto-Vecchio", "Bastia"))+
    guides(fill = guide_legend(title = "Phylum"))+
    fig+
    theme(axis.text.x = element_text(size=13, angle=90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.text = element_text(size=16))
  return(plot.bar)
}


#RESS

df<-sample_data(RESS)%>%data.frame()

sample.per.loc.RESS<-df %>% group_by(Locality, Type_habitat_broad)%>%summarise(samples=sample_name)
sample.per.loc.RESS$combined<-paste0(sample.per.loc.RESS$Locality,"_",sample.per.loc.RESS$Type_habitat_broad)
otu.RESS<-as.data.frame(otu_table(RESS))


otu.pooled.RESS<-as.data.frame(matrix(nrow=495,ncol=24))
colnames(otu.pooled.RESS)[1:12]<-paste0(unique(df$Locality),"_Natural")
colnames(otu.pooled.RESS)[13:24]<-paste0(unique(df$Locality),"_Artificial")
rownames(otu.pooled.RESS)<-taxa_names(RESS)

loc.hab.ress<-unique(sample.per.loc.RESS$combined)

for (i in loc.hab.ress){
  sub<-subset (sample.per.loc.RESS,combined %in% i)
  try<- otu.RESS %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled.RESS[,i]<-try[,i]
}

df$combined.loc.hab<-paste0(df$Locality,"_",df$Type_habitat_broad)

df.sub<-subset(df, !duplicated(combined.loc.hab ))
rownames(df.sub)<-df.sub$combined.loc.hab

# How many species by each locality and habitat?
species.count<-otu.pooled.RESS
species.count[species.count>1]<-1
species.count<-colSums(species.count)

# remake ps.subset
RESS.sub<-phyloseq(otu_table(as.matrix(otu.pooled.RESS),taxa_are_rows = TRUE),
                    tax_table(tax_table(RESS)),taxa_names(taxa_names(RESS)),
                    sample_data(df.sub))


barplot.RESS.pooled<-barplot_habitat.loc(RESS.sub, "Phylum")



# NIS 

sample.per.loc<-metaNIS %>% group_by(Locality, Type_habitat_broad)%>%summarise(samples=sample_name)
sample.per.loc$combined<-paste0(sample.per.loc$Locality,"_",sample.per.loc$Type_habitat_broad)
otu<-as.data.frame(otu_table(NIS))


otu.pooled<-as.data.frame(matrix(nrow=43,ncol=24))
colnames(otu.pooled)[1:12]<-paste0(unique(metaNIS$Locality),"_Natural")
colnames(otu.pooled)[13:24]<-paste0(unique(metaNIS$Locality),"_Artificial")
rownames(otu.pooled)<-taxa_names(NIS)

loc.hab<-unique(sample.per.loc$combined)

for (i in loc.hab){
  sub<-subset (sample.per.loc,combined %in% i)
  try<- otu %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled[,i]<-try[,i]
}

metaNIS$combined.loc.hab<-paste0(metaNIS$Locality,"_",metaNIS$Type_habitat_broad)

metaNIS.sub<-subset(metaNIS, !duplicated(combined.loc.hab ))
rownames(metaNIS.sub)<-metaNIS.sub$combined.loc.hab

# How many species by each locality and habitat?
species.count<-otu.pooled
species.count[species.count>1]<-1
species.count<-colSums(species.count)

# remake ps.subset
NIS.sub<-phyloseq(otu_table(as.matrix(otu.pooled),taxa_are_rows = TRUE),
                    tax_table(tax_table(NIS)),taxa_names(taxa_names(NIS)),
                    sample_data(metaNIS.sub))

barplot.NIS.pooled<-barplot_habitat.loc(NIS.sub, "Phylum")


barplot_pooled<-ggarrange(barplot.RESS.pooled,
                          barplot.NIS.pooled, nrow=1, ncol=2, common.legend = TRUE)


## Add species numbers on barplot

# Numeber of unique species per each locality and habitat: 
unique_localities<-unique(sample_data(RESS)$Locality)

# NoNIS: Number of unique species in ports by locality
for (locality in unique_localities){
  loc.hab<-subset_samples(RESS, Locality%in%locality)
  loc.hab<-subset_samples(loc.hab, Type_habitat_broad%in%"Artificial")
  loc.hab<-prune_taxa(taxa_sums(loc.hab)!=0,loc.hab)
  print(locality)
  print("Port")
  print(ntaxa(loc.hab))
}

# NoNIS: Numera of unique species in ports by locality
for (locality in unique_localities){
  loc.hab<-subset_samples(RESS, Locality%in%locality)
  loc.hab<-subset_samples(loc.hab, Type_habitat_broad%in%"Natural")
  loc.hab<-prune_taxa(taxa_sums(loc.hab)!=0,loc.hab)
  print(locality)
  print("Natural")
  print(ntaxa(loc.hab))
}

# NIS: Numera of unique species in ports by locality
for (locality in unique_localities[-5]){
  loc.hab<-subset_samples(NIS, Locality%in%locality)
  loc.hab<-subset_samples(loc.hab, Type_habitat_broad%in%"Natural")
  loc.hab<-prune_taxa(taxa_sums(loc.hab)!=0,loc.hab)
  print(locality)
  print("Natural")
  print(ntaxa(loc.hab))
}
# NIS: Numera of unique species in ports by locality

for (locality in unique_localities[]){
  loc.hab<-subset_samples(NIS, Locality%in%locality)
  loc.hab<-subset_samples(loc.hab, Type_habitat_broad%in%"Artificial")
  loc.hab<-prune_taxa(taxa_sums(loc.hab)!=0,loc.hab)
  print(locality)
  print("Ports")
  print(ntaxa(loc.hab))
}


## a. Lower richness of species in ports than outside ----
## Alpha diversity analyses performed on whole community (COMM) and RESS species (RESS) one after the other. Analyses performed on whole are flagged with COMM (C) and those by RESS species as RESS (R)

## COMM ----

### C.a.1. Alpha diversity (richness) all localities together on full comunity (COMM)
### C.a.2. Alpha diversity by locality on full comunity (COMM)
### C.a.3. Alpha diversity by phylum and locality on full comunity (COMM)

## One tail tests ----
# Boxplot of diversity

# by habitat

COMM_rich<-estimate_richness(COMM, measures="Observed")
COMM_rich$sample_name<-rownames(COMM_rich)
COMM.meta.rich<-left_join(metaCOMM,COMM_rich,by="sample_name")

# How many species by localities ?

#### C.a.1. Alpha diversity (richness) all localities together ----
alpha_tests.new<-function(ps, levels){
  
  list.checks<-list()
  
  rich<-estimate_richness(ps, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  
  # Check for identical values
  if (length(unique(meta.rich$Observed)) == 1) {
    message("All 'Observed' values are identical. Skipping normality test.")
    shapiro <- list(p.value = 0)  # Force non-normal result to trigger Wilcoxon
  } else {
    shapiro <- shapiro.test(meta.rich$Observed)
  }
  
  homogeneity_habitat<-leveneTest( 
    Observed~Type_habitat_broad,data=meta.rich)
  
  list.checks<-list(shapiro,homogeneity_habitat)#,homogeneity_region)
  
  names(list.checks)<-c("Normality","Homogeneity_by_habitat")
  
  print(list.checks)
  
  meta.rich$Type_habitat_broad<-  factor(meta.rich$Type_habitat_broad,
                                         levels = levels)
  # Extract single values for checks
  normality_p_value <- list.checks$Normality$p.value
  
  homogeneity_habitat_p_value <- list.checks$Homogeneity_by_habitat[1, "Pr(>F)"] 
  
  
  # Test by habitat 
  
  if(normality_p_value > 0.05 & homogeneity_habitat_p_value> 0.05) {
    
    test.result.habitat<-t.test(Observed~Type_habitat_broad, data=meta.rich, 
                                paired = FALSE,
                                alternative=c("less"))
    
    pvalue<-test.result.habitat$p.value
    stat<-test.result.habitat$statistic
    test<-"t.test"
    variable.tested<-"Habitat"
    
    result.habitat<-data.frame(pvalue,stat,test,variable.tested)
    
    
  } else {
    
    test.result.habitat<-wilcox.test(Observed~Type_habitat_broad, data=meta.rich, paired = FALSE, exact=FALSE, alternative=c("less"))
    
    pvalue<-test.result.habitat$p.value
    stat<-test.result.habitat$statistic
    test<-"wilcox"
    variable.tested<-"Habitat"
    
    result.habitat<-data.frame(pvalue,stat,test,variable.tested)
    
  } 
  return(result.habitat)
}

# One tail test to test if ricness is lower in ports than in natural habitats 

C.alpha.test.sp.greater.artificial<-alpha_tests.new(COMM,levels=c("Artificial","Natural")) # is richness lower in ports

print(C.alpha.test.sp.greater.artificial)


## Two tailed test ----
### to check different species richness between hapitats across localities and phyla
alpha_tests.two.tail<-function(ps){
  
  list.checks<-list()
  
  rich<-estimate_richness(ps, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  
  # Check for identical values
  if (length(unique(meta.rich$Observed)) == 1) {
    message("All 'Observed' values are identical. Skipping normality test.")
    shapiro <- list(p.value = 0)  # Force non-normal result to trigger Wilcoxon
  } else {
    shapiro <- shapiro.test(meta.rich$Observed)
  }
  
  homogeneity_habitat<-leveneTest( 
    Observed~Type_habitat_broad,data=meta.rich)
  
  list.checks<-list(shapiro,homogeneity_habitat)#,homogeneity_region)
  
  names(list.checks)<-c("Normality","Homogeneity_by_habitat")
  
  print(list.checks)
  
  # Extract single values for checks
  normality_p_value <- list.checks$Normality$p.value
  
  homogeneity_habitat_p_value <- list.checks$Homogeneity_by_habitat[1, "Pr(>F)"] 
  
  
  # Test by habitat 
  
  if(normality_p_value > 0.05 & homogeneity_habitat_p_value> 0.05) {
    
    test.result.habitat<-t.test(Observed~Type_habitat_broad, data=meta.rich)
    
    pvalue<-test.result.habitat$p.value
    stat<-test.result.habitat$statistic
    test<-"t.test"
    variable.tested<-"Habitat"
    
    result.habitat<-data.frame(pvalue,stat,test,variable.tested)
    
    
  } else {
    
    test.result.habitat<-wilcox.test(Observed~Type_habitat_broad, data=meta.rich, paired = FALSE, exact=FALSE)
    
    pvalue<-test.result.habitat$p.value
    stat<-test.result.habitat$statistic
    test<-"wilcox"
    variable.tested<-"Habitat"
    
    result.habitat<-data.frame(pvalue,stat,test,variable.tested)
    
  } 
  return(result.habitat)
}

# First perform two tail test 
C.alpha.test.sp<-alpha_tests.two.tail(COMM)  
rownames(C.alpha.test.sp)<-"All locations"
print(C.alpha.test.sp)


#### C.a.2. Alpha diversity by locality ----

# Greater richness in Natural habitats 

unique_localities<-unique(as.data.frame(sample_data(COMM))$Locality)

# Different richness inside and outside ports

C.result.locality<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(COMM, Locality == locality)
  
  result.test<-alpha_tests.two.tail(ps.sub)
  
  rownames(result.test)<-locality
  
  C.result.locality<-rbind(C.result.locality,result.test)
  
}


print(C.result.locality)

# save

C.result.alpha.two.tail<-rbind(C.alpha.test.sp,
                               C.result.locality)


#### C.a.3a. Alpha diversity by phylum  ----
unique_phylum<-unique(as.data.frame(tax_table(COMM))$Phylum_Worms_Accepted)

C.alpha.test.phyla<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(COMM, Phylum_Worms_Accepted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.two.tail(ps.sub.phyl)
  
  result.test$Phlyum<-phylum
  
  C.alpha.test.phyla<-rbind(C.alpha.test.phyla,result.test)
}

C.alpha.test.phyla

# Make multiple plots for localities and for phylum showing richness 

plot.alpha.multiple<-function(ps.sub){
  
  COMM_rich<-estimate_richness(ps.sub, measures="Observed")
  COMM_rich$sample_name<-rownames(COMM_rich)
  COMM.meta.rich<-left_join(metaCOMM,COMM_rich,by="sample_name")
  
  ggplot(COMM.meta.rich,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                            color=Type_habitat_broad))+
    geom_boxplot()+
    geom_point()+
    xlab("")+ylab("")+
    # labs(title = "Richness of Species by habitat type")+
    # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
    #                    label.y = 100)+   # show only *, **, etc.
    #hide.ns = TRUE) +     # hide 'ns' when not significant
    ggtitle(as.data.frame(tax_table(ps.sub))$Phylum_Worms_Accepted)+
    scale_fill_manual(values=c("darkblue","yellow2"))+
    scale_color_manual(values=c("dodgerblue","gold1"))+
    stat_compare_means(method = "t.test",
                       label = "p.signif")+   # show only *, **, etc.
    
    #facet_grid(~Locality, scale="free_x")+
    # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
    # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
    fig+
    theme(axis.text.x = element_blank(),
          legend.position = "top",
          legend.text=element_text(size=9),
          legend.title=element_text(size=7),
          axis.titile.x = element_text(size=7),
          axis.title.y = element_text(size=7))
  
}




#  localities
theme_set(theme_minimal())

plot.rich.loc<-list()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(COMM, Locality == locality)
  plot<-plot.alpha.multiple(ps.sub)
  
  plot.rich.loc[locality]<-list(plot)
}


plot.rich.loc<-ggarrange(plot.rich.loc$Ajaccio,
                         plot.rich.loc$Bastia,
                         plot.rich.loc$Bonifaccio,
                         plot.rich.loc$Calvi,
                         plot.rich.loc$Nice,
                         plot.rich.loc$`Port-la-Nouvell`,
                         plot.rich.loc$`Seyne-sur-mer`,
                         plot.rich.loc$Sete_1,
                         plot.rich.loc$Sete_2,
                         plot.rich.loc$`Saint-Tropez`,
                         plot.rich.loc$`Porto-Vecchio`,
                         plot.rich.loc$Toulon,
                         nrow=1,ncol=12, common.legend = TRUE)


# Phyla 
theme_set(theme_minimal())

plot.rich.phyla<-list()

for (phylum in unique_phylum){
  
  ps.sub<-subset_taxa(COMM, Phylum_Worms_Accepted == phylum) #repeat the same loop by changing everytime the subset 
  plot<-plot.alpha.multiple(ps.sub)
  
  plot.rich.phyla[phylum]<-list(plot)
}

# Keep only those phyla with significantly different richness
plot.rich.phyla.comb<-ggarrange(plot.rich.phyla$Arthropoda,
                                plot.rich.phyla$Cnidaria,
                                plot.rich.phyla$Vertebrata,
                                plot.rich.phyla$Bryozoa,
                                plot.rich.phyla$Annelida,
                                plot.rich.phyla$Tunicata,
                                plot.rich.phyla$Porifera,
                                plot.rich.phyla$Nemertea,
                                plot.rich.phyla$Rotifera,
                                nrow=1,ncol=9, common.legend = TRUE)

# plot.rich.phyla$Echinodermata,
# plot.rich.phyla$Mollusca,
# plot.rich.phyla$Entoprocta,
# plot.rich.phyla$Porifera,
# plot.rich.phyla$Xenacoelomorpha,
# plot.rich.phyla$Platyhelminthes,
# plot.rich.phyla$Chaetognatha,
# plot.rich.phyla$Phoronida,
# plot.rich.phyla$Nematoda,
# plot.rich.phyla$Placozoa,
# 
# ggsave(plot=plot.rich.phyla.comb, "Whole_community/boxplot_alpha_phyla.pdf",
#        width = 13, height = 6, dpi=300)


## RESS (R) ----


### a.1. Alpha diversity (richness) all localities together on RESS species 
### a.2. Alpha diversity by locality on RESS species 
### a.3. Alpha diversity by phylum and locality on RESS species 

# Boxplot of diversity

# by habitat

rich<-estimate_richness(RESS, measures="Observed")
rich$sample_name<-rownames(rich)
meta<-sample_data(RESS) %>% data.frame()
meta.rich<-left_join(meta,rich,by="sample_name")

#### R.a.1. Alpha diversity (richness) all localities together ----

# tEST IF RICHNESS IS LOWER IN PORTS
alpha.test.sp.greater.artificial<-alpha_tests.new(RESS,levels=c("Artificial","Natural"))

print(alpha.test.sp.greater.artificial)

#### R.a.2. Alpha diversity by locality ----

# Greater richness in Natural habitats 

# Different richness in Artificial habitats

result.locality.RESS.lower.Artificial<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(RESS, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Artificial","Natural"))
  rownames(result.test)<-locality
  
  result.locality.RESS.lower.Artificial<-rbind(result.locality.RESS.lower.Artificial,result.test)
  
}

print(result.locality.RESS.Artificial)


# Plot richness across localities 
theme_set(theme_minimal())

RESS.richness<-ggplot(meta.rich,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                               color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  # stat_compare_means(method = "t.test",label = "p.signif",   # show only *, **, etc.
  # hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

RESS.richness


RESS.richness.allLocs<-ggplot(meta.rich,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                       color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
                     label.y = 100)+   # show only *, **, etc.
  #hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  #facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

RESS.richness.allLocs


## b. RESS communities in ports are dissimilar from those outside ----
## Anlyses performed on COMM and RESS community 

## COMM ----

COMM_PCoA<-ordinate(COMM, method = "PCoA", distance="jaccard") 

sites_PCoA.COMM<-as.data.frame(COMM_PCoA$vectors)
sites_PCoA.COMM$sample_name<-rownames(sites_PCoA.COMM)

sites_PCoA.COMM<-left_join(metaCOMM,sites_PCoA.COMM, by="sample_name")


# Plot PCoA 

theme_set(theme_minimal())

COMM_pcoa.plot<-plot_ordination(COMM, COMM_PCoA, color= "Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
  scale_color_manual(values = c("darkblue","yellow3"))+ 
  labs(title = paste0(" PCoA - COMM"))+
  ggConvexHull:: geom_convexhull(data=sites_PCoA.COMM,aes(x = Axis.1, y = Axis.2, color=Type_habitat_broad,
                                                    fill=Type_habitat_broad),alpha=0.0)+fig

COMM_pcoa.plot   


# Dispersio plot:

centroid.COMM<-sites_PCoA.COMM%>% group_by(Type_habitat_broad)%>%summarise(centroid_Axis.1=mean(Axis.1),centroid_Axis.2=mean(Axis.2))
sites_PCoA.COMM <-  left_join(sites_PCoA.COMM, centroid.COMM, by="Type_habitat_broad")


COMM_pcoa.plot.centroids<-ggplot(data=sites_PCoA.COMM, aes(Axis.1, Axis.2, color=Type_habitat_broad))+
  geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
  geom_point(data=centroid.COMM, aes(x=centroid_Axis.1,y=centroid_Axis.2, color=Type_habitat_broad))+
  geom_segment(aes(xend = centroid_Axis.1, yend = centroid_Axis.2),
               alpha = 0.5)+
  scale_color_manual(values = c("darkblue","yellow3"))+ 
  labs(title = paste0(" PCoA - COMM"))

### C.b.1. Beta diversity test on all localities together ----

jaccard.COMM<-distance(COMM, method = "jaccard", binary=TRUE)
metaCOMM[,"Type_habitat_broad"]<-as.factor(metaCOMM[,"Type_habitat_broad"])

# check beta dispersion - betadisperser test for homogeneity of variance
dispr.COMM<-betadisper(jaccard.COMM,df$Type_habitat_broad) 
permutest(dispr.COMM) 

# Welch MANOVA - Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)

WDS<-MicEco::WdS.test(dm = jaccard.COMM, f = metaCOMM$Type_habitat_broad, nrep =999) 
print(WDS) 

#### C.b.1.1. Indicator species - COMM ####
run_indicspecies<-function(ps){
  
  df<- sample_data(ps) %>% data.frame() # extract the dataframe
  
  otu<-as.data.frame(otu_table(ps))
  
  otu$ids<-rownames(otu)
  
  tax<-as.data.frame(tax_table(ps))
  
  tax$ids<-rownames(tax)
  
  tax<-tax[,c("FARTA","ids")]
  
  otu<-left_join(otu,tax,by="ids") %>% column_to_rownames(var="FARTA")
  
  otu<-select(otu,-ids)
  
  otu<-t(otu) # indicspecies requires sites as rows and species as columns
  
  groups<-df$Type_habitat_broad
  
  indicspecs<-multipatt(otu, groups, 
                        control = how(nperm=999))
  return(indicspecs)
}

set.seed(1234)
indicspecs.COMM<-run_indicspecies(COMM)
indicspecs.lines.COMM<-capture.output(print(summary(indicspecs.COMM))) #57 species

# Classification of indicator species 

indicspecs.COMM.sign<-indicspecs.COMM$sign %>% subset(p.value<=0.05)
indicspecs.COMM.sign$FARTA<-rownames(indicspecs.COMM.sign)
COMM.tax.indicator<- as.data.frame(tax_table(COMM))%>% subset(FARTA %in% rownames(indicspecs.COMM.sign))
COMM.tax.indicator<-left_join(indicspecs.COMM.sign, COMM.tax.indicator, by="FARTA")
COMM.tax.indicator<-COMM.tax.indicator[,c(1:13,38)]

## Add indicator to the PCOA plot 

COMM_indic.axis.top<-indicspecs.COMM.sign %>%
  arrange(desc(stat))

COMM_indic.axis.topN <- subset(COMM_indic.axis.top, s.Natural%in%1)
COMM_indic.axis.topA<-subset(COMM_indic.axis.top, s.Artificial%in%1)
COMM_indic.axis.topA<-COMM_indic.axis.topA[1:20,]
COMM_indic.axis.topNA<-rbind(COMM_indic.axis.topA,COMM_indic.axis.topN)

COMM_indic.axis<-plot_ordination(COMM, COMM_PCoA, type="species")   
COMM_indic.axis<-COMM_indic.axis$data
COMM_indic.axis<-COMM_indic.axis %>% subset(FARTA %in% COMM_indic.axis.top$FARTA)
COMM_indic.axis.top<-COMM_indic.axis %>% subset(FARTA %in% COMM_indic.axis.topNA$FARTA)

COMM_plot.sp.idicator<-plot_ordination(COMM, COMM_PCoA, color= "Type_habitat_broad")+
  xlab("PCoA Axis 1 [5%]")+
  ylab("PCoA Axis 2 [4.3%]")+
  
  geom_point(size=0, alpha=0.1, aes(color=Type_habitat_broad))+
  scale_color_manual(values = c("darkblue","yellow3"))+
  scale_fill_manual(values = c("darkblue","yellow3"))+
  
  ggConvexHull:: geom_convexhull(aes(fill=Type_habitat_broad,
                                     group=Type_habitat_broad),alpha=0.2)+
  
  geom_point(data = COMM_indic.axis.top,inherit.aes = FALSE,
             aes(x = Axis.1, y= Axis.2),
             color = "grey30",
             alpha = 0.3,
             size=3)+
  ggtitle("PCoA - COMM")+
  geom_point(data = COMM_indic.axis,inherit.aes = FALSE,
             aes(x = Axis.1, y= Axis.2),
               color = "grey56",
             alpha = 0.7,
             size=3)+
  
  geom_point(data = COMM_indic.axis.top,inherit.aes = FALSE,
             aes(x = Axis.1, y= Axis.2),
             color = "black",
             alpha = 0.7,
             size=3)+
  
  # geom_segment(data = indic.axis.top,inherit.aes = FALSE,
  #              aes(x = 0, y = 0, xend = Axis.1, yend = Axis.2),
  #              color = "black",
  #              alpha = 0.7,
  #              arrow = arrow(length = unit(0.25, "cm")))+
  ggrepel::geom_text_repel(data=COMM_indic.axis.top,
                           inherit.aes = FALSE,
                           aes(x = Axis.1,
                               y = Axis.2,
                               label=FARTA),color="black",fontface="italic",
                           size = 4,
                           alpha=1,
                           max.overlaps = Inf,       # allow more labels, if needed
                           box.padding = 0.6,          # increase space around text box
                           point.padding = 0.5,        # increase space from data point
                           segment.curvature = 0.2,    # add curvature to leader lines (optional)
                           segment.ncp = 3,        # number of control points (for curvature)
                           segment.angle = 20,         # angle for segment bending
                           force = 2,                  # increase repelling force
                           force_pull = 0.1,           # lower pull towards original location
                           max.time = 2              # give it more time to resolve overlaps
                           # alpha = 1                   # set to max (your alpha=2 was invalid, range is 0–1)
  ) + 
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.position = "top")


COMM_plot.sp.idicator

## RESS ----

RESS_PCoA<-ordinate(RESS, method = "PCoA", distance="jaccard") 

sites_PCoA.RESS<-as.data.frame(RESS_PCoA$vectors)
sites_PCoA.RESS$sample_name<-rownames(sites_PCoA.RESS)

RESS.pcoa.plot<-left_join(df,sites_PCoA.RESS, by="sample_name")


# Plot PCoA 

theme_set(theme_minimal())
    
RESS.pcoa.plot<-plot_ordination(RESS, PCoA, color= "Type_habitat_broad")+
      geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
      scale_color_manual(values = c("darkblue","yellow3"))+ 
      labs(title = paste0(" PCoA "))+
      fig

RESS.pcoa.plot   

# Dispersio plot:

RESS.pcoa.plot
centroid.RESS<-sites_PCoA.RESS%>% group_by(Type_habitat_broad)%>%summarise(centroid_Axis.1=mean(Axis.1),centroid_Axis.2=mean(Axis.2))

sites_PCoA.RESS <-  left_join(sites_PCoA.RESS, centroid.RESS, by="Type_habitat_broad")
RESS_pcoa.plot.centroids<-ggplot(data=sites_PCoA.COMM, aes(Axis.1, Axis.2, color=Type_habitat_broad))+
  geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
  geom_point(data=centroid.RESS, aes(x=centroid_Axis.1,y=centroid_Axis.2, color=Type_habitat_broad))+
  geom_segment(aes(xend = centroid_Axis.1, yend = centroid_Axis.2),
               alpha = 0.5)+
  scale_color_manual(values = c("darkblue","yellow3"))+ 
  labs(title = paste0(" PCoA - RESS"))

### R.b.1. Beta diversity test on all localities together ----

jaccard.RESS<-distance(RESS, method = "jaccard", binary=TRUE)
df<- sample_data(RESS) %>% data.frame() # extract the dataframe
df[,"Type_habitat_broad"]<-as.factor(df[,"Type_habitat_broad"])

# check beta dispersion - betadisperser test for homogeneity of variance
dispr.RESS<-betadisper(jaccard.RESS,df$Type_habitat_broad) 
permutest(dispr.RESS) 

# Welch MANOVA - Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)

WDS<-MicEco::WdS.test(dm = jaccard.RESS, f = df$Type_habitat_broad, nrep =999) 
print(WDS)


# c. Biotic homogenization in ports is greater than in natural habitats ----
## Analyses repeated on COMM and LOCAL


## COMM ----

dispr.COMM.dist<-as.data.frame(dispr.COMM$distances)

colnames(dispr.COMM.dist)[1]<-"Jaccard_distance"

dispr.COMM.dist$sample_name<-rownames(dispr.COMM.dist)

metaCOMM.sub<-metaCOMM %>% select (Type_habitat_broad, sample_name) 

dispr.COMM.dist<-left_join(dispr.COMM.dist, metaCOMM.sub, by="sample_name")

rownames(dispr.COMM.dist)<-NULL

dispr.COMM.dist<- dispr.COMM.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.COMM.dist$Jaccard_distance)

dispr.COMM.dist$Type_habitat_broad<-  factor(dispr.COMM.dist$Type_habitat_broad,
                                             levels = c("Artificial", "Natural"))


# is dispersion lower n ports ? Is homogenization greater in ports?
t.test(Jaccard_distance~Type_habitat_broad, data=dispr.COMM.dist, paired=FALSE, alternative=c("less")) #0.0009

# Plot

comm_homogen.plot<-ggplot(dispr.COMM.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                                               fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  ylab("Jaccard distance")+
  xlab("")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))


## RESS ----

dispr.RESS.dist<-as.data.frame(dispr.RESS$distances)

colnames(dispr.RESS.dist)[1]<-"Jaccard_distance"

dispr.RESS.dist$sample_name<-rownames(dispr.RESS.dist)

df.sub<-df %>% select (Type_habitat_broad, sample_name) 

dispr.RESS.dist<-left_join(dispr.RESS.dist, df.sub, by="sample_name")

rownames(dispr.RESS.dist)<-NULL

dispr.RESS.dist<- dispr.RESS.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.RESS.dist$Jaccard_distance)

dispr.RESS.dist$Type_habitat_broad<-  factor(dispr.RESS.dist$Type_habitat_broad,
                                             levels = c("Artificial", "Natural"))

t.test(Jaccard_distance~Type_habitat_broad, data=dispr.RESS.dist, paired=FALSE, alternative=c("less")) # 0.001

# Plot

RESS.homogen.plot<-ggplot(dispr.RESS.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                        fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  ylab("Jaccard distance")+
  xlab("")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))

#d.  Venn diagram : unique and shared species inside and outside ----

# Find all the Native that are present both in the artificial and natural habitats and those that are shared between the two.

###  COMM----
# Combine taxa table otu table and metadata

COMM_tax<-as.data.frame(tax_table(COMM))
COMM_tax$OTU<-rownames(COMM_tax)

COMM_otu<-as.data.frame(t(otu_table(COMM)))
COMM_otu$sample_name<-rownames(COMM_otu)

# Convert the binary tableOTU# Convert the binary table to a two-column dataset
long_otu.COMM <- COMM_otu %>%
  pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  select(sample_name, OTU)

COMM_metaVenn<-sample_data(COMM)%>% data.frame()

COMM_metaVenn<-COMM_metaVenn %>% select(sample_name,Type_habitat_broad)

COMM_meta_long_otu<-left_join(COMM_metaVenn,long_otu.COMM,by="sample_name")

venn_data.COMM <- COMM_meta_long_otu %>%
  group_by(Type_habitat_broad) %>%
  summarise(OTUs = list(unique(OTU)))

# Prepare the OTU lists for the Venn diagram
# venn_list.COMM <- list(
#   Outside = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Natural"][[1]],
#   Inside = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Artificial"][[1]]
# )

venn_list.COMM <- list(
  Natural = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Natural"][[1]],
  Artificial = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Artificial"][[1]]
)

# Create and plot the Venn diagram with shared OTUs included
library(VennDiagram)
venn.plot.COMM <- venn.diagram(
  x = venn_list.COMM,
  category.names = c("Natural", "Artificial"),
  filename = NULL, # Don't save to a file
  output = TRUE,
  col = "black",
  fill = c("yellow2", "blue"),
  alpha = 0.7,
  cex = 2.5,
  fontface = "bold",
  fontfamily = "calibri",
  cat.cex = 2,
  # cat.fontface = "normal",
  cat.fontfamily = "calibri"
)


grid.draw(venn.plot.COMM)


COMM_shared<-intersect(venn_list.COMM$Natural,venn_list.COMM$Artificial)

COMM_unique.artificial<-venn_list.COMM$Artificial[!venn_list.COMM$Artificial %in% venn_list.COMM$Natural]

comm_unique.natural<-venn_list.COMM$Natural[!venn_list.COMM$Natural %in% venn_list.COMM$Artificial]

# combine with taxonomy info

COMM_shared<- COMM_tax %>% subset(OTU %in% COMM_shared )
COMM_shared$Venn<-"Shared"

COMM_unique.artificial<- COMM_tax %>% subset(OTU %in% COMM_unique.artificial )
COMM_unique.artificial$Venn<-"Unique_Artificial"

comm_unique.natural<- COMM_tax %>% subset(OTU %in% comm_unique.natural )
comm_unique.natural$Venn<-"Unique_Natural"

COMM_unique.shared<-rbind(COMM_unique.artificial,comm_unique.natural,COMM_shared) %>% 
  select(FARTA,Venn, Phylum_Worms_Accepted,
         Class_Worms_Accepted,
         Order_Worms_Accepted,
         Family_Worms_Accepted,
         Genus_Worms_Accepted,
         Species_Worms_Accepted)



### RESS ----

# Find all the Native that are present both in the artificial and natural habitats and those that are shared between the two.

# Combine taxa table otu table and metadata
tax<-as.data.frame(tax_table(RESS))
tax$OTU<-rownames(tax)

otu<-as.data.frame(t(otu_table(RESS)))
otu$sample_name<-rownames(otu)

# Convert the binary tableOTU# Convert the binary table to a two-column dataset
long_otu.RESS <- otu %>%
  pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  select(sample_name, OTU)

metaVenn<-sample_data(RESS)%>% data.frame()

metaVenn<-metaVenn %>% select(sample_name,Type_habitat_broad)

meta_long_otu<-left_join(metaVenn,long_otu.RESS,by="sample_name")

venn_data.RESS <- meta_long_otu %>%
  group_by(Type_habitat_broad) %>%
  summarise(OTUs = list(unique(OTU)))

# Prepare the OTU lists for the Venn diagram
venn_list.RESS <- list(
  Natural = venn_data.RESS$OTUs[venn_data.RESS$Type_habitat_broad == "Natural"][[1]],
  Artificial = venn_data.RESS$OTUs[venn_data.RESS$Type_habitat_broad == "Artificial"][[1]]
)
# Create and plot the Venn diagram with shared OTUs included
library(VennDiagram)
venn.plot.RESS <- venn.diagram(
  x = venn_list.RESS,
  category.names = c("Natural", "Artificial"),
  filename = NULL, # Don't save to a file
  output = TRUE,
  col = "black",
  fill = c("yellow2", "blue"),
  alpha = 0.7,
  cex = 2.5,
  fontface = "bold",
  fontfamily = "calibri",
  cat.cex = 2,
  # cat.fontface = "normal",
  cat.fontfamily = "calibri"
)


grid.draw(venn.plot.RESS)


shared<-intersect(venn_list.RESS$Natural,venn_list.RESS$Artificial)

unique.artificial<-venn_list.RESS$Artificial[!venn_list.RESS$Artificial %in% venn_list.RESS$Natural]

unique.natural<-venn_list.RESS$Natural[!venn_list.RESS$Natural %in% venn_list.RESS$Artificial]

# combine with taxonomy info

shared<- tax %>% subset(OTU %in% shared )
shared$Venn<-"Shared"

unique.artificial<- tax %>% subset(OTU %in% unique.artificial )
unique.artificial$Venn<-"Unique_Artificial"

unique.natural<- tax %>% subset(OTU %in% unique.natural )
unique.natural$Venn<-"Unique_Natural"

unique.shared<-rbind(unique.artificial,unique.artificial,shared) %>% 
  select(FARTA,Venn, Phylum_Worms_Accepted,
         Class_Worms_Accepted,
         Order_Worms_Accepted,
         Family_Worms_Accepted,
         Genus_Worms_Accepted,
         Species_Worms_Accepted)



#e. Beta diversity components - COMM -----
### e.1. Result of beta diversity components: test between nestedness and turnover ----

# Generate a table with the mean jaccard distance, nestedness and turnover for each artificial sample vs each natural sample in each location

final.combine.COMM<-data.frame()

for (i in unique_localities) {
  
  loc<-subset_samples(COMM, Locality %in% i)
  meta.loc<-subset(metaCOMM, Locality %in% i)
  
  
  # calculate the beta distance and components within each locality ans subset the object 
  betapair<- beta.pair(as.data.frame(t(otu_table(loc))), index.family="jaccard") 
  beta.nest<-betapair$beta.jne
  beta.tu<-betapair$beta.jtu
  beta.full<-betapair$beta.jac
  
  # select only the artificial and natural samples in two distinct carachters
  artif.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Artificial"]
  natural.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Natural"]
  
  # from the Jaccard distance matrix select as rows the artificial samples and by column the natural ones 
  beta.full.sub<-as.matrix(beta.full) %>% data.frame()
  beta.full.sub<- beta.full.sub[,c(natural.samples)]
  beta.full.sub<-beta.full.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in% artif.samples) 
  rownames(beta.full.sub)<-NULL
  beta.full.sub<- beta.full.sub %>% column_to_rownames(var="sample")
  
  # from the nestedness matrix select as rows the artificial samples and by column the natural ones 
  beta.nest.sub<-as.matrix(beta.nest) %>% data.frame()
  mean(beta.nest)
  mean(beta.nest.sub)
  beta.nest.sub<- beta.nest.sub[,c(natural.samples)]
  beta.nest.sub<-beta.nest.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in%artif.samples ) 
  rownames(beta.nest.sub)<-NULL
  beta.nest.sub<- beta.nest.sub %>% column_to_rownames(var="sample")
  
  # from the turnover matrix select as rows the artificial samples and by column the natural ones 
  beta.tu.sub<-as.matrix(beta.tu) %>% data.frame()
  mean(beta.tu)
  mean(beta.tu.sub)
  beta.tu.sub<- beta.tu.sub[,c(natural.samples)]
  beta.tu.sub<-beta.tu.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in%artif.samples ) 
  rownames(beta.tu.sub)<-NULL
  beta.tu.sub<- beta.tu.sub %>% column_to_rownames(var="sample")
  
  # calculate the mean distance, nestedness and turnover for each artificial sample vs natural samples
  beta.tu.sub$mean_turnover<- rowMeans(beta.tu.sub)
  
  beta.nest.sub$mean_nest<- rowMeans(beta.nest.sub)
  
  beta.full.sub$mean_full<- rowMeans(beta.full.sub)
  
  
  combine<-as.data.frame(cbind(beta.tu.sub$mean_turnover,beta.nest.sub$mean_nest,beta.full.sub$mean_full))
  colnames(combine)<-c("mean_turnover","mean_nestedness","mean_jaccard")
  rownames(combine)<-rownames(beta.tu.sub)
  
  combine$sample<-rownames(combine)
  
  combine<- pivot_longer(combine, -sample, names_to = "Group")
  combine$Locality<-i
  
  final.combine.COMM<-rbind(combine,final.combine.COMM)
  
}


theme_set(theme_minimal())
COMM_beta.components<-ggplot(final.combine.COMM , aes(Group, value, fill=Group))+
  geom_boxplot()+
  ylab("Distance")+xlab("")+
  ggtitle("Beta diversity components - COMM",
          "Artificial vs Natural sites")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))



#f. Traits explain the differences between habitats - COMM ----

## Dissimilarities between communities inside and outside ports are accounted for by different species traits

COMM_tax<-as.data.frame(tax_table(COMM))

COMM_tax %>% group_by(Trait1) %>% summarise(freq=(n()*100)/537) 
COMM_tax %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/537)


COMM_tax.trait2<-subset(COMM_tax, !Class_Worms_Accepted %in% c("Teleostei",
                                                               "Elasmobranchii"))

COMM_tax.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/387)# all mobile species without fishes (benthic and zooplankton)

# Vertical distribution of fish
COMM_tax.fish<-subset(COMM_tax, Class_Worms_Accepted %in% c("Teleostei",
                                                            "Elasmobranchii"))

COMM_tax.fish %>% group_by(Fishbase_depth_broad) %>% summarise(freq=(n()*100)/150)

# Perform dbRDA 

jaccard.COMM

COMM_otu<-as.data.frame(otu_table(COMM)) 

COMM.dbrda<-capscale(jaccard.COMM~Type_habitat_broad +
                       Condition(Locality), # set locality as condition since it not the focus of the analysis 
                     distance="jaccard", data=metaCOMM, 
                     comm = as.data.frame(t(COMM_otu)))


anova.cca(COMM.dbrda)
anova.cca(COMM.dbrda, by ="terms",permutations = 999)
RsquareAdj(COMM.dbrda) 

summ.comm<-summary(COMM.dbrda)
summ.comm$cont$importance # importanc eof axes

# Test if the CAP1 is influenced by the traits 
cap.COMM<-as.data.frame(scores(COMM.dbrda, display = "species"))
cap.COMM$OTU_seed<-rownames(cap.COMM)

# Create a table where to store the results:

COMM_traits.cap1.res<-as.data.frame(matrix(0,nrow=1,ncol=5))
colnames(COMM_traits.cap1.res)<-c("Trait","test","stat","p value")

# Combine with the trait taxa table

trait.COMM<-as.data.frame(tax_table(COMM)) %>% select(OTU_seed, Trait1,
                                                      Trait2_Worms,
                                                      Fishbase_depth_broad,
                                                      Fishbase_Fishing_vulnerability_numeric,
                                                      Phylum_Worms_Accepted,
                                                      Class_Worms_Accepted,
                                                      NINIS_Native,
                                                      FARTA)


cap.COMM<-left_join(cap.COMM,trait.COMM, by="OTU_seed")%>% column_to_rownames(var="OTU_seed")
cap.COMM$OTU_seed<-rownames(cap.COMM)

# Motility trait 
cap.COMM$Trait1<-as.factor(cap.COMM$Trait1)
cap.COMM.1<-subset(cap.COMM, !Trait1%in%"")

shapiro.test(cap.COMM.1$CAP1)

wilcox.test(cap.COMM.1$CAP1~cap.COMM.1$Trait1) 

COMM_traits.cap1.res[1,"test"]<-"wilcox"
COMM_traits.cap1.res[1,"Trait"]<-"Mobility"
COMM_traits.cap1.res[1,"stat"]<-wilcox.test(cap.COMM.1$CAP1~cap.COMM.1$Trait1)$statistic
COMM_traits.cap1.res[1,"p value"]<-wilcox.test(cap.COMM.1$CAP1~cap.COMM.1$Trait1)$p.value

# Ecology Trait

cap.COMM$Trait2_Worms<-as.factor(cap.COMM$Trait2_Worms)
cap.COMM.2<-subset(cap.COMM, !Trait2_Worms%in%"")
wilcox.test(cap.COMM.2$CAP1~cap.COMM.2$Trait2_Worms) # 5.628e-06

COMM_traits.cap1.res[2,"test"]<-"wilcox"
COMM_traits.cap1.res[2,"Trait"]<-"Distribution in water column"
COMM_traits.cap1.res[2,"stat"]<-wilcox.test(cap.COMM.2$CAP1~cap.COMM.2$Trait2_Worms)$statistic
COMM_traits.cap1.res[2,"p value"]<-wilcox.test(cap.COMM.2$CAP1~cap.COMM.2$Trait2_Worms)$p.value

# Ecology for fish
cap.COMM$Fishbase_depth_broad<-as.factor(cap.COMM$Fishbase_depth_broad)

cap.COMM.3<-cap.COMM %>% subset(!Fishbase_depth_broad=="") 
kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad) # 0.03  
dunn.test::dunn.test(x=cap.COMM.3$CAP1, g=cap.COMM.3$Fishbase_depth_broad, method="bh")

COMM_traits.cap1.res[3,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[3,"Trait"]<-"Distribution in water column (Fishes)"
COMM_traits.cap1.res[3,"p value"]<-kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad)$p.value
COMM_traits.cap1.res[3,"stat"]<-kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad)$statistic
colnames(COMM_traits.cap1.res)[5]<-"significant_pairwise_comparison"
COMM_traits.cap1.res[3,"significant_pairwise_comparison"]<-"reef-associated vs bethopelagic,reef-associated vs demersal, slightly pelagic vs benthopelagic p=0.06"


#  PHylum
View(cap.COMM %>% group_by(Phylum_Worms_Accepted)%>% summarise(freq=(n()/536)*100))
cap.COMM$Phylum_Worms_Accepted<-as.factor(cap.COMM$Phylum_Worms_Accepted)
kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Worms_Accepted) # 0.002351

COMM_res.phyl<-FSA::dunnTest(x=cap.COMM$CAP1,g=cap.COMM$Phylum_Worms_Accepted, method = "bh")
COMM_res.phyl<-COMM_res.phyl$res
COMM_res.phyl<-COMM_res.phyl$Comparison[COMM_res.phyl$P.adj<=0.05] # "Arthropoda - Tunicata"    "Echinodermata - Tunicata" "Arthropoda - Vertebrata" 

COMM_traits.cap1.res[4,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[4,"Trait"]<-"Phylum"
COMM_traits.cap1.res[4,"stat"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Worms_Accepted)$statistic
COMM_traits.cap1.res[4,"p value"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Worms_Accepted)$p.value
COMM_traits.cap1.res[4,"significant_pairwise_comparison"]<-c("Arthropoda - Tunicata,Echinodermata - Tunicata,Arthropoda - Vertebrata" )


#Class

cap.COMM$Class_Worms_Accepted<-as.factor(cap.COMM$Class_Worms_Accepted)
kruskal.test(cap.COMM$CAP1,cap.COMM$Class_Worms_Accepted) #  0.04

COMM_res.class<-FSA::dunnTest(x=cap.COMM$CAP1,g=cap.COMM$Class_Worms_Accepted, method = "bh")
COMM_res.class<-COMM_res.class$res
COMM_res.class<-COMM_res.class$Comparison[COMM_res.class$P.adj<=0.05] # "Ascidiacea - Copepoda" "Copepoda - Teleostei" 

COMM_traits.cap1.res[5,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[5,"Trait"]<-"Class"
COMM_traits.cap1.res[5,"stat"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Class_Worms_Accepted)$statistic
COMM_traits.cap1.res[5,"p value"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Class_Worms_Accepted)$p.value
COMM_traits.cap1.res[5,"significant_pairwise_comparison"]<-c("Ascidiacea - Copepoda, Copepoda - Teleostei" )


# fishing vulnerability 

cap.COMM.6<-cap.COMM %>% subset(!Fishbase_Fishing_vulnerability_numeric=="") 
cap.COMM.6$Fishbase_Fishing_vulnerability_numeric<-as.numeric(cap.COMM.6$Fishbase_Fishing_vulnerability_numeric)
COMM_cor.tets<-cor.test(cap.COMM.6$CAP1,cap.COMM.6$Fishbase_Fishing_vulnerability_numeric, method = "spearman") # 0.06

COMM_traits.cap1.res[6,"test"]<-"Spearman_correlation"
COMM_traits.cap1.res[6,"Trait"]<-"Fishing_vulnerability"
COMM_traits.cap1.res[6,"stat"]<-COMM_cor.tets$statistic
COMM_traits.cap1.res[6,"p value"]<-COMM_cor.tets$p.value


COMM_fv.plot<-ggplot(cap.COMM.6, aes(CAP1,Fishbase_Fishing_vulnerability_numeric))+
  geom_point()+
  geom_smooth(method="lm", color="red")

# Status

cap.COMM$NINIS_Native<-as.character(cap.COMM$NINIS_Native)
cap.COMM$NINIS_Native[cap.COMM$NINIS_Native!="NIS"]<-"RESS" # resident species
cap.COMM$NINIS_Native[ is.na(cap.COMM$NINIS_Native)]<-"RESS" 
unique(cap.COMM$NINIS_Native)

cap.COMM$NINIS_Native<-as.factor(cap.COMM$NINIS_Native)

shapiro.test(cap.COMM$CAP1)
wilcox.test(cap.COMM$CAP1~cap.COMM$NINIS_Native) 

COMM_traits.cap1.res[7,"test"]<-"Wilcox"
COMM_traits.cap1.res[7,"Trait"]<-"NIS"
COMM_traits.cap1.res[7,"stat"]<-wilcox.test(cap.COMM$CAP1~cap.COMM$NINIS_Native)$statistic
COMM_traits.cap1.res[7,"p value"]<-wilcox.test(cap.COMM$CAP1~cap.COMM$NINIS_Native)$p.value
COMM_traits.cap1.res[7,"significant_pairwise_comparison"]<-NA

#### Plot dbRDA 

### Generate Figure 6 - COMM 

COMM_indic.axis.top<-indicspecs.COMM.sign %>%
  arrange(desc(stat))

COMM_indic.axis.dbRDA<-cap.COMM %>% subset(FARTA %in% COMM_indic.axis.top$FARTA)
COMM_indic.axis.dbRDA.top<-COMM_indic.axis.dbRDA %>% subset(FARTA %in% COMM_indic.axis.topNA$FARTA)

# plot indicator species on dbRDA
dbRDA.COMM.species<-ggplot() +
  
  geom_point(data = COMM_indic.axis.dbRDA,inherit.aes = FALSE,
             aes(x = CAP1, y= MDS1),
             color = "grey30",
             alpha = 0.3,
             size=3)+
  
  geom_point(data = COMM_indic.axis.dbRDA.top,inherit.aes = FALSE,
             aes(x = CAP1, y= MDS1),
             color = "grey5",
             alpha = 0.7,
             size=3)+
  ggrepel::geom_text_repel(data=COMM_indic.axis.dbRDA.top,
                           inherit.aes = FALSE,
                           aes(x = CAP1,
                               y = MDS1,
                               label=FARTA),color="black",fontface="italic", size=3)+
  
  xlab("CAP1")+ylab("MDS1")+
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))

dbRDA.COMM.species

# Generate bubblu plot with funcitonal traits 
cap.COMM.summary.Trait1<-cap.COMM.1 %>% group_by(Trait1)%>% summarise(mean_CAP1=mean(CAP1),
                                                                      mean_MDS1=mean(MDS1))
cap.COMM.summary.Trait1$Trait1<-as.character(cap.COMM.summary.Trait1$Trait1)
cap.COMM.summary.Trait1$Trait1[cap.COMM.summary.Trait1$Trait1=="Mobile"] <- "mobile"
cap.COMM.summary.Trait1$Trait1[cap.COMM.summary.Trait1$Trait1=="Sessile"] <- "sessile"

cap.COMM.summary.Trait2<-cap.COMM.2 %>% group_by(Trait2_Worms)%>% summarise(mean_CAP1=mean(CAP1),
                                                                            mean_MDS1=mean(MDS1))

cap.COMM.summary.Fish_depth<-cap.COMM.3 %>% group_by(Fishbase_depth_broad)%>% summarise(mean_CAP1=mean(CAP1),
                                                                                        mean_MDS1=mean(MDS1))
cap.COMM.summary.Phylum<-cap.COMM %>% group_by(Phylum_Worms_Accepted)%>%
  subset(Phylum_Worms_Accepted%in% c("Vertebrata","Arthropoda","Echinodermata", "Tunicata"))%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

cap.COMM.summary.Class<-cap.COMM %>% group_by(Class_Worms_Accepted)%>%
  subset(Class_Worms_Accepted%in% c("Ascidiacea","Teleostei","Copepoda"))%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

cap.COMM.summary.nis<-cap.COMM %>% group_by(NINIS_Native)%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))


# By each dataset: buble plot
bubble_plot<-ggplot()+
  
  geom_hline(yintercept = 0, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Trait1,
             aes(x = mean_CAP1,
                 y = 0),
             color = "darkgoldenrod3",
             alpha=5,
             size = 3)+ 
  # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Trait1,aes(
    x = mean_CAP1,
    y = 0,
    label = Trait1),
    size = 4, segment.color = "darkgoldenrod3",
    color="darkgoldenrod3",
    angle=45) +
  
  # Depth non fish 
  geom_hline(yintercept = -0.03, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Trait2, 
             aes(x = mean_CAP1,
                 y = -0.03),
             color = "steelblue",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Trait2,
                           aes(x = mean_CAP1,
                               y = -0.03,
                               label = Trait2_Worms),
                           size = 4, segment.color = "steelblue",
                           color="steelblue",
                           angle=45) +
  
  
  
  #Depth - fish
  geom_hline(yintercept = -0.06, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Fish_depth, 
             aes(x = mean_CAP1,
                 y = -0.06),
             color = "#DC143C",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Fish_depth,
                           aes(x = mean_CAP1,
                               y = -0.06,
                               label = Fishbase_depth_broad),
                           size = 4, segment.color = "#DC143C",
                           color="#DC143C",
                           angle=45) +
  
  # Phylum
  geom_hline(yintercept = -0.09, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Phylum, 
             aes(x = mean_CAP1,
                 y = -0.09),
             color = "grey5",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Phylum,
                           aes(x = mean_CAP1,
                               y = -0.09,
                               label = Phylum_Worms_Accepted),
                           size = 4, segment.color = "grey5",
                           color="grey5",
                           angle=45) +
  
  
  # Class
  geom_hline(yintercept = -0.12, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Class, 
             aes(x = mean_CAP1,
                 y = -0.12),
             color = "grey35",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Class,
                           aes(x = mean_CAP1,
                               y = -0.12,
                               label = Class_Worms_Accepted),
                           size = 4, segment.color = "grey35",
                           color="grey35",
                           angle=45) +
  
  
  # Status
  geom_hline(yintercept = -0.15, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.nis, 
             aes(x = mean_CAP1,
                 y = -0.15),
             color = "#117733",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.nis,
                           aes(x = mean_CAP1,
                               y = -0.15,
                               label = NINIS_Native),
                           size = 4, segment.color = "#117733",
                           color="#117733",
                           angle=45) +
  
  # Remove y-axis (since everything is on y = 0)
  xlab("CAP1")+
  theme_minimal() +
  theme(axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=15),
        axis.title.x = element_text(size=17))

# ggsave(plot=bubble_plot, "Whole_community/Figure_6B.pdf",
#        width=12,height = 6, dpi = 300)




#e. Beta diversity components - RESS -----
### e.2. Result of beta diversity components: test between nestedness and turnover ----

# Generate a table with the mean jaccard distance, nestedness and turnover for each artificial sample vs each natural sample in each location

final.combine.RESS<-data.frame()

for (i in unique_localities) {
  
  loc<-subset_samples(RESS, Locality %in% i)
  meta.loc<-subset(df, Locality %in% i)

    
  # calculate the beta distance and components within each locality ans subset the object 
  betapair<- beta.pair(as.data.frame(t(otu_table(loc))), index.family="jaccard") 
  beta.nest<-betapair$beta.jne
  beta.tu<-betapair$beta.jtu
  beta.full<-betapair$beta.jac
  
  # select only the artificial and natural samples in two distinct carachters
  artif.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Artificial"]
  natural.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Natural"]
  
  # from the Jaccard distance matrix select as rows the artificial samples and by column the natural ones 
  beta.full.sub<-as.matrix(beta.full) %>% data.frame()
  beta.full.sub<- beta.full.sub[,c(natural.samples)]
  beta.full.sub<-beta.full.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in% artif.samples) 
  rownames(beta.full.sub)<-NULL
  beta.full.sub<- beta.full.sub %>% column_to_rownames(var="sample")
  
  # from the nestedness matrix select as rows the artificial samples and by column the natural ones 
  beta.nest.sub<-as.matrix(beta.nest) %>% data.frame()
  mean(beta.nest)
  mean(beta.nest.sub)
  beta.nest.sub<- beta.nest.sub[,c(natural.samples)]
  beta.nest.sub<-beta.nest.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in%artif.samples ) 
  rownames(beta.nest.sub)<-NULL
  beta.nest.sub<- beta.nest.sub %>% column_to_rownames(var="sample")
  
  # from the turnover matrix select as rows the artificial samples and by column the natural ones 
  beta.tu.sub<-as.matrix(beta.tu) %>% data.frame()
  mean(beta.tu)
  mean(beta.tu.sub)
  beta.tu.sub<- beta.tu.sub[,c(natural.samples)]
  beta.tu.sub<-beta.tu.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in%artif.samples ) 
  rownames(beta.tu.sub)<-NULL
  beta.tu.sub<- beta.tu.sub %>% column_to_rownames(var="sample")
  
  # calculate the mean distance, nestedness and turnover for each artificial sample vs natural samples
  beta.tu.sub$mean_turnover<- rowMeans(beta.tu.sub)
  
  beta.nest.sub$mean_nest<- rowMeans(beta.nest.sub)
  
  beta.full.sub$mean_full<- rowMeans(beta.full.sub)
  
  
  combine<-as.data.frame(cbind(beta.tu.sub$mean_turnover,beta.nest.sub$mean_nest,beta.full.sub$mean_full))
  colnames(combine)<-c("mean_turnover","mean_nestedness","mean_jaccard")
  rownames(combine)<-rownames(beta.tu.sub)
  
  combine$sample<-rownames(combine)
  
  combine<- pivot_longer(combine, -sample, names_to = "Group")
  combine$Locality<-i
  
  final.combine.RESS<-rbind(combine,final.combine.RESS)
  
}


theme_set(theme_minimal())
beta.components<-ggplot(final.combine.RESS , aes(Group, value, fill=Group))+
  geom_boxplot()+
  ylab("Distance")+xlab("")+
  ggtitle("Beta diversity components",
          "Artificial vs Natural sites")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))



#f. Traits explain the differences between habitats - RESS ----
## Dissimilarities between communities inside and outside ports are accounted for by different species traits

tax<-as.data.frame(tax_table(RESS))

tax %>% group_by(Trait1) %>% summarise(freq=(n()*100)/495)
tax %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/495)


tax.trait2<-subset(tax, !Class_Worms_Accepted %in% c("Teleostei",
                                                    "Elasmobranchii"))

tax.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/347)

## what percent of mobile species is benthic or pelagic?

tax.trait2.mob<-subset(tax.trait2, Trait1 %in% "Mobile")
tax.trait2.mob %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/184)


# Vertical distribution of fish
tax.fish<-subset(tax, Class_Worms_Accepted %in% c("Teleostei",
                                                              "Elasmobranchii"))

tax.fish %>% group_by(Fishbase_depth_broad) %>% summarise(freq=(n()*100)/148)


# Perform dbRDA 

jaccard.RESS

community<-as.data.frame(otu_table(RESS)) 

RESS.dbrda<-capscale(jaccard.RESS~Type_habitat_broad +
                        Condition(Locality), # set locality as condition since it not the focus of the analysis 
                      distance="jaccard", data=df, 
                      comm = as.data.frame(t(community)))


anova.cca(RESS.dbrda)
anova.cca(RESS.dbrda, by ="terms",permutations = 999)
RsquareAdj(RESS.dbrda)

# Test if the CAP1 is influenced by the traits 
cap.RESS<-as.data.frame(scores(RESS.dbrda, display = "species"))
cap.RESS$OTU_seed<-rownames(cap.RESS)

# Create a table where to store the results:

traits.cap1.res<-as.data.frame(matrix(0,nrow=1,ncol=5))
colnames(traits.cap1.res)<-c("Trait","test","stat","p value")

# Combine with the trait taxa table

trait.RESS<-as.data.frame(tax_table(RESS)) %>% select(OTU_seed, Trait1,
                                                        Trait2_Worms,
                                                        Fishbase_depth_broad,
                                                        Fishbase_Fishing_vulnerability_numeric,
                                                        Phylum_Worms_Accepted,
                                                        Class_Worms_Accepted)
                                                        #FARTA)


cap.RESS<-left_join(cap.RESS,trait.RESS, by="OTU_seed")%>% column_to_rownames(var="OTU_seed")
cap.RESS$OTU_seed<-rownames(cap.RESS)

# Motility trait 

cap.RESS$Trait1<-as.factor(cap.RESS$Trait1)
cap.RESS.1<-subset(cap.RESS, !Trait1%in%"")
  
shapiro.test(cap.RESS.1$CAP1)

wilcox.test(cap.RESS.1$CAP1~cap.RESS.1$Trait1) # 0.04

traits.cap1.res[1,"test"]<-"wilcox"
traits.cap1.res[1,"Trait"]<-"Mobility"
traits.cap1.res[1,"stat"]<-wilcox.test(cap.noNIS.1$CAP1~cap.noNIS.1$Trait1)$statistic
traits.cap1.res[1,"p value"]<-wilcox.test(cap.noNIS.1$CAP1~cap.noNIS.1$Trait1)$p.value

# Ecology 
cap.RESS$Trait2_Worms<-as.factor(cap.RESS$Trait2_Worms)
cap.RESS.2<-subset(cap.RESS, !Trait2_Worms%in%"")
wilcox.test(cap.RESS.2$CAP1~cap.RESS.2$Trait2_Worms) # 1.037e-07

traits.cap1.res[2,"test"]<-"wilcox"
traits.cap1.res[2,"Trait"]<-"Distribution in water column"
traits.cap1.res[2,"stat"]<-wilcox.test(cap.RESS.2$CAP1~cap.RESS.2$Trait2_Worms)$statistic
traits.cap1.res[2,"p value"]<-wilcox.test(cap.RESS.2$CAP1~cap.RESS.2$Trait2_Worms)$p.value

# Ecology Fish 
cap.RESS$Fishbase_depth_broad<-as.factor(cap.RESS$Fishbase_depth_broad)

cap.RESS.3<-cap.RESS %>% subset(!Fishbase_depth_broad=="") 
kruskal.test(cap.RESS.3$CAP1~cap.RESS.3$Fishbase_depth_broad) # 0.02 
dunn.test::dunn.test(x=cap.RESS.3$CAP1, g=cap.RESS.3$Fishbase_depth_broad, method="bh")

traits.cap1.res[3,"test"]<-"Kruskal-Wallis"
traits.cap1.res[3,"Trait"]<-"Distribution in water column (Fishes)"
traits.cap1.res[3,"stat"]<-kruskal.test(cap.RESS.3$CAP1~cap.RESS.3$Fishbase_depth_broad)$p.value
traits.cap1.res[3,"p value"]<-kruskal.test(cap.RESS.3$CAP1~cap.RESS.3$Fishbase_depth_broad)$statistic
colnames(traits.cap1.res)[5]<-"significant_pairwise_comparison"
traits.cap1.res[3,"significant_pairwise_comparison"]<-"pelagic vs benthopelagic, reef-associated vs bethopelagic,reef-associated vs demersal, slightly pelagic vs demersal p=0.06"


# Phylum
cap.RESS$Phylum_Worms_Accepted<-as.factor(cap.RESS$Phylum_Worms_Accepted)
kruskal.test(cap.RESS$CAP1,cap.RESS$Phylum_Worms_Accepted) # 0.001

res.phyl<-FSA::dunnTest(x=cap.RESS$CAP1,g=cap.RESS$Phylum_Worms_Accepted, method = "bh")
res.phyl<-res.phyl$res
res.phyl<-res.phyl$Comparison[res.phyl$P.adj<=0.05] # Arthropoda vs Vertebrata

traits.cap1.res[4,"test"]<-"Kruskal-Wallis"
traits.cap1.res[4,"Trait"]<-"Phylum"
traits.cap1.res[4,"stat"]<-kruskal.test(cap.RESS$CAP1,cap.RESS$Phylum_Worms_Accepted)$statistic
traits.cap1.res[4,"p value"]<-kruskal.test(cap.RESS$CAP1,cap.RESS$Phylum_Worms_Accepted)$p.value
traits.cap1.res[4,"significant_pairwise_comparison"]<-"Arthropoda vs Vertebrata"


# Class
cap.RESS$Class_Worms_Accepted<-as.factor(cap.RESS$Class_Worms_Accepted)
kruskal.test(cap.RESS$CAP1,cap.RESS$Class_Worms_Accepted) # 0.03

res.class<-FSA::dunnTest(x=cap.RESS$CAP1,g=cap.RESS$Class_Worms_Accepted, method = "bh")
res.class<-res.class$res
res.class<-res.class$Comparison[res.class$P.adj<=0.05] # Copepoda vs Teleostei

traits.cap1.res[5,"test"]<-"Kruskal-Wallis"
traits.cap1.res[5,"Trait"]<-"Class"
traits.cap1.res[5,"stat"]<-kruskal.test(cap.RESS$CAP1,cap.RESS$Class_Worms_Accepted)$statistic
traits.cap1.res[5,"p value"]<-kruskal.test(cap.RESS$CAP1,cap.RESS$Class_Worms_Accepted)$p.value
traits.cap1.res[5,"significant_pairwise_comparison"]<-"Copepoda vs Teleostei"


# Fishing vulnerability

cap.RESS.6<-cap.RESS %>% subset(!Fishbase_Fishing_vulnerability_numeric=="") 
cap.RESS.6$Fishbase_Fishing_vulnerability_numeric<-as.numeric(cap.RESS.6$Fishbase_Fishing_vulnerability_numeric)
cor.tets<-cor.test(cap.RESS.6$CAP1,cap.RESS.6$Fishbase_Fishing_vulnerability_numeric, method = "spearman")

traits.cap1.res[6,"test"]<-"Spearman_correlation"
traits.cap1.res[6,"Trait"]<-"Fishing_vulnerability"
traits.cap1.res[6,"stat"]<-cor.tets$statistic
traits.cap1.res[6,"p value"]<-cor.tets$p.value


#### Plot dbRDA

cap.RESS.summary.Trait1<-cap.RESS.1 %>% group_by(Trait1)%>% summarise(mean_CAP1=mean(CAP1),
                                                                      mean_MDS1=mean(MDS1))
cap.RESS.summary.Trait1$Trait1<-as.character(cap.RESS.summary.Trait1$Trait1)
cap.RESS.summary.Trait1$Trait1[cap.RESS.summary.Trait1$Trait1=="Mobile"] <- "mobile"
cap.RESS.summary.Trait1$Trait1[cap.RESS.summary.Trait1$Trait1=="Sessile"] <- "sessile"

cap.RESS.summary.Trait2<-cap.RESS.2 %>% group_by(Trait2_Worms)%>% summarise(mean_CAP1=mean(CAP1),
                                                                              mean_MDS1=mean(MDS1))

cap.RESS.summary.Fish_depth<-cap.RESS.3 %>% group_by(Fishbase_depth_broad)%>% summarise(mean_CAP1=mean(CAP1),
                                                                                          mean_MDS1=mean(MDS1))
cap.RESS.summary.Phylum<-cap.RESS %>% group_by(Phylum_Worms_Accepted)%>%
  subset(Phylum_Worms_Accepted%in% c("Vertebrata","Arthropoda"))%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

dbRDA.RESS.species<-ggplot() +
  
  geom_point(data=cap.RESS.summary.Trait1,
             aes(x = mean_CAP1,
                 y = mean_MDS1,
                 color = Trait1,
                 label=Trait1),
             alpha=5,
             size = 2) +
  ggrepel::geom_text_repel(data=cap.RESS.summary.Trait1,
                           aes(x = mean_CAP1,
                               y = mean_MDS1,
                               label = Trait1,
                               color = Trait1),
                           size = 5,
                           alpha=2) +
  geom_segment(data = cap.RESS.summary.Trait1,
               aes(x = 0, y = 0, xend = mean_CAP1, yend = mean_MDS1,
                   color = Trait1),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  
  geom_point(data=cap.RESS.summary.Trait2,
             aes(x = mean_CAP1,
                 y = mean_MDS1,
                 color = Trait2_Worms,
                 label=Trait2_Worms),
             alpha=5,
             size = 2) +
  ggrepel::geom_text_repel(data=cap.RESS.summary.Trait2,
                           aes(x = mean_CAP1,
                               y = mean_MDS1,
                               label = Trait2_Worms,
                               color = Trait2_Worms),
                           size = 5,
                           alpha=2) +
  geom_segment(data = cap.RESS.summary.Trait2,
               aes(x = 0, y = 0, xend = mean_CAP1, yend = mean_MDS1,
                   color = Trait2_Worms),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  
  geom_point(data=cap.RESS.summary.Fish_depth,
             aes(x = mean_CAP1,
                 y = mean_MDS1,
                 color = Fishbase_depth_broad,
                 label=Fishbase_depth_broad),
             alpha=5,
             size = 2)+
  
  ggrepel::geom_text_repel(data=cap.RESS.summary.Fish_depth,
                           aes(x = mean_CAP1,
                               y = mean_MDS1,
                               label = Fishbase_depth_broad,
                               color = Fishbase_depth_broad),
                           size = 5,
                           alpha=2) +
  geom_segment(data = cap.RESS.summary.Fish_depth,
               aes(x = 0, y = 0, xend = mean_CAP1, yend = mean_MDS1,
                   color = Fishbase_depth_broad),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  ggrepel::geom_text_repel(data=cap.RESS.summary.Phylum,
                           aes(x = mean_CAP1,
                               y = mean_MDS1,
                               label = Phylum_Worms_Accepted,
                               color = Phylum_Worms_Accepted),
                           size = 5,
                           alpha=2) +
  geom_segment(data = cap.RESS.summary.Phylum,
               aes(x = 0, y = 0, xend = mean_CAP1, yend = mean_MDS1,
                   color = Phylum_Worms_Accepted),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  
  scale_color_manual(values = c(
    "grey46",#Arthropoda
    "grey57",#Vertebrata
    "#DC143C",#benthopelagic
    "steelblue",#benthos
    "#FF0000",#demersal
    "#8A2BE2",#mobile
    "firebrick", #pelagic
    "brown1", #reef associated
    "#800080", #sessile
    "lightblue3"))+ #zooplankton
  
  
  xlab("CAP1")+ylab("MDS1")+
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))



# Plot sites and taxa together and colo by mobility
sites.RESS<-as.data.frame(scores(RESS.dbrda, display = "sites"))
sites.RESS$sample_name<-rownames(sites.RESS)
sites.RESS<-left_join(sites.RESS, df[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.RESS<-ggplot() +
  geom_point(data=sites.RESS, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.noNIS.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.RESS,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                      fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
  theme_minimal()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))

combo.dbrda.plot<-ggarrange(dbrda.RESS,dbRDA.RESS.species,nrow = 1, ncol = 2)




# g. Generalities of NIS -----

taxNIS<-as.data.frame(tax_table(NIS))
length(taxNIS$FARTA[taxNIS$NINIS_Native%in%"NIS"]) #43
length(taxNIS$FARTA[taxNIS$Status_Med_Sea%in%"Non-indigenous-Med"]) #34
length(taxNIS$FARTA[taxNIS$Status_Med_Sea%in%"Non-indigenous-non-Med"]) #8
length(taxNIS$FARTA[taxNIS$Status_NIS_table%in%"Cryptogenic"]) #7

# checke the number of RESS assigned to NINIS or to Native
tax<-as.data.frame(tax_table(RESS))
length(tax$FARTA[tax$NINIS_Native%in%"NINIS"]) #120 
length(tax$FARTA[tax$NINIS_Native%in%"Native"]) #326
length(taxNIS$FARTA[!tax$NINIS_Native%in% c("NINIS","Native")]) #49

## h. Proportion of NIS in the community ----

### h.1. barplot with proportion NIS vs RESS inside and outside by locality ----

counts<-as.data.frame(otu_table(combined_taxa_ps_summer.sp))

counts$OTU_seed<-rownames(counts)

tax.all<-as.data.frame(tax_table(combined_taxa_ps_summer.sp))

tax.count<-left_join(counts,tax.all,by="OTU_seed")

rownames(tax.count)<-tax.count$OTU_seed

# write.csv2(tax.count, file="OTU_Taxonomy_table_with_NIS_info_SUMMER_SPECIES_only.csv")

# Obtain the proportions of Non-native species over total species in each site.

get_total<-function(data, colname){
  
  sample.cols <- unlist(lapply(data, is.numeric), use.names = FALSE)  
  
  sample.cols<-colnames(data[ , sample.cols])
  
  total<-as.data.frame(unlist(lapply(data[ , sample.cols], sum)))%>% 
    rownames_to_column(var="sample_name")
  
  colnames(total)[2]<-{{colname}}
  
  return(total)
  
}


total<-get_total(tax.count,"total_number_species")

tax.countNIS.all<-subset(tax.count, NINIS_Native %in% c("NIS")) # here I include all the NIS both from and not the Mediterranean Sea

total.NIS.all<-get_total(tax.countNIS.all,"total_number_NIS")

proportion.table<-list(total,total.NIS.all) %>% reduce(left_join, by="sample_name")

proportion.table<-proportion.table %>%  
  mutate(percentage_NIS_all=
           (proportion.table$total_number_NIS/proportion.table$total_number_species)*100)

# Mean proportion of NIS in the samples

mean(proportion.table$percentage_NIS_all) 
sd(proportion.table$percentage_NIS_all) 

# Bind the sample data to this table
df.NIS<-sample_data(combined_taxa_ps_summer.sp)%>% data.frame()

meta.with.proportions<-left_join(df.NIS,proportion.table,by="sample_name")

rownames(meta.with.proportions)<-meta.with.proportions$sample_name

meta.with.proportions %>% group_by(Locality) %>% summarise(mean.NIS<-mean(percentage_NIS_all),
                                                           sd.NIS<-sd(percentage_NIS_all)) 

# PLot the prcentage of species (Non-Indigenous) by habitat type and locality 

summary_NIS.bysample<-meta.with.proportions %>% group_by(Type_habitat_broad,Locality) %>%
  summarise(percent_NIS_all=percentage_NIS_all,
            sample_name=sample_name) %>% as.data.frame()

rownames(summary_NIS.bysample)<-summary_NIS.bysample$sample_name


summary_NIS<-meta.with.proportions %>% group_by(Type_habitat_broad,Locality) %>% 
  summarise(mean_NIS.ALL=mean(percentage_NIS_all),
            sd_NIS.ALL=sd(percentage_NIS_all)) %>% as.data.frame()

meta.with.proportions %>% group_by(Type_habitat_broad) %>% 
  summarise(mean_NIS.ALL=mean(percentage_NIS_all),
            sd_NIS.ALL=sd(percentage_NIS_all)) %>% as.data.frame()

summary_NIS$Locality[summary_NIS$Locality=="Bonifaccio"]<-"Bonifacio"

theme_set(theme_minimal())

NIS_all<-ggplot(summary_NIS, aes(fill=Type_habitat_broad, colour = Type_habitat_broad))+
  
  geom_bar(aes(Type_habitat_broad,mean_NIS.ALL),stat="identity")+
  
  geom_errorbar(aes(x=Type_habitat_broad,
                    ymin=mean_NIS.ALL-sd_NIS.ALL,
                    ymax=mean_NIS.ALL+sd_NIS.ALL),
                width=0.4,  alpha=0.9, size=1)+
  ylab("% of NIS")+xlab("Type of habitat")+ggtitle("Percentage of NIS in the community","")+
  facet_wrap(~Locality,  scale="free_x")+
  scale_fill_manual(values = c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold2"))+
  theme(strip.background = element_blank(), strip.text= element_text(size=12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme (axis.text.x = element_text(size=15))+
  theme (axis.text.y = element_text(size=15))+
  theme (axis.title.y= element_text(size=16))+
  theme(title = element_text(size=17, vjust = 0.5))+
  theme (legend.position = "none")

NIS_all


### h.2. test on the proportions ----

shapiro.test(meta.with.proportions$percentage_NIS_all)

meta.with.proportions$Type_habitat_broad<-factor(meta.with.proportions$Type_habitat_broad, 
                                                 levels=c("Artificial","Natural"))

test<- wilcox.test(meta.with.proportions$percentage_NIS_all~meta.with.proportions$Type_habitat_broad, alternative="greater") # t = 5.1316, df = 62.749, p-value = 1.498e-06 
p_val<-test$p.value
stat<-test$statistic
res.prop<-data.frame(p_val,stat)
rownames(res.prop)<-"All_localities"


# By locality - Artificial 

unique_localities<-unique(sample_data(NIS)$Locality)

result.data.finAll<-data.frame()

for (i in unique_localities[-5]) {
  
  loc<-subset(meta.with.proportions, Locality == i)
  
  shapiro<-shapiro.test(loc$percentage_NIS_all)
  
  homogeneity_habitat<-leveneTest( 
    percentage_NIS_all~Type_habitat_broad,data=loc)
  
  list.checks<-list(shapiro,homogeneity_habitat)
  
  # Extract single values for checks
  
  normality_p_value <- as.numeric(list.checks[[1]][2][1])
  
  homogeneity_habitat_p_value <- as.numeric(list.checks[[2]][3][[1]][1])
  
  loc$Type_habitat_broad<-factor(loc$Type_habitat_broad, 
                                 levels=c("Artificial","Natural"))
  
  if(normality_p_value > 0.05 & homogeneity_habitat_p_value> 0.05) {
    
    
    test<-  t.test(loc$percentage_NIS_all~loc$Type_habitat_broad, alternative="greater")
    p_val<-test$p.value
    stat<-test$statistic
    result.data<-data.frame(p_val,stat,i)
    rownames(result.data)<-i
    result.data.finAll<-rbind(result.data,result.data.finAll)
  } else {
    
    test<-  wilcox.test(loc$percentage_NIS_all~loc$Type_habitat_broad, alternative="greater")
    p_val<-test$p.value
    stat<-test$statistic
    result.data<-data.frame(p_val,stat,i)
    rownames(result.data)<-i
    result.data.finAll<-rbind(result.data,result.data.finAll)
  }
}

res.prop$i<-"All_localities"

# Localities - Natural 
result.data.fin.Nat<-data.frame()

for (i in unique_localities[-5]) {
  
  loc<-subset(meta.with.proportions, Locality == i)
  
  shapiro<-shapiro.test(loc$percentage_NIS_all)
  
  homogeneity_habitat<-leveneTest( 
    percentage_NIS_all~Type_habitat_broad,data=loc)
  
  list.checks<-list(shapiro,homogeneity_habitat)
  
  # Extract single values for checks
  
  normality_p_value <- as.numeric(list.checks[[1]][2][1])
  
  homogeneity_habitat_p_value <- as.numeric(list.checks[[2]][3][[1]][1])
  
  loc$Type_habitat_broad<-factor(loc$Type_habitat_broad, 
                                 levels=c("Natural","Artificial"))
  
  if(normality_p_value > 0.05 & homogeneity_habitat_p_value> 0.05) {
    
    
    test<-  t.test(loc$percentage_NIS_all~loc$Type_habitat_broad, alternative="greater")
    p_val<-test$p.value
    stat<-test$statistic
    result.data<-data.frame(p_val,stat,i)
    rownames(result.data)<-i
    result.data.fin.Nat<-rbind(result.data,result.data.fin.Nat)
  } else {
    
    test<-  wilcox.test(loc$percentage_NIS_all~loc$Type_habitat_broad, alternative="greater")
    p_val<-test$p.value
    stat<-test$statistic
    result.data<-data.frame(p_val,stat,i)
    rownames(result.data)<-i
    result.data.fin.Nat<-rbind(result.data,result.data.fin.Nat)
  }
}
res.prop$Greater<-"Artificial"
result.data.fin.Nat$Greater<-"Natural"
result.data.finAll$Greater<-"Artificial"

result.data.finAll<-rbind(res.prop,result.data.finAll,result.data.fin.Nat)


#### h.2.d Test if the proportion of NIS varies between localities ------

meta.with.proportions

shapiro.test(meta.with.proportions$percentage_NIS_all)

kruskal.test(meta.with.proportions$percentage_NIS_all~meta.with.proportions$Locality)

theme_set(theme_minimal())
boxplot.proportions<-meta.with.proportions %>%
  ggplot(aes(Locality,percentage_NIS_all))+
  geom_boxplot(color="grey8")+ 
  fig +
  theme(axis.text.x =  element_text(size=15,angle=90),
        axis.title.x = element_blank())+
  # ggtitle("Percentage of NIS across localities","Kruskal-Wallis test, P value = 0.07")+
  # facet_grid(~Type_habitat_broad, scale="free_x")+
  ylab ("Percentage NIS")

boxplot.proportions
# ggsave(plot=boxplot.proportions, "NIS/NIS_proportion_across_localities_no_habitat_distinction.pdf",
#        width=8,height=6,dpi=300)

# Not significant difference in prop of NIS across localities, neither for artificial or natural sites

## i. Greater richness of NIS in ports than outside -----

### i.1. NIS alpha diversity is greater inside than outside + boxplot ----

richNIS<-estimate_richness(NIS, measures="Observed")

richNIS$sample_name<-rownames(richNIS)

metaNIS<-sample_data(NIS) %>% data.frame()

meta.rich.NIS<-left_join(metaNIS,richNIS,by="sample_name")

# Test greater richness of NIS inside than outside 
shapiro.test(meta.rich.NIS$Observed)

wilcox.test(Observed~Type_habitat_broad, data=meta.rich.NIS, alternative="greater") 


# Test and generate boxplot throughout localities splittniig inside and outside

## Only inside 
rich.NIS.in<-meta.rich.NIS[meta.rich.NIS$Type_habitat_broad %in% "Artificial",]
rich.NIS.in %>% group_by(Locality) %>% summarise(mean_NIS=mean(Observed),
                                                 sd_NIS=sd(Observed))

shapiro.test(rich.NIS.in$Observed)

kruskal.test(Observed~Locality, data=rich.NIS.in) # 0.03621

tbl.rni<-as.data.frame(dunn.test::dunn.test(x=rich.NIS.in$Observed,g=rich.NIS.in$Locality, method = "bh"))
tbl.rni$Habitat<-"Artificial"

## Only outside

rich.NIS.out<-meta.rich.NIS[meta.rich.NIS$Type_habitat_broad %in% "Natural",]

rich.NIS.out %>% group_by(Locality) %>% summarise(mean_NIS=mean(Observed),
                                                           sd_NIS=sd(Observed))


shapiro.test(rich.NIS.out$Observed)

kruskal.test(Observed~Locality, data=rich.NIS.out) # 0.02347

tbl.rno<-as.data.frame(dunn.test::dunn.test(x=rich.NIS.out$Observed,g=rich.NIS.out$Locality, method = "bh"))
tbl.rno$Habitat<-"Natural"


NIS.rich.plot<-ggplot()+
  geom_boxplot(aes(rich.NIS.in$Locality,rich.NIS.in$Observed), fill="blue", alpha=0.7)+
  geom_boxplot(aes(rich.NIS.out$Locality,rich.NIS.out$Observed), fill="yellow3", alpha=0.7)+
  theme_minimal()+
  xlab("Locality")+ylab("Richness of NIS")+
  scale_x_discrete(labels=c("Ajaccio","Bastia","Bonifacio",
                            "Calvi","Nice","Port-la-Nouvelle",
                            "Porto-Vecchio", "Saint-Tropez",
                            "Sete_1","Sete_2","Seyne_sur_mer",
                            "Toulon"))+
  theme(axis.title.x = element_text(size = 15))+
  theme(axis.title.y = element_text(size = 15))+
  theme(axis.text.x = element_text(size = 14, angle=45,vjust = 0.5))+
  theme(axis.text.y = element_text(size = 14))

# Perform test on NIS across all localities> GREATER nis IN PORTS

nis.alpha.artif.greater<-alpha_tests.new(NIS, levels=c("Natural","Artificial"))  #lower NIS in natural habitats

print(nis.alpha.artif.greater)
nis.alpha.artif.greater$Greater<-"Greater NIS in Artificial"

# Greater Artificial - Note : I excluded Nice and Port la Nouvelle becasue 
# the former did not have any NIS outside, and for the latter the levene test was givin Nas as result blockng the loop. 
# I performed PLN separately

unique_localitiesNIS<-unique(sample_data(NIS)$Locality)
nis.alpha.artif.greater.locality<-data.frame()

for (locality in unique_localitiesNIS[c(1:4,7:11)]){
  
  ps.sub<-subset_samples(NIS, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub,levels=c("Natural","Artificial"))
  
  rownames(result.test)<-locality
  
  nis.alpha.artif.greater.locality<-rbind(nis.alpha.artif.greater.locality,result.test)
  
}


nis.alpha.artif.greater.locality$Greater<-"Greater NIS in Artificial"

print(nis.alpha.artif.greater.locality)

###### Peroform Port la Nouvelle : #

plv<-subset_samples(NIS, Locality == "Port-la-Nouvelle")

plv.rich<-estimate_richness(plv, measures="Observed")

plv.rich$sample_name<-rownames(plv.rich)

metaNIS.plv<-sample_data(plv) %>% data.frame()

metaNIS.plv<-left_join(metaNIS.plv,plv.rich,by="sample_name")

# Test greater richness of NIS inside than outside 
shapiro.test(metaNIS.plv$Observed)

metaNIS.plv$Type_habitat_broad<-factor(metaNIS.plv$Type_habitat_broad, levels = c("Artificial","Natural"))
t.pvl<-t.test(Observed~Type_habitat_broad, data=metaNIS.plv, alternative="greater") 


nis.alpha.res<-rbind(nis.alpha.artif.greater, nis.alpha.artif.greater.locality)
nis.alpha.res[12,1]<-t.pvl$p.value
nis.alpha.res[12,2]<-t.pvl$statistic
nis.alpha.res[12,3]<-"t.test"
nis.alpha.res[12,4]<-"Habitat"
nis.alpha.res[12,5]<-"Artificial"

rownames(nis.alpha.res)[11]<-"Port-la-Nouvelle"

# Richness plot all locality

richness.NIS.allLocs<-ggplot(meta.rich.NIS,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                       color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
                     label.y = 20)+   # show only *, **, etc.
  #hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  #facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

richness.NIS.allLocs


## l. Greater number of unique NIS found only inside port ----

### l.1. Venn diagram  -----

# Find all the NIS that are present both in the artificial and natural habitats and those that are shared between the two.

# Combine taxa table otu table and metadata

tax<-as.data.frame(tax_table(NIS))
tax$OTU<-rownames(tax)

otu<-as.data.frame(t(otu_table(NIS)))
otu$sample_name<-rownames(otu)


# Convert the binary tableOTU# Convert the binary table to a two-column dataset
long_otu <- otu %>%
  pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  select(sample_name, OTU)

metaNIS

metaNIS.Venn<-metaNIS %>% select(sample_name,Type_habitat_broad)

meta_long_otu<-left_join(metaNIS.Venn,long_otu,by="sample_name")

venn_data <- meta_long_otu %>%
  group_by(Type_habitat_broad) %>%
  summarise(OTUs = list(unique(OTU)))

# Prepare the OTU lists for the Venn diagram
venn_list <- list(
  Natural = venn_data$OTUs[venn_data$Type_habitat_broad == "Natural"][[1]],
  Artificial = venn_data$OTUs[venn_data$Type_habitat_broad == "Artificial"][[1]]
)


# Create and plot the Venn diagram with shared OTUs included
library(VennDiagram)
venn.plot <- venn.diagram(
  x = venn_list,
  category.names = c("Natural", "Artificial"),
  filename = NULL, # Don't save to a file
  output = TRUE,
  col = "black",
  fill = c("yellow2", "blue"),
  alpha = 0.7,
  cex = 2.5,
  fontface = "bold",
  fontfamily = "calibri",
  cat.cex = 2,
  cat.fontface = "plain",
  cat.fontfamily = "calibri"
)


grid.draw(venn.plot)


shared<-intersect(venn_list$Natural,venn_list$Artificial)

unique.artificial<-venn_list$Artificial[!venn_list$Artificial %in% venn_list$Natural]

unique.natural<-venn_list$Natural[!venn_list$Natural %in% venn_list$Artificial]

# combine with taxonomy info

shared<- tax %>% subset(OTU %in% shared )

unique.artificial<- tax %>% subset(OTU %in% unique.artificial )

unique.natural<- tax %>% subset(OTU %in% unique.natural )

# combine with taxonomy info

shared$Venn<-"Shared"

unique.artificial$Venn<-"Unique_Artificial"

unique.natural$Venn<-"Unique_Natural"

unique.shared.NIS<-rbind(unique.artificial,unique.natural,shared) %>% 
  select(FARTA,Venn, Phylum_Worms_Accepted,
         Class_Worms_Accepted,
         Order_Worms_Accepted,
         Family_Worms_Accepted,
         Genus_Worms_Accepted,
         Species_Worms_Accepted,
         Status_Med_Sea, 
         Status_NIS_table_collaps,
         Origin_NIS_table_collaps,
         Trait1,
         Trait2_Worms,
         Trait4_larvae)


# Plot for poster

theme_set(theme_minimal())

mobility<-ggplot(unique.shared.NIS
                 , aes(x=Trait1, fill = Trait1, color=Trait1))+
  geom_bar()+
  scale_fill_manual(values = c("orangered1","orange4"))+
  scale_color_manual(values = c("orangered1","orange4"))+
  ylab("Number of Species")+
  facet_grid(~Venn)+
  theme(strip.background = element_blank(), strip.text= element_text(size=12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme (axis.text.x = element_text(size=15))+
  theme (axis.text.y = element_text(size=15))+
  theme (axis.title.y= element_text(size=16))+
  theme(legend.text = element_text(size=12))+
  theme(legend.position = "none")+
  theme(title = element_text(size=17, vjust = 0.5))


depth<- ggplot(subset(unique.shared.NIS, !Trait2_Worms%in%""), 
               aes(x=Trait2_Worms, fill = Trait2_Worms, color=Trait2_Worms))+
  geom_bar()+
  scale_fill_manual(values = c("lightblue1","paleturquoise4"))+
  scale_color_manual(values = c("lightblue1","paleturquoise4"))+
  ylab("Number of Species")+
  facet_grid(~Venn)+
  theme(strip.background = element_blank(), strip.text= element_text(size=12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme (axis.text.x = element_text(size=15))+
  theme (axis.text.y = element_text(size=15))+
  theme (axis.title.y= element_text(size=16))+
  theme(legend.text = element_text(size=12))+
  theme(legend.position = "none")+
  theme(title = element_text(size=17, vjust = 0.5))

larvae<-ggplot(subset(unique.shared.NIS, !Trait4_larvae%in%"")
               , aes(x=Trait4_larvae, fill = Trait4_larvae, color=Trait4_larvae))+
  geom_bar()+
  scale_fill_manual(values = c("lightgoldenrod","olivedrab"))+
  scale_color_manual(values = c("lightgoldenrod","olivedrab"))+
  ylab("Number of Species")+
  facet_grid(~Venn)+
  theme(strip.background = element_blank(), strip.text= element_text(size=12)) +
  theme(panel.border = element_rect(colour = "black", fill = NA))+
  theme (axis.text.x = element_text(size=15))+
  theme (axis.text.y = element_text(size=15))+
  theme (axis.title.y= element_text(size=16))+
  theme(legend.text = element_text(size=12))+
  theme(legend.position = "none")+
  theme(title = element_text(size=17, vjust = 0.5))



mobility_depth<-ggarrange(mobility,
          depth,
          larvae, nrow=1,ncol=3)


## m. Beta diversity of NIS ---- 

#### m.1. Beta diversity test on all localities together ----

jaccard.NIS<-distance(NIS, method = "jaccard", binary=TRUE)
metaNIS<- sample_data(NIS) %>% data.frame() # extract the dataframe
metaNIS[,"Type_habitat_broad"]<-as.factor(metaNIS[,"Type_habitat_broad"])

# check beta dispersion - betadisperser test for homogeneity of variance

dispr.NIS<-betadisper(jaccard.NIS,metaNIS$Type_habitat_broad) 
permutest(dispr.NIS) # P value 0.18   
boxplot(dispr.NIS)

# Welch MANOVA - Robust distance-based multivariate analysis of variance (https://doi.org/10.1186/s40168-019-0659-9)

adonis2(jaccard.NIS~metaNIS$Type_habitat_broad) # 0.001 


#### m.1.1. Indicator species - NIS ####

set.seed(1234)
indicspecs.NIS<-run_indicspecies(NIS)
indicspecs.lines.nNIS<-capture.output(print(summary(indicspecs.NIS))) 
# writeLines(indicspecs.lines.noNIS,"Local/Indicspecs_full_output_ALL.txt")  


#### m.2. NIS dissimialrity components -----

unique_localities<-unique(metaNIS$Locality)
final.combine.NIS<-data.frame()

for (i in unique_localities) {
  
  loc<-subset_samples(NIS, Locality %in% i)
  meta.loc<-subset(metaNIS, Locality %in% i)
  
  betapair<- beta.pair(as.data.frame(t(otu_table(loc))), index.family="jaccard")
  beta.nest<-betapair$beta.jne
  beta.tu<-betapair$beta.jtu
  beta.full<-betapair$beta.jac
  
  # artif.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Natural"]
  # natural.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Artificial"]
  
  artif.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Artificial"]
  natural.samples<-meta.loc$sample_name[meta.loc$Type_habitat_broad=="Natural"]
  
  beta.full.sub<-as.matrix(beta.full) %>% data.frame()
  mean(beta.full)
  mean(beta.full.sub)
  beta.full.sub<- beta.full.sub[,c(natural.samples)]
  beta.full.sub<-beta.full.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in% artif.samples) 
  rownames(beta.full.sub)<-NULL
  beta.full.sub<- beta.full.sub %>% column_to_rownames(var="sample")
  
  
  beta.nest.sub<-as.matrix(beta.nest) %>% data.frame()
  mean(beta.nest)
  mean(beta.nest.sub)
  beta.nest.sub<- beta.nest.sub[,c(natural.samples)]
  beta.nest.sub<-beta.nest.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in% artif.samples) 
  rownames(beta.nest.sub)<-NULL
  beta.nest.sub<- beta.nest.sub %>% column_to_rownames(var="sample")
  
  
  beta.tu.sub<-as.matrix(beta.tu) %>% data.frame()
  mean(beta.tu)
  mean(beta.tu.sub)
  beta.tu.sub<- beta.tu.sub[,c(natural.samples)]
  beta.tu.sub<-beta.tu.sub %>% rownames_to_column(var="sample") %>%
    subset(sample %in% artif.samples) 
  rownames(beta.tu.sub)<-NULL
  beta.tu.sub<- beta.tu.sub %>% column_to_rownames(var="sample")
  
  beta.tu.sub$mean_turnover<- rowMeans(beta.tu.sub)
  
  beta.nest.sub$mean_nest<- rowMeans(beta.nest.sub)
  
  beta.full.sub$mean_full<- rowMeans(beta.full.sub)
  
  
  combine<-as.data.frame(cbind(beta.tu.sub$mean_turnover,beta.nest.sub$mean_nest,beta.full.sub$mean_full))
  colnames(combine)<-c("mean_turnover","mean_nestedness","mean_jaccard")
  rownames(combine)<-rownames(beta.tu.sub)
  
  combine$sample<-rownames(combine)
  
  combine<- pivot_longer(combine, -sample, names_to = "Group")
  combine$Locality<-i
  
  final.combine.NIS<-rbind(combine,final.combine.NIS)
  
}


theme_set(theme_minimal())

NIS.plot.components<-ggplot(final.combine.NIS , aes(Group, value, fill=Group))+
  geom_boxplot()+
  ylab("Distance")+xlab("")+
  ggtitle("Beta diversity components",
          "Artificial vs Natural sites")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))


## NIS in ports are more homogeneous than those outside? ----
dispr.NIS.dist<-as.data.frame(dispr.NIS$distances)

colnames(dispr.NIS.dist)[1]<-"Jaccard_distance"

dispr.NIS.dist$sample_name<-rownames(dispr.NIS.dist)

metaNIS.sub<-metaNIS %>% select (Type_habitat_broad, sample_name, Region) 

dispr.NIS.dist<-left_join(dispr.NIS.dist, metaNIS.sub, by="sample_name")

rownames(dispr.NIS.dist)<-NULL

dispr.NIS.dist<- dispr.NIS.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.NIS.dist$Jaccard_distance)

dispr.NIS.dist$Type_habitat_broad<-  factor(dispr.NIS.dist$Type_habitat_broad,
                                             levels = c("Artificial", "Natural"))

wilcox.test(Jaccard_distance~Type_habitat_broad, data=dispr.NIS.dist, paired=FALSE, alternative=c("less")) #0.16

# Plot

NIS_homogen.plot<-ggplot(dispr.NIS.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                                               fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_jitter(size=3)+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold3"))+
  ylab("Jaccard distance")+
  xlab("")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))



## n.2. NIS Heatmap ----


# Step 1: Detemine the number of replicates per site in which the NIS is found

#obtain the samples names in each locality

site.list<-metaNIS%>% group_by(Locality) %>% summarise(sites=as.vector(sample_name))

replicates.final<-data.frame(matrix(NA, nrow = 43, ncol = 0))

for (locality in unique(site.list$Locality)) {
  
  Locality<-subset(site.list, Locality %in% locality)
  
  sub<-as.data.frame(tax.countNIS.all[,Locality$sites])
  
  replicates<-as.data.frame(rowSums(sub))
  
  colnames(replicates)<-locality
  
  replicates.final<-cbind(replicates,replicates.final)
  
}

replicates.final$OTU_seed<-rownames(replicates.final)

tax.sub<-tax[,c("OTU_seed","FARTA","AphiaID","Phylum_Worms_Accepted",
                "collapsed_Markers",
                "Accepted_names_NIS_table",
                "Origin_NIS_table",
                "Status_NIS_table",
                "Recorded_in_North_East_Atlantic_Ocean_NIS_table",
                "Recorded_in_Baltic_Sea_NIS_table",
                "Recorded_in_Black_Sea_NIS_table", 
                "Status_Med_Sea",
                "Trait1",
                "Trait4_larvae")]

replicates.final<-left_join(replicates.final,tax.sub,by ="OTU_seed")

rownames(replicates.final)<-replicates.final$OTU_seed


# add column with "shared","unique artificial" and "unique natural" form the Venn diagram of NIS

replicates.final<-replicates.final %>% mutate(Habitat_Venn=case_when(OTU_seed%in%shared$OTU_seed~"Shared",
                                                                     OTU_seed%in%unique.artificial$OTU_seed~"Inside",
                                                                     OTU_seed%in%unique.natural$OTU_seed~"Outside"))

# write.csv2(replicates.final,file="NIS/All_NIS_replicates_localities.csv")

# Make heatmap with species and localities and numb of replicates

# load libraries
library(pheatmap)
library(grid)

# load matrix

mat.to.heatmap<-read_xlsx("NIS/matrix_NIS_med_sea_to_heatmap.xlsx") %>%column_to_rownames(var="FARTA")

# Order all the species by Phylum

mat.to.heatmap<- mat.to.heatmap %>% arrange(Phylum_Worms_Accepted)
# italiacize the species names 

species_names_italics<-lapply(rownames(mat.to.heatmap), function(x)bquote(italic(.(x))))

# Define colors 
heat_color<-colorRampPalette(c("seashell1","firebrick"))(100)


# Prepare annotation for row (species)
annotation_row <- data.frame(Habitat=mat.to.heatmap$Habitat_Venn,
                             Mobility=mat.to.heatmap$Trait1,
                             Status_Med_Sea = mat.to.heatmap$Status_Med_Sea,
                             Phylum= mat.to.heatmap$Phylum_Worms_Accepted)

rownames(annotation_row) <- rownames(mat.to.heatmap)

colnames(mat.to.heatmap)
mat.to.heatmap<-mat.to.heatmap[,-c(13:16)]

desired_order <- c("Port-la-Nouvelle", "Sete_1", "Sete_2", "Seyne-sur-mer", "Toulon", "Saint-Tropez",
                   "Nice", "Calvi", "Ajaccio", "Bonifacio", "Porto-Vecchio", "Bastia")

# Reorder columns of mat.to.heatmap
mat.to.heatmap <- mat.to.heatmap[, desired_order]

# Generate unique colors for Phylum (bright)
phylum_levels <- unique(annotation_row$Phylum_Worms_Accepted)
phylum_colors <- setNames(rainbow(length(phylum_levels)), phylum_levels)

# Black and white for Status
status_levels <- unique(annotation_row$Status_Med_Sea)
status_colors <- setNames(c("black", "white")[1:length(status_levels)], status_levels)

# Custom colors for Mobility
mobility_colors <- c("Mobile" = "deeppink", "Sessile" = "greenyellow")

# Custom colors for Habitat
# Blend yellow2 and blue for shared
shared_color <- colorRampPalette(c("yellow", "blue"))(3)[2]
habitat_colors <- c("Inside" = "blue", "Outside" = "yellow", "Shared" = shared_color)

# Define annotation colors list
annotation_colors <- list(
  Phylum_Worms_Accepted = phylum_colors,
  Status_Med_Sea = status_colors,
  Mobility = mobility_colors,
  Habitat = habitat_colors
)




# Create the heatmap
heatmap<-pheatmap(
  mat = mat.to.heatmap,
  color = heat_color,
  scale = "none",  # Disable scaling
  cluster_rows = FALSE,  # Keep species in original order
  cluster_cols = FALSE,  # Keep sites in original order
  annotation_row = annotation_row,  # Add status annotation
  annotation_colors = annotation_colors,
  legend = TRUE,  # Show the legend
  display_numbers = FALSE,  # Disable showing numbers on heatmap
  annotation_legend = TRUE,  # Show annotation legend
  fontsize_row = 10 ,
  fontsize_col = 12,
  labels_row = as.expression(species_names_italics))

heatmap


## o. Characterization of NIS by traits ----

NIS.tax<-as.data.frame(tax_table(NIS))


NIS.tax %>% group_by(Trait1)%>% summarise(freq=(n()*100)/42)
NIS.tax %>% group_by(Trait1_3rd_lev_fish)%>% summarise(freq=(n()*100)/42)


NIS.tax.trait2<-subset(NIS.tax, !Class_Worms_Accepted %in% c("Teleostei",
                                                             "Elasmobranchii"))

NIS.tax.trait2 %>% group_by(Trait2_Worms)%>% summarise(freq=(n()*100)/40)


NIS.tax %>% group_by(Trait4_larvae)%>% summarise(freq=(n()*100)/42)


### p. Traits do not account for different compositin of NIS inside and outside ports ----

### p.1. dbRDA ----

jaccard.matrix.NIS<-phyloseq::distance(NIS, method = "jaccard", binary=TRUE)

community.NIS<-as.data.frame(otu_table(NIS)) 


NIS.dbrda<-capscale(jaccard.matrix.NIS~Type_habitat_broad +
                      # percentage_NIS_all.sessile + 
                      Condition(Locality), 
                    distance="jaccard", data=metaNIS, 
                    comm = as.data.frame(t(community.NIS)))

anova.cca(NIS.dbrda)
anova.cca(NIS.dbrda, by ="terms",permutations = 999)
anova.cca(NIS.dbrda, by ="margin", permutations = 999)
RsquareAdj(NIS.dbrda)

sum<-summary(NIS.dbrda)
sum$cont$importance

# Test if the CAP1 is influenced by the traits 
cap.NIS<-as.data.frame(scores(NIS.dbrda, display = "species"))
cap.NIS$OTU_seed<-rownames(cap.NIS)


# Combine with the trait taxa table

trait.NIS<-as.data.frame(tax_table(NIS)) %>% select(FARTA,
                                                    OTU_seed,
                                                    Trait2_Worms,
                                                    Trait1,
                                                    Trait4_larvae,
                                                    Phylum_Worms_Accepted,
                                                    Class_Worms_Accepted)

cap.NIS<-left_join(cap.NIS,trait.NIS, by="OTU_seed") %>% column_to_rownames(var="OTU_seed")


# cap.NIS$Trait1_3rd_lev_fish<-as.factor(cap.NIS$Trait1_3rd_lev_fish)
cap.NIS$Trait1<-as.factor(cap.NIS$Trait1)
cap.NIS$Trait2_Worms<-as.factor(cap.NIS$Trait2_Worms)
cap.NIS$Phylum_Worms_Accepted<-as.factor(cap.NIS$Phylum_Worms_Accepted)

# Simple test of variables on CAP1 

# noNIS
cap.NIS2<-subset(cap.NIS, !Trait2_Worms%in%"")
cap.NIS3<-subset(cap.NIS, !Trait1%in%"")
cap.NIS4<-subset(cap.NIS, !Trait4_larvae%in%"")

# shapiro.test(cap.NIS1$CAP1)
shapiro.test(cap.NIS4$CAP1)
shapiro.test(cap.NIS2$CAP1)
shapiro.test(cap.NIS3$CAP1)
shapiro.test(cap.NIS$CAP1)


wilcox.test(cap.NIS3$CAP1~cap.NIS3$Trait1)
wilcox.test(cap.NIS2$CAP1~cap.NIS2$Trait2_Worms)
# wilcox.test(cap.NIS2$CAP1~cap.NIS2$Trait4_larvae) # Signif

kruskal.test(cap.NIS$CAP1~cap.NIS$Phylum_Worms_Accepted)
kruskal.test(cap.NIS$CAP1~cap.NIS$Class_Worms_Accepted)

# Plot site only dbRDA for NIS 

# Plot sites and taxa together and colo by mobility
sites.NIS<-as.data.frame(scores(NIS.dbrda, display = "sites"))
sites.NIS$sample_name<-rownames(sites.NIS)
sites.NIS<-left_join(sites.NIS, df[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.NIS<-ggplot() +
  geom_point(data=sites.NIS, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.NIS,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                     fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
  theme_minimal()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))


## NIS PCoA


NIS_PCoA<-ordinate(NIS, method = "PCoA", distance="jaccard") 

sites_PCoA.NIS<-as.data.frame(NIS_PCoA$vectors)
sites_PCoA.NIS$sample_name<-rownames(sites_PCoA.NIS)
NIS_PCoA$vectors
sites_PCoA.NIS<-left_join(metaNIS,sites_PCoA.NIS, by="sample_name")


# Plot PCoA 

theme_set(theme_minimal())

NIS_pcoa.plot<-plot_ordination(NIS, NIS_PCoA, color= "Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
  scale_color_manual(values = c("darkblue","yellow3"))+ 
  labs(title = paste0(" PCoA - COMM"))+
  ggConvexHull:: geom_convexhull(data=sites_PCoA.NIS,aes(x = Axis.1, y = Axis.2, color=Type_habitat_broad,
                                                     fill=Type_habitat_broad),alpha=0.0)+
fig


# Dispersio plot:
centroid.NIS<-sites_PCoA.NIS%>% group_by(Type_habitat_broad)%>%summarise(centroid_Axis.1=mean(Axis.1),centroid_Axis.2=mean(Axis.2))
sites_PCoA.NIS <-  left_join(sites_PCoA.NIS, centroid.NIS, by="Type_habitat_broad")

NIS_pcoa.plot.centroids<-ggplot(data=sites_PCoA.NIS, aes(Axis.1, Axis.2, color=Type_habitat_broad))+
  geom_point(size=4, alpha=1, aes(color=Type_habitat_broad))+
  geom_point(data=centroid.NIS, aes(x=centroid_Axis.1,y=centroid_Axis.2, color=Type_habitat_broad))+
  geom_segment(aes(xend = centroid_Axis.1, yend = centroid_Axis.2),
               alpha = 0.5)+
  scale_color_manual(values = c("darkblue","yellow3"))+ 
  labs(title = paste0(" PCoA - NIS"))
                                                          fill=Type_habitat_broad),alpha=0.0)

