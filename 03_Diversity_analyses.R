# Diversity analyses (Alpha and Beta) and Funcitonal traits analyses as reported in the article "An ecological perspective on the portuarization syndrome: 
# eDNA reveals species-rich communities, non-indigenous species hotspots, and biotic homogenization in ports"

# The phyloseq object called in this script has been obtained by the combination 
# of the OTU and TAXA generated for three different markers 12S,18S, 
# COI accoridng to the "Preparation table script" also present in the github repository

# The analyses performed are organized as follow:

# 1) General information about dataset and communities
# 2) section a. - Alpha diversity analyses on whole community (COMM) dataset, and on native species (NATS) and Uncertain group datasets 
# 3) section b., c., e., f. - Beta diversity analyses on whole community (COMM) dataset, and on native species (NATS) and Uncertain group datasets 
# including analyses of beta diversity components (turnover and nestedness), biotic homogenization, and correlation between species scores of 1st dbRDA axis and funcional and taxonommical traits
# 4) section d. Venn diagrams for unique and shared species in ports and in natural habitats
# 5) section g-p: same diversity analyses of secitons a-d only on NIS community



# Author: Ginevra Lilli 

# Load libraries -----------------------------------------------------------------
pkgs <- c("phyloseq","ggplot2","dplyr","vegan","tidyr",
          "cowplot","gridExtra","broom","worrms","car","purrr",
        "tibble","ggpubr","MicEco","readxl","indicspecies",
        "adespatial","stringr","pheatmap","grid","ggforce","betapart","glmmTMB")


lapply(pkgs, require,character.only=TRUE)

# color palettes ---------------------------------------------------------------
c25<-c("dodgerblue2",   "#E31A1C" ,      "green4",        "#6A3D9A",      
"#FF7F00"     ,  "black"  ,       "gold1" ,        "skyblue2",     
"#FB9A99"      , "palegreen2",    "#CAB2D6",       "#FDBF6F",      
"gray70"        ,"khaki2",        "maroon",        "orchid1",      
"deeppink1"     ,"blue1",         "steelblue4",    "darkturquoise",
"green1"        ,"yellow4",       "yellow3",       "darkorange4",  
"brown" )

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

# combined_taxa_ps_summer.sp<-readRDS("ps_538_species_summer2021.RData")

combined_taxa_ps_summer.sp<-readRDS("~ps_536_species.RData")

# Load the Species occurrence table 
# 
# traits.ninis.df<-readxl::read_excel("C:/Users/Ginevra Lilli/Dropbox/Postdoc_CEAB/MarEEE-Med/Writing/Article/Manuscrip_MarEEEMed/MS_Reviewed/Final_to_Ecography/Final/reviewed/Ecography/Species_occurrence_table.xlsx")

# Generate the three datasets:
# Whole community (COMM=NIS + NATIVE + UNCERTAIN), 
# NATS (Native community accoridng to litterature)
# NIS (includin NIS already recorded in the Mediterraenan Sea NIS-MED, first record in Mediterraenan Sea but already recorded in Europe NIS-EU, and Cryptogenic)
# Uncertain (species that could represented putative new NIS, but also errors in the classification )

COMM <- combined_taxa_ps_summer.sp

metaCOMM<-sample_data(COMM)%>% data.frame()

NAT <-subset_taxa(combined_taxa_ps_summer.sp, STATUS %in% c("Native"))

NAT<-prune_samples(sample_sums(NAT)!=0,NAT) 

metaNAT<-sample_data(NAT)%>% data.frame()

NIS<-subset_taxa(combined_taxa_ps_summer.sp, STATUS %in% c("NIS-MED","NIS-EU","Cryptogenic"))

NIS<-prune_samples(sample_sums(NIS)!=0,NIS) 

metaNIS<-sample_data(NIS)%>% data.frame()

UNCERT<-subset_taxa(combined_taxa_ps_summer.sp, STATUS %in% c("Uncertain",""))

UNCERT<-prune_samples(sample_sums(UNCERT)!=0,UNCERT) 

metaUNCERT<-sample_data(UNCERT)%>% data.frame()


# for plots 
unique_phylum<-unique(as.data.frame(tax_table(COMM))$Phylum_Chordata_splitted)
phylum_colors <- setNames(c25[1:length(unique_phylum)], unique_phylum)

# 0. General info about dataset and community ----

# Determine how many species were classified with 1, 2 or 3 markers and which taxa 


COMM.tax<-as.data.frame(tax_table(COMM))

freq.collaps.markers<- COMM.tax %>% group_by(collapsed_Markers, Phylum_Chordata_splitted) %>% summarise(freq=n())

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
  "18S"=155,
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
desired_order <- c(  "COI","18S","12S", "12S ,COI","18S ,COI","12S ,18S", "12S ,18S ,COI")

# Apply that order to the factor
freq.collaps.markers$collapsed_Markers <- factor(freq.collaps.markers$collapsed_Markers, 
                                                 levels = desired_order)


markers.plot.by.taxa<-ggplot(freq.collaps.markers,aes(x=collapsed_Markers, y=perc.by.taxa, fill=Phylum_Chordata_splitted))+ 
  geom_bar(stat="identity", position = "stack")+
  # ggtitle("Number of Species classified by combining the 3 markers","Classificatio at Phylum level")+
  scale_fill_manual(values = phylum_colors)+fig+
  theme(legend.position = "right")




# Add total phya abundances bar in the bar plot 
total.phy<-COMM.tax %>% group_by(Phylum_Chordata_splitted) %>% summarise(n=n())
total.phy$percnt<-(total.phy$n/sum(total.phy$n))*100
total.phy$collapsed_Markers<-"All data"

colnames(freq.collaps.markers)
colnames(total.phy)[3]<-"perc.by.taxa"
total.phy<-total.phy[,-2]
total.phy$tot.spec.by.mark<-NA
total.phy$freq<-NA

freq.collaps.markers<-rbind(freq.collaps.markers,total.phy)


# Save figures:
ggsave(plot=markers.plot.by.taxa,"Final_files_pubblication/Figures/markers.plot.by.taxa.pdf",
       dpi=300,height = 8,width=9)

# Average number of species per locality----- 

sample.per.loc<-metaCOMM %>% group_by(Locality)%>%summarise(samples=sample_name)
otu<-as.data.frame(otu_table(COMM))

otu.pooled<-as.data.frame(matrix(nrow=536,ncol=12))
colnames(otu.pooled)[1:12]<-unique(metaCOMM$Locality)
rownames(otu.pooled)<-taxa_names(COMM)

loc<-unique(sample.per.loc$Locality)

for (i in loc){
  sub<-subset (sample.per.loc,Locality %in% i)
  try<- otu %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled[,i]<-try[,i]
}

otu.pooled[otu.pooled>1]<-1 # transform in presence absence 
mean(colSums(otu.pooled))# 125.25
sd(colSums(otu.pooled))#22.47069


# How many species by each locality and habitat?
species.count.COMM<-otu.pooled
species.count.COMM[species.count.COMM>1]<-1
species.count.COMM<-colSums(species.count.COMM)

# Barplot taxa by localities pooled ----

theme_set(theme_minimal())

# barplot by habitat and locality


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
                               "Porto-Vecchio", "Bastia", "All_localities"),
                     labels=c("Port-la-Nouvelle", "Sete_1", "Sete_2", "Seyne-sur-mer",
                              "Toulon",
                              "Saint-Tropez", "Nice", "Calvi", "Ajaccio", "Bonifacio", 
                              "Porto-Vecchio", "Bastia", "All_localities"))+
    guides(fill = guide_legend(title = "Phylum"))+
    fig+
    theme(axis.text.x = element_text(size=13, angle=90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.text = element_text(size=16))
  return(plot.bar)
}




### NAT - barplot habitat x loclaity ----

df<-sample_data(NAT)%>%data.frame()

sample.per.loc.NAT<-df %>% 
  group_by(Locality, Type_habitat_broad)%>%summarise(samples=sample_name)

sample.per.loc.NAT$combined<-paste0(sample.per.loc.NAT$Locality,"_",
                                    sample.per.loc.NAT$Type_habitat_broad)

otu.NAT<-as.data.frame(otu_table(NAT))

otu.pooled.NAT<-as.data.frame(matrix(nrow=354,ncol=24))
colnames(otu.pooled.NAT)[1:12]<-paste0(unique(df$Locality),"_Natural")
colnames(otu.pooled.NAT)[13:24]<-paste0(unique(df$Locality),"_Artificial")
rownames(otu.pooled.NAT)<-taxa_names(NAT)

loc.hab.NAT<-unique(sample.per.loc.NAT$combined)

for (i in loc.hab.NAT){
  sub<-subset (sample.per.loc.NAT,combined %in% i)
  try<- otu.NAT %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled.NAT[,i]<-try[,i]
}

# Add column where you report the sum of species in all ports and all nats

otu.pooled.NAT$All_Artificial<-rowSums(otu.pooled.NAT[,grep("Artificial",names(otu.pooled.NAT))])
otu.pooled.NAT$All_Natural<-rowSums(otu.pooled.NAT[,grep("Natural",names(otu.pooled.NAT))])

df$combined.loc.hab<-paste0(df$Locality,"_",df$Type_habitat_broad)

df.sub<-subset(df, !duplicated(combined.loc.hab ))
df.sub[nrow(df.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Artificial","All_localities","Artificial")
df.sub[nrow(df.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Natural","All_localities", "Natural")

rownames(df.sub)<-df.sub$combined.loc.hab

# How many species by each locality and habitat?
species.count<-otu.pooled.NAT
species.count[species.count>1]<-1
species.count<-colSums(species.count)

species.count.NAT<-as.data.frame(species.count)
colnames(species.count.NAT)<-"Native_species_count"


# remake ps.subset
NAT.sub<-phyloseq(otu_table(as.matrix(otu.pooled.NAT),taxa_are_rows = TRUE),
                    tax_table(tax_table(NAT)),taxa_names(taxa_names(NAT)),
                    sample_data(df.sub))


barplot.NAT.pooled<-barplot_habitat.loc(NAT.sub, "Phylum_Chordata_splitted")


### NIS - barplot x locality x habitat -----

sample.per.loc.NIS<-metaNIS %>% group_by(Locality, Type_habitat_broad)%>%summarise(samples=sample_name)
sample.per.loc.NIS$combined<-paste0(sample.per.loc.NIS$Locality,"_",sample.per.loc.NIS$Type_habitat_broad)
otu.NIS<-as.data.frame(otu_table(NIS))


otu.pooled.NIS<-as.data.frame(matrix(nrow=50,ncol=24))
colnames(otu.pooled.NIS)[1:12]<-paste0(unique(metaNIS$Locality),"_Natural")
colnames(otu.pooled.NIS)[13:24]<-paste0(unique(metaNIS$Locality),"_Artificial")
rownames(otu.pooled.NIS)<-taxa_names(NIS)

loc.hab.NIS<-unique(sample.per.loc.NIS$combined)

for (i in loc.hab.NIS){
  sub<-subset (sample.per.loc.NIS,combined %in% i)
  try<- otu.NIS %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled.NIS[,i]<-try[,i]
}

otu.pooled.NIS$All_Artificial<-rowSums(otu.pooled.NIS[,grep("Artificial",names(otu.pooled.NIS))])
otu.pooled.NIS$All_Natural<-rowSums(otu.pooled.NIS[,grep("Natural",names(otu.pooled.NIS))])

metaNIS$combined.loc.hab<-paste0(metaNIS$Locality,"_",metaNIS$Type_habitat_broad)

metaNIS.sub<-subset(metaNIS, !duplicated(combined.loc.hab ))
metaNIS.sub[nrow(metaNIS.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Artificial","All_localities","Artificial")
metaNIS.sub[nrow(metaNIS.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Natural","All_localities", "Natural")

rownames(metaNIS.sub)<-metaNIS.sub$combined.loc.hab

# How many species by each locality and habitat?
species.count.NIS<-otu.pooled.NIS
species.count.NIS[species.count.NIS>1]<-1
species.count.NIS<-colSums(species.count.NIS)

species.count.NIS<-as.data.frame(species.count.NIS)
colnames(species.count.NIS)<-"NIS_species_count"

# remake ps.subset
NIS.sub<-phyloseq(otu_table(as.matrix(otu.pooled.NIS),taxa_are_rows = TRUE),
                    tax_table(tax_table(NIS)),taxa_names(taxa_names(NIS)),
                    sample_data(metaNIS.sub))

barplot.NIS.pooled<-barplot_habitat.loc(NIS.sub, "Phylum_Chordata_splitted")


barplot_pooled<-ggarrange(barplot.NAT.pooled,
                                  barplot.NIS.pooled, nrow=1, ncol=2, common.legend = TRUE)


ggsave(plot=barplot_pooled, "Final_files_pubblication/Figures/barplot_NATIVE_NIS.pdf", dpi=300,
       width = 15,height=8)

ggsave(plot=barplot_pooled, "Final_files_pubblication/Figures/barplot_NATIVE_NIS.jpeg", dpi=300,
       width = 15,height=8)


### UNCERTAIN- barplot x habitat x locality ----

sample.per.loc.UNC<-metaUNCERT %>% group_by(Locality, Type_habitat_broad)%>%summarise(samples=sample_name)
sample.per.loc.UNC$combined<-paste0(sample.per.loc.UNC$Locality,"_",sample.per.loc.UNC$Type_habitat_broad)
otu.UNC<-as.data.frame(otu_table(UNCERT))


otu.pooled.UNC<-as.data.frame(matrix(nrow=132,ncol=24))
colnames(otu.pooled.UNC)[1:12]<-paste0(unique(metaNIS$Locality),"_Natural")
colnames(otu.pooled.UNC)[13:24]<-paste0(unique(metaNIS$Locality),"_Artificial")
rownames(otu.pooled.UNC)<-taxa_names(UNCERT)

loc.hab<-unique(sample.per.loc.UNC$combined)

for (i in loc.hab){
  sub<-subset (sample.per.loc.UNC,combined %in% i)
  try<- otu.UNC %>% select (sub$samples)
  try[,i]<-rowSums(try)
  otu.pooled.UNC[,i]<-try[,i]
}


otu.pooled.UNC$All_Artificial<-rowSums(otu.pooled.UNC[,grep("Artificial",names(otu.pooled.UNC))])
otu.pooled.UNC$All_Natural<-rowSums(otu.pooled.UNC[,grep("Natural",names(otu.pooled.UNC))])

metaUNCERT$combined.loc.hab<-paste0(metaUNCERT$Locality,"_",metaUNCERT$Type_habitat_broad)

metaUNCERT.sub<-subset(metaUNCERT, !duplicated(combined.loc.hab ))
metaUNCERT.sub[nrow(metaUNCERT.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Artificial","All_localities","Artificial")
metaUNCERT.sub[nrow(metaUNCERT.sub)+1,c("combined.loc.hab", "Locality","Type_habitat_broad")]<-c("All_Natural","All_localities", "Natural")

rownames(metaUNCERT.sub)<-metaUNCERT.sub$combined.loc.hab


# How many species by each locality and habitat?
species.count.UNC<-otu.pooled.UNC
species.count.UNC[species.count.UNC>1]<-1
species.count.UNC<-colSums(species.count.UNC)

species.count.UNC<-as.data.frame(species.count.UNC)
colnames(species.count.UNC)<-"Uncertain_species_count"

# remake ps.subset
UNC.sub<-phyloseq(otu_table(as.matrix(otu.pooled.UNC),taxa_are_rows = TRUE),
                  tax_table(tax_table(UNCERT)),taxa_names(taxa_names(UNCERT)),
                  sample_data(metaUNCERT.sub))

barplot.UNC.pooled<-barplot_habitat.loc(UNC.sub, "Phylum_Chordata_splitted")

ggsave(plot=barplot.UNC.pooled, "Final_files_pubblication/Figures/barplot_UNCERTAIN.pdf",
       width = 10, height=8, dpi=300)


## Save number of species tables 
species.count.NAT$ID<-rownames(species.count.NAT)
species.count.NIS$ID<-rownames(species.count.NIS)
species.count.UNC$ID<-rownames(species.count.UNC)

species.count<-left_join(species.count.NAT,species.count.NIS,by="ID")
species.count<-left_join(species.count,species.count.UNC,by="ID")

write.csv(species.count, file="Final_files_pubblication/Tables/species.counts.csv")

## a. Lower richness of species in ports than outside ----
## Alpha diversity analyses performed on whole community (COMM) and RESS species (RESS) one after the other. Analyses performed on whole are flagged with COMM (C) and those by RESS species as RESS (R)

## COMM ----

### C.a.1. Alpha diversity (richness) all localities together on full comunity (COMM)
### C.a.2. Alpha diversity by locality on full comunity (COMM)
### C.a.3. Alpha diversity by phylum and locality on full comunity (COMM)

## One tail tests ----

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
    
    test.result.habitat<-t.test(Observed~Type_habitat_broad, data=meta.rich) # anova for 3 sites
    
    pvalue<-test.result.habitat$p.value
    stat<-test.result.habitat$statistic
    test<-"t.test"
    variable.tested<-"Habitat"
    
    result.habitat<-data.frame(pvalue,stat,test,variable.tested)
    
    
  } else {
    
    test.result.habitat<-wilcox.test(Observed~Type_habitat_broad, data=meta.rich, paired = FALSE, exact=FALSE) # kruskall-wallis
    
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
C.alpha.test.sp$test_type<-"two_tails"


# First perform one tail test 
C.alpha.test.sp.ot<-alpha_tests.new(COMM, levels = c("Artificial","Natural"))  
rownames(C.alpha.test.sp.ot)<-"All locations"
print(C.alpha.test.sp.ot)
C.alpha.test.sp.ot$test_type<-"one_tail_lower_Port"

# combine:
C.alpha.test.sp<-rbind(C.alpha.test.sp,C.alpha.test.sp.ot)

##### C.a.2. Alpha diversity by locality ----

unique_localities<-unique(as.data.frame(sample_data(COMM))$Locality)

# Different richness inside and outside ports - two tail test

C.result.locality<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(COMM, Locality == locality)
  
  result.test<-alpha_tests.two.tail(ps.sub)
  
  rownames(result.test)<-locality
  
  C.result.locality<-rbind(C.result.locality,result.test)
  
}


print(C.result.locality)

C.result.locality$BH.adjusted.p.values<-p.adjust(C.result.locality$pvalue, method = "BH")
C.result.locality$Holm.adjusted.p.values<-p.adjust(C.result.locality$pvalue, method = "holm")

C.result.locality$test_type<-"two_tails"

# save

C.result.alpha.two.tail<-rbind(C.alpha.test.sp,
                               C.result.locality)


# Lower richness inside ports - one tail test

C.result.locality.ot<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(COMM, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels = c("Artificial", "Natural"))
  
  rownames(result.test)<-locality
  
  C.result.locality.ot<-rbind(C.result.locality.ot,result.test)
  
}

C.result.locality.ot$test_type<-"one_tail_lower_Port"


C.result.locality.ot$BH.adjusted.p.values<-p.adjust(C.result.locality.ot$pvalue, method = "BH")
C.result.locality.ot$Holm.adjusted.p.values<-p.adjust(C.result.locality.ot$pvalue, method = "holm")




# Lower richness inside ports - one tail test

C.result.locality.lower.outside<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(COMM, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels = c("Natural", "Artificial"))
  
  rownames(result.test)<-locality
  
  C.result.locality.lower.outside<-rbind(C.result.locality.lower.outside,result.test)
  
}

C.result.locality.lower.outside$test_type<-"one_tail_lower_Natural"

C.result.locality.lower.outside$BH.adjusted.p.values<-p.adjust(C.result.locality.lower.outside$pvalue, method = "BH")
C.result.locality.lower.outside$Holm.adjusted.p.values<-p.adjust(C.result.locality.lower.outside$pvalue, method = "holm")



# combine result tables

combined_community_div<-rbind(C.alpha.test.sp,C.result.locality,C.result.locality.ot)

#### C.a.3a. Alpha diversity by phylum  ----
unique_phylum<-unique(as.data.frame(tax_table(COMM))$Phylum_Chordata_splitted)

C.alpha.test.phyla<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(COMM, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.two.tail(ps.sub.phyl)
  
  result.test$Phlyum<-phylum
  
  C.alpha.test.phyla<-rbind(C.alpha.test.phyla,result.test)
}

C.alpha.test.phyla
C.alpha.test.phyla$test_type<-"two_tails"

# One tail

C.alpha.test.phyla.ot<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(COMM, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.new(ps.sub.phyl, levels = c("Artificial","Natural"))
  
  result.test$Phlyum<-phylum
  
  C.alpha.test.phyla.ot<-rbind(C.alpha.test.phyla.ot,result.test)
}

C.alpha.test.phyla.ot
C.alpha.test.phyla.ot$test_type<-"one_tail_lower_Port"

combined_community_div$Phlyum<-"All_phyla"
combined_community_div<-rbind(combined_community_div, C.alpha.test.phyla, C.alpha.test.phyla.ot)

# Save table
write.csv(combined_community_div,file="Final_files_pubblication/Tables/Alpha_diversity_tests_outputs_COMM.csv")


#### C.a.4.a Plot richness across the whole dateset - COMM ----

COMM.meta.rich

COMM.richness.allLocs<-ggplot(COMM.meta.rich,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                               color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
  #                    label.y = 100)+   # show only *, **, etc.
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

COMM.richness.allLocs

#### C.a.4.b Plot richness across localities and phyla - COMM ----

plot.alpha.multiple<-function(ps.sub){
  
  COMM_rich<-estimate_richness(ps.sub, measures="Observed")
  COMM_rich$sample_name<-rownames(COMM_rich)
  COMM.meta.rich<-left_join(metaCOMM,COMM_rich,by="sample_name")
  
  ggplot(COMM.meta.rich,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                            color=Type_habitat_broad))+
    geom_boxplot()+
    geom_point()+
    xlab("")+ylab("")+
    # scale_y_continuous(limits = c(15, 100)) +
    # labs(title = "Richness of Species by habitat type")+
    # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
    #                    label.y = 100)+   # show only *, **, etc.
    #hide.ns = TRUE) +     # hide 'ns' when not significant
    ggtitle(as.data.frame(tax_table(ps.sub))$Phylum_Chordata_splitted)+
    scale_fill_manual(values=c("darkblue","yellow2"))+
    scale_color_manual(values=c("dodgerblue","gold1"))+
    # stat_compare_means(method = "t.test",
    #                    label = "p.signif")+   # show only *, **, etc.
    # 
    # facet_grid(~Locality, scale="free_x")+
    # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
    # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
    fig+
    theme(axis.text.x = element_blank(),
          legend.position = "top",
          legend.text=element_text(size=9),
          legend.title=element_text(size=7),
          axis.title.x = element_text(size=7),
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



plot.rich.loc<-ggarrange(plot.rich.loc$Ajaccio+ggtitle("Ajaccio") ,
                         plot.rich.loc$Bastia+ggtitle("Bastia") ,
                         plot.rich.loc$Bonifaccio+ggtitle("Bonifacio") ,
                         plot.rich.loc$Calvi+ggtitle("Calvi") ,
                         plot.rich.loc$Nice+ggtitle("Nice") ,
                         plot.rich.loc$`Port-la-Nouvelle`+ggtitle("Port-la-Nouvelle") ,
                         plot.rich.loc$`Seyne-sur-mer`+ggtitle("Seyne-sur-mer") ,
                         plot.rich.loc$Sete_1+ggtitle("Sete_1") ,
                         plot.rich.loc$Sete_2+ggtitle("Sete_2") ,
                         plot.rich.loc$`Saint-Tropez`+ggtitle("Saint-Tropez") ,
                         plot.rich.loc$`Porto-Vecchio`+ggtitle("Porto-Vecchio") ,
                         plot.rich.loc$Toulon+ggtitle("Toulon") ,
                         nrow=1,ncol=12, common.legend = TRUE)


ggsave(plot=plot.rich.loc,"Final_files_pubblication/Figures/alpha_boxplot_by_locality.pdf",
       dpi=300, width=24, height=6)


# Phyla 
theme_set(theme_minimal())

plot.rich.phyla<-list()

for (phylum in unique_phylum){
  
  ps.sub<-subset_taxa(COMM, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  plot<-plot.alpha.multiple(ps.sub)
  
  plot.rich.phyla[phylum]<-list(plot)
}

# Keep only those phyla with significantly different richness
plot.rich.phyla.comb<-ggarrange(plot.rich.phyla$Arthropoda,#
                                plot.rich.phyla$Cnidaria,#
                                plot.rich.phyla$Vertebrata,#
                                plot.rich.phyla$Bryozoa,#
                                plot.rich.phyla$Annelida,#
                                plot.rich.phyla$Tunicata,#
                                plot.rich.phyla$Porifera,#
                                plot.rich.phyla$Nemertea,#
                                plot.rich.phyla$Rotifera,#
                                plot.rich.phyla$Mollusca,
                                nrow=1,ncol=10, common.legend = TRUE)

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
ggsave(plot=plot.rich.phyla.comb, "Final_files_pubblication/Figures/alpha_boxplot_alpha_phyla.pdf",
       width = 20, height = 6, dpi=300)


##### GENERALIZED LINEAR MIXED MODELS - ALTERNATIVE TO WILCOX TO TAKE INTO ACCOUNT NESTED STRUCTURE OF DATA ----

# Perform generalized linear mixed models for the whole dataset and repeated for each phyla
shapiro.test(COMM.meta.rich$Observed)

# Set the levels to have as reference the Artificial - Port one
COMM.meta.rich$Type_habitat_broad<-  factor(COMM.meta.rich$Type_habitat_broad,
                                            levels = c("Artificial", "Natural"))


mod<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
             family=nbinom2(),
             data=COMM.meta.rich)

# Test overdispersion 

res<-DHARMa::simulateResiduals(mod)
disp.test<-DHARMa::testDispersion(res)
disp.test$p.value

# mod2<-glmmTMB(Observed~Type_habitat_broad*Locality, 
#                            family=nbinom2,
#                            data=COMM.meta.rich) # this overfit so better to not use it even if the AIC is th elowest. Too little replicates per each habotat to rely on it.

# mod2<-glmmTMB(Observed~Type_habitat_broad*Locality, 
#              family=nbinom2,
#              data=COMM.meta.rich) # this overfit so better to not use it even if the AIC is th elowest. Too little replicates per each habotat to rely on it.
# 
# 
# mod3<-glmmTMB(Observed~Type_habitat_broad+(Type_habitat_broad|Locality), 
#               family=nbinom2,
#               data=COMM.meta.rich)


summary.mod.COMM<-summary(mod)
confint(mod)
# Wald intervals (fast)

summary.mod.COMM<-summary(mod2)
confint(mod)

# summarmod# summary(mod2)
# summary(mod3)

## get he coefficients

summary.mod.COMM<-as.data.frame(summary.mod.COMM$coefficients$cond)[2,]
summary.mod.COMM$Phylum<-"All"
summary.mod.COMM$ref.habitat<-"Port"

# By each phylum 
unique_phylum<-unique(as.data.frame(tax_table(COMM))$Phylum_Chordata_splitted)

glmm.phyl<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(COMM, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  rich<-estimate_richness(ps.sub.phyl, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps.sub.phyl) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  mod.phyl<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                       family=poisson,
                       data=meta.rich)
  
  
  
  result.test<-as.data.frame(summary(mod.phyl)$coefficients$cond)
  result.test<-result.test[2,]
  result.test$Phylum<-phylum
  
  dispersion<-DHARMa::simulateResiduals(mod.phyl)
  disp.test<-DHARMa::testDispersion(dispersion)
  result.test$overdisp.test.pvalue<-disp.test$p.value
  
    
  glmm.phyl<-rbind(glmm.phyl,result.test)
}

glmm.phyl

glmm.phyl$ref.habitat<-"Port"

glmm.COMM<-rbind(summary.mod.COMM,glmm.phyl)

write.csv(glmm.COMM,file="Final_files_pubblication/Tables/glmmTMB_COMM.csv")


## Confidence intervals: 

conf.glmm.phyl<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(COMM, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  rich<-estimate_richness(ps.sub.phyl, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps.sub.phyl) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  result.test<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                       family=nbinom2,
                       data=meta.rich)
  
  
  
  confint<-as.data.frame(confint(result.test))
  confint<-confint[2,]
  confint$Phylum<-phylum
  
  
  conf.glmm.phyl<-rbind(conf.glmm.phyl,confint)
}

View(conf.glmm.phyl)
## NAT (N) ----

### a.1. Alpha diversity (richness) all localities together on NATIVE species 
### a.2. Alpha diversity by locality on NATIVE species 
### a.3. Alpha diversity by phylum and locality on NATIVE species 

unique_localities<-unique(sample_data(NAT)$Locality)

# Boxplot of diversity

# by habitat

rich.NAT<-estimate_richness(NAT, measures="Observed")
rich.NAT$sample_name<-rownames(rich.NAT)
meta.rich.NAT<-left_join(metaNAT,rich.NAT,by="sample_name")

#### N.a.1. Alpha diversity (richness) all localities together ----

# Test different richness in ports and in natural habitats
# First perform two tail test 
N.alpha.test.sp<-alpha_tests.two.tail(NAT)  
rownames(N.alpha.test.sp)<-"All locations"
print(N.alpha.test.sp)
N.alpha.test.sp$test_type<-"two_tails"

# Then perform one tail test - lower richness in ports
N.alpha.test.sp.ot<-alpha_tests.new(NAT, levels = c("Artificial","Natural"))  
rownames(N.alpha.test.sp.ot)<-"All locations"
print(N.alpha.test.sp.ot)
N.alpha.test.sp.ot$test_type<-"one_tail_lower_Port"

# combine:
N.alpha.test.sp<-rbind(N.alpha.test.sp,N.alpha.test.sp.ot)

#### N.a.2. Alpha diversity by locality ----

# Different richness in Artificial habitats
result.locality.NAT<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(NAT, Locality == locality)
  
  result.test<-alpha_tests.two.tail(ps.sub)
  rownames(result.test)<-locality
  
  result.locality.NAT<-rbind(result.locality.NAT,result.test)
  
}

print(result.locality.NAT)
result.locality.NAT$test_type<-"two_tails"
result.locality.NAT$adjusted.pvalu.holm<-p.adjust(result.locality.NAT$pvalue,method="holm")
result.locality.NAT$adjusted.pvalu.BH<-p.adjust(result.locality.NAT$pvalue,method="BH")

# Lower richness in ports habitats 

result.locality.NAT.lower.Artificial<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(NAT, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Artificial","Natural"))
  rownames(result.test)<-locality
  
  result.locality.NAT.lower.Artificial<-rbind(result.locality.NAT.lower.Artificial,result.test)
  
}

print(result.locality.NAT.lower.Artificial)
result.locality.NAT.lower.Artificial$test_type<-"one_tail_lower_Port"
result.locality.NAT.lower.Artificial$adjusted.pvalu.holm<-p.adjust(result.locality.NAT.lower.Artificial$pvalue,method="holm")
result.locality.NAT.lower.Artificial$adjusted.pvalu.BH<-p.adjust(result.locality.NAT.lower.Artificial$pvalue,method="BH")

combined_NATIVE_div<-rbind(N.alpha.test.sp,result.locality.NAT,result.locality.NAT.lower.Artificial)
combined_NATIVE_div$Phlyum<-"All_Phyla"



# Lower richness in natural habitats 

result.locality.NAT.lower.Natural<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(NAT, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Natural","Artificial"))
  rownames(result.test)<-locality
  
  result.locality.NAT.lower.Natural<-rbind(result.locality.NAT.lower.Natural,result.test)
  
}

print(result.locality.NAT.lower.Natural)
result.locality.NAT.lower.Natural$test_type<-"one_tail_lower_Natural"

result.locality.NAT.lower.Natural$adjusted.pvalu.holm<-p.adjust(result.locality.NAT.lower.Natural$pvalue,method="holm")
result.locality.NAT.lower.Natural$adjusted.pvalu.BH<-p.adjust(result.locality.NAT.lower.Natural$pvalue,method="BH")



combined_NATIVE_div<-rbind(N.alpha.test.sp,result.locality.NAT,result.locality.NAT.lower.Artificial)
combined_NATIVE_div$Phlyum<-"All_Phyla"




#### N.a.3. Alpha diversity by phylum  ----
unique_phylum<-unique(as.data.frame(tax_table(NAT))$Phylum_Chordata_splitted)

N.alpha.test.phyla<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(NAT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.two.tail(ps.sub.phyl)
  
  result.test$Phlyum<-phylum
  
  N.alpha.test.phyla<-rbind(N.alpha.test.phyla,result.test)
}

N.alpha.test.phyla
N.alpha.test.phyla$test_type<-"two_tails"

# One tail

N.alpha.test.phyla.ot<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(NAT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.new(ps.sub.phyl, levels = c("Artificial","Natural"))
  
  result.test$Phlyum<-phylum
  
  N.alpha.test.phyla.ot<-rbind(N.alpha.test.phyla.ot,result.test)
}

N.alpha.test.phyla.ot
N.alpha.test.phyla.ot$test_type<-"one_tail_lower_Port"

combined_NATIVE_div<-rbind(combined_NATIVE_div, N.alpha.test.phyla, N.alpha.test.phyla.ot)

# Save table
write.csv(combined_NATIVE_div,file="Final_files_pubblication/Tables/Alpha_diversity_tests_outputs_NATIVE.csv")


#### N.a.4. Plot richness across all dataset ---- 

NAT.richness.allLocs<-ggplot(meta.rich.NAT,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                       color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
  #                    label.y = 100)+   # show only *, **, etc.
  # #hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  #facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

NAT.richness.allLocs

#### GENERALIZED MIXED MODELS FOR NATIVE COMMUNITY ----

meta.rich.NAT

shapiro.test(meta.rich.NAT$Observed)

# Perform generalized linear mixed models for the whole dataset and repeated for each phyla

# Set the levels to have as reference the Artificial - Port one
meta.rich.NAT$Type_habitat_broad<-  factor(meta.rich.NAT$Type_habitat_broad,
                                            levels = c("Artificial", "Natural"))


mod.NAT<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
             family=nbinom2,
             data=meta.rich.NAT)


# mod2<-glmmTMB(Observed~Type_habitat_broad*Locality, 
#              family=nbinom2,
#              data=COMM.meta.rich) # this overfit so better to not use it even if the AIC is th elowest. Too little replicates per each habotat to rely on it.
# 
# 
# mod3<-glmmTMB(Observed~Type_habitat_broad+(Type_habitat_broad|Locality), 
#               family=nbinom2,
#               data=COMM.meta.rich)


summary.mod.NAT<-summary(mod.NAT)
confint(mod.NAT)

res<-DHARMa::simulateResiduals(mod.NAT)
DHARMa::testDispersion(res)

## get he coefficients

summary.mod.NAT<-as.data.frame(summary.mod.NAT$coefficients$cond)[2,]
summary.mod.NAT$Phylum<-"All"
summary.mod.NAT$ref.habitat<-"Port"

# By each phylum 
unique_phylum<-unique(as.data.frame(tax_table(NAT))$Phylum_Chordata_splitted)

glmm.phyl.NAT<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(NAT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  rich<-estimate_richness(ps.sub.phyl, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps.sub.phyl) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  result.test<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                       family=nbinom2,
                       data=meta.rich)
  
  
  
  result.test<-as.data.frame(summary(result.test)$coefficients$cond)
  result.test<-result.test[2,]
  result.test$Phylum<-phylum
  
  
  glmm.phyl.NAT<-rbind(glmm.phyl.NAT,result.test)
}

glmm.phyl.NAT

glmm.phyl.NAT$ref.habitat<-"Port"

glmm.NAT<-rbind(summary.mod.NAT,glmm.phyl.NAT)

write.csv(glmm.NAT,file="Final_files_pubblication/Tables/glmmTMB_NAT.csv")



## UNCERTAIN (U) ----

### a.1. Alpha diversity (richness) all localities together on UNCERTAIN species 
### a.2. Alpha diversity by locality on UNCERTAIN species 
### a.3. Alpha diversity by phylum and locality on UNCERTAIN species 

# Boxplot of diversity
unique_localities<-unique(sample_data(UNCERT)$Locality)

# by habitat

rich.UNC<-estimate_richness(UNCERT, measures="Observed")
rich.UNC$sample_name<-rownames(rich.UNC)
meta.rich.UNC<-left_join(metaNAT,rich.UNC,by="sample_name")

#### U.a.1. Alpha diversity (richness) all localities together ----

# Test different richness in ports and in natural habitats
# First perform two tail test 
U.alpha.test.sp<-alpha_tests.two.tail(NAT)  
rownames(U.alpha.test.sp)<-"All locations"
print(U.alpha.test.sp)
U.alpha.test.sp$test_type<-"two_tails"


# Then perform one tail test - lower richness in ports
U.alpha.test.sp.ot<-alpha_tests.new(UNCERT, levels = c("Artificial","Natural"))  
rownames(U.alpha.test.sp.ot)<-"All locations"
print(U.alpha.test.sp.ot)
U.alpha.test.sp.ot$test_type<-"one_tail_lower_Port"

# combine:
U.alpha.test.sp<-rbind(U.alpha.test.sp,U.alpha.test.sp.ot)

#### U.a.2. Alpha diversity by locality ----

# Different richness in Artificial habitats
result.locality.UNC<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(UNCERT, Locality == locality)
  
  result.test<-alpha_tests.two.tail(ps.sub)
  rownames(result.test)<-locality
  
  result.locality.UNC<-rbind(result.locality.UNC,result.test)
  
}

print(result.locality.UNC)
result.locality.UNC$test_type<-"two_tails"



result.locality.UNC$adjust.pval.holm<-p.adjust(result.locality.UNC$pvalue,method="holm")
result.locality.UNC$adjust.pval.BH<-p.adjust(result.locality.UNC$pvalue,method="BH")



# Lower richness in port habitats 

result.locality.UNC.lower.port<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(UNCERT, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Artificial","Natural"))
  rownames(result.test)<-locality
  
  result.locality.UNC.lower.port<-rbind(result.locality.UNC.lower.port,result.test)
  
}

print(result.locality.UNC.lower.port)
result.locality.UNC.lower.port$test_type<-"one_tail_lower_port"

result.locality.UNC.lower.port$adjust.pval.holm<-p.adjust(result.locality.UNC.lower.port$pvalue,method="holm")
result.locality.UNC.lower.port$adjust.pval.BH<-p.adjust(result.locality.UNC.lower.port$pvalue,method="BH")



# Lower richness in natural habitats 

result.locality.UNC.lower.Natural<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(UNCERT, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Natural","Artificial"))
  rownames(result.test)<-locality
  
  result.locality.UNC.lower.Natural<-rbind(result.locality.UNC.lower.Natural,result.test)
  
}

print(result.locality.UNC.lower.Natural)
result.locality.UNC.lower.Natural$test_type<-"one_tail_lower_Natural"

result.locality.UNC.lower.Natural$adjust.pval.holm<-p.adjust(result.locality.UNC.lower.Natural$pvalue,method="holm")
result.locality.UNC.lower.Natural$adjust.pval.BH<-p.adjust(result.locality.UNC.lower.Natural$pvalue,method="BH")



# Lower richness in natural habitats 

result.locality.UNC.lower.Natural<-data.frame()

for (locality in unique_localities){
  
  ps.sub<-subset_samples(UNCERT, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Natural","Artificial"))
  rownames(result.test)<-locality
  
  result.locality.UNC.lower.Natural<-rbind(result.locality.UNC.lower.Natural,result.test)
  
}

print(result.locality.UNC.lower.Natural)
result.locality.UNC.lower.Natural$test_type<-"one_tail_lower_Natural"

result.locality.UNC.lower.Natural$adjust.pval.holm<-p.adjust(result.locality.UNC.lower.Natural$pvalue,method="holm")
result.locality.UNC.lower.Natural$adjust.pval.BH<-p.adjust(result.locality.UNC.lower.Natural$pvalue,method="BH")



combined_UNCERT_div<-rbind(U.alpha.test.sp,result.locality.UNC,result.locality.UNC.lower.Artificial)
combined_UNCERT_div$Phlyum<-"All_Phyla"




result.locality.UNC$adjust.pval.holm<-p.adjust(result.locality.UNC$pvalue,method="holm")
result.locality.UNC$adjust.pval.BH<-p.adjust(result.locality.UNC$pvalue,method="BH")

#### U.a.3. Alpha diversity by phylum  ----
unique_phylum<-unique(as.data.frame(tax_table(UNCERT))$Phylum_Chordata_splitted)

U.alpha.test.phyla<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(UNCERT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.two.tail(ps.sub.phyl)
  
  result.test$Phlyum<-phylum
  
  U.alpha.test.phyla<-rbind(U.alpha.test.phyla,result.test)
}

U.alpha.test.phyla
U.alpha.test.phyla$test_type<-"two_tails"

# One tail

U.alpha.test.phyla.ot<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(UNCERT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.new(ps.sub.phyl, levels = c("Artificial","Natural"))
  
  result.test$Phlyum<-phylum
  
  U.alpha.test.phyla.ot<-rbind(U.alpha.test.phyla.ot,result.test)
}

U.alpha.test.phyla.ot
U.alpha.test.phyla.ot$test_type<-"one_tail_lower_Port"

combined_UNCERT_div<-rbind(combined_UNCERT_div, U.alpha.test.phyla, U.alpha.test.phyla.ot)

# Save table
write.csv(combined_UNCERT_div,file="Final_files_pubblication/Tables/Alpha_diversity_tests_outputs_UNCERTAIN.csv")


#### U.a.4. Plot richness across all dataset ---- 

UNC.richness.allLocs<-ggplot(meta.rich.UNC,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                               color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
  #                    label.y = 100)+   # show only *, **, etc.
  # #hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  #facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

UNC.richness.allLocs

#### GENERALIZED LINEAR MIXE DMODEL - UNCETAIN -----


shapiro.test(meta.rich.UNC$Observed)
# Perform generalized linear mixed models for the whole dataset and repeated for each phyla

# Set the levels to have as reference the Artificial - Port one
meta.rich.UNC$Type_habitat_broad<-  factor(meta.rich.UNC$Type_habitat_broad,
                                           levels = c("Artificial", "Natural"))


mod.UNC<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                 family=nbinom2,
                 data=meta.rich.UNC)


# mod2<-glmmTMB(Observed~Type_habitat_broad*Locality, 
#              family=nbinom2,
#              data=COMM.meta.rich) # this overfit so better to not use it even if the AIC is th elowest. Too little replicates per each habotat to rely on it.
# 
# 
# mod3<-glmmTMB(Observed~Type_habitat_broad+(Type_habitat_broad|Locality), 
#               family=nbinom2,
#               data=COMM.meta.rich)


summary.mod.UNC<-summary(mod.UNC)
confint(mod.UNC)

res<-DHARMa::simulateResiduals(mod.UNC)
DHARMa::testDispersion(res)

## get he coefficients

summary.mod.UNC<-as.data.frame(summary.mod.UNC$coefficients$cond)[2,]
summary.mod.UNC$Phylum<-"All"
summary.mod.UNC$ref.habitat<-"Port"

# By each phylum 
unique_phylum<-unique(as.data.frame(tax_table(UNCERT))$Phylum_Chordata_splitted)

glmm.phyl.UNC<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(UNCERT, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  rich<-estimate_richness(ps.sub.phyl, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps.sub.phyl) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  result.test<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                       family=nbinom2,
                       data=meta.rich)
  
  
  
  result.test<-as.data.frame(summary(result.test)$coefficients$cond)
  result.test<-result.test[2,]
  result.test$Phylum<-phylum
  
  
  glmm.phyl.UNC<-rbind(glmm.phyl.UNC,result.test)
}

glmm.phyl.UNC

glmm.phyl.UNC$ref.habitat<-"Port"

glmm.UNC<-rbind(summary.mod.UNC,glmm.phyl.UNC)

write.csv(glmm.UNC,file="Final_files_pubblication/Tables/glmmTMB_UNC.csv")

## b. Communities in ports are dissimilar from those outside ----
## Analyses performed on COMM and NAT community 

### b.1. Unconstrained ordination - PCoA ----

#### COMM ----

COMM_PCoA<-ordinate(COMM, method = "PCoA", distance="jaccard") 

# Make color variable
loc.cols <- c(
  "Ajaccio"          = "#E69F00",
  "Bastia"           = "#56B4E9",
  "Porto-Vecchio"    = "#009E73",
  "Bonifaccio"       = "#F0E442",
  "Calvi"            = "#0072B2",
  "Nice"             = "#D55E00",
  "Port-la-Nouvelle" = "#CC79A7",
  "Seyne-sur-mer"    = "#999999",
  "Sete_1"           = "#117733",
  "Sete_2"           = "#332288",
  "Saint-Tropez"     = "#88CCEE",
  "Toulon"           = "#44AA99"
)

# Plot MDS 

theme_set(theme_minimal())

COMM_PCoA.plot<-plot_ordination(COMM, COMM_PCoA, color= "Locality", 
                                shape="Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Locality))+
  scale_color_manual(values = loc.cols)+ 
  labs(title = paste0(" PCoA - COMM"))+fig+
  theme(legend.position = "none")

COMM_PCoA.plot

ggsave(plot=COMM_PCoA.plot, "Final_files_pubblication/Figures/Unconstrained_MDS_loc_habitat.LEGEND.pdf",
       dpi=300,width=6,height = 6)

ggsave(plot=COMM_PCoA.plot, "Final_files_pubblication/Figures/Unconstrained_MDS_loc_habitat_COMM.pdf",
       dpi=300,width=6,height = 6)



#### NAT ----

NAT_PCoA<-ordinate(NAT, method = "PCoA", distance="jaccard") 


# Plot PCoA 

theme_set(theme_minimal())

NAT_PCoA.plot<-plot_ordination(NAT, NAT_PCoA, color= "Locality", 
                               shape="Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Locality))+
  scale_color_manual(values = loc.cols)+ 
  labs(title = paste0(" PCoA - NAT"))+fig+
  theme(legend.position = "none")

NAT_PCoA.plot

ggsave(plot=NAT_PCoA.plot, "Final_files_pubblication/Figures/Unconstrained_MDS_loc_habitat_NAT.pdf",
       dpi=300,width=6,height = 6)


#### UNCERTAIN ----

UNC_PCoA<-ordinate(UNCERT, method = "PCoA", distance="jaccard") 

# Plot PCoA 

theme_set(theme_minimal())

UNC_PCoA.plot<-plot_ordination(UNCERT, UNC_PCoA, color= "Locality", 
                               shape="Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Locality))+
  scale_color_manual(values = loc.cols)+ 
  labs(title = paste0(" PCoA - UNCERTAIN"))+fig+
  theme(legend.position = "none")

UNC_PCoA.plot

ggsave(plot=UNC_PCoA.plot, "Final_files_pubblication/Figures/Unconstrained_MDS_loc_habitat_UNCERT.pdf",
       dpi=300,width=6,height = 6)




### b.2. Beta diversity - tests ----

#### COMM - C.b.2. dbRDA ----

jaccard.COMM<-distance(COMM, method = "jaccard", binary=TRUE)
metaCOMM[,"Type_habitat_broad"]<-as.factor(metaCOMM[,"Type_habitat_broad"])

# check beta dispersion - betadisperser test for homogeneity of variance
fst.rare<-as.dist(list_D.1s$meanD_tots)
adonis2(fst.rare~Province, data=meta[!meta$Sample_Name %in% ""])

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


# Plot simple dbRDA 

# Reporte the importance of the axes 
summ.comm<-summary(COMM.dbrda)
summ.comm$cont$importance[2, "CAP1"] # importanc eof axes
summ.comm$cont$importance[2, "MDS1"] # importanc eof axes

# Plot dbRDA simple by habitats
# Plot site only dbRDA for NIS 
# Reporte the importance of the axes 


# Plot sites and taxa together and colo by mobility
sites.COMM<-as.data.frame(scores(COMM.dbrda, display = "sites"))
sites.COMM$sample_name<-rownames(sites.COMM)
sites.COMM<-left_join(sites.COMM, df[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.COMM<-ggplot() +
  geom_point(data=sites.COMM, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.COMM,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                     fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
  xlab(paste0("CAP 1 (", round(summ.comm$cont$importance[2,"CAP1"]*100, 1),"%)"))+
  ylab(paste0("PCoA 1 (", round(summ.comm$cont$importance[2,"MDS1"]*100, 1),"%)"))+
  theme_minimal()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))



ggsave(plot=dbrda.COMM, "Final_files_pubblication/Figures/dbRDA_COMM_simple.pdf",
       dpi=300, width=7,height=6)


#### C.b.3. Indicator species ----
run_indicspecies<-function(ps){
  
  df<- sample_data(ps) %>% data.frame() # extract the dataframe
  
  otu<-as.data.frame(otu_table(ps))
  
  otu$ids<-rownames(otu)
  
  tax<-as.data.frame(tax_table(ps))
  
  tax$ids<-rownames(tax)
  
  tax<-tax[,c("Species","ids")]
  
  otu<-left_join(otu,tax,by="ids") %>% column_to_rownames(var="Species")
  
  otu<-select(otu,-ids)
  
  otu<-t(otu) # indicspecies requires sites as rows and species as columns
  
  groups<-df$Type_habitat_broad
  
  indicspecs<-multipatt(otu, groups, 
                        control = how(nperm=999))
  return(indicspecs)
}

set.seed(1234)
indicspecs.COMM<-run_indicspecies(COMM)
indicspecs.lines.COMM<-capture.output(print(summary(indicspecs.COMM))) 


# Classification of indicator species 

indicspecs.COMM.sign<-indicspecs.COMM$sign %>% subset(p.value<=0.05)
indicspecs.COMM.sign$Species<-rownames(indicspecs.COMM.sign)
COMM.tax.indicator<- as.data.frame(tax_table(COMM))%>% subset(Species %in% rownames(indicspecs.COMM.sign))
COMM.tax.indicator<-left_join(indicspecs.COMM.sign, COMM.tax.indicator, by="Species")
colnames(COMM.tax.indicator)
COMM.tax.indicator<-COMM.tax.indicator[,c(1:17,25)]

write.csv(COMM.tax.indicator, file="Final_files_pubblication/Tables/Indicator_species_full_list.csv")


### N.b.1. Beta diversity test - NAT ----

jaccard.NAT<-distance(NAT, method = "jaccard", binary=TRUE)
 
metaNAT

# check beta dispersion - betadisperser test for homogeneity of variance
dispr.NAT<-betadisper(jaccard.NAT,metaNAT$Type_habitat_broad) 
permutest(dispr.NAT) 


#### NAT - C.b.2 dbRDDA -----

jaccard.NAT

community<-as.data.frame(otu_table(NAT)) 

NAT.dbrda<-capscale(jaccard.NAT~Type_habitat_broad +
                      Condition(Locality), # set locality as condition since it not the focus of the analysis 
                    distance="jaccard", data=metaNAT, 
                    comm = as.data.frame(t(community)))


anova.cca(NAT.dbrda)
anova.cca(NAT.dbrda, by ="terms",permutations = 999)
RsquareAdj(NAT.dbrda)


### U.b.1. Beta diversity - UNCERTAIN ----
# Test
jaccard.UNC<-distance(UNCERT, method = "jaccard", binary=TRUE)

metaUNCERT

# check beta dispersion - betadisperser test for homogeneity of variance
dispr.UNC<-betadisper(jaccard.UNC,metaUNCERT$Type_habitat_broad) 
permutest(dispr.UNC) 

#### UNC - C.b.2 - dbRDA ----
# Perform dbRDA 

jaccard.UNC

community.UNC<-as.data.frame(otu_table(UNCERT)) 

UNCERT.dbrda<-capscale(jaccard.UNC~Type_habitat_broad +
                         Condition(Locality), # set locality as condition since it not the focus of the analysis 
                       distance="jaccard", data=metaUNCERT, 
                       comm = as.data.frame(t(community.UNC)))


anova.cca(UNCERT.dbrda)
anova.cca(UNCERT.dbrda, by ="terms",permutations = 999)
RsquareAdj(UNCERT.dbrda)

# c. Biotic homogenization in ports is greater than in natural habitats ----
## Analyses repeated on COMM, NATIVE and UNCERTAIN species here. 


## COMM ----

dispr.COMM.dist<-as.data.frame(dispr.COMM$distances)

colnames(dispr.COMM.dist)[1]<-"Jaccard_distance"

dispr.COMM.dist$sample_name<-rownames(dispr.COMM.dist)

metaCOMM.sub<-metaCOMM %>% select (Type_habitat_broad, sample_name, Locality) 

dispr.COMM.dist<-left_join(dispr.COMM.dist, metaCOMM.sub, by="sample_name")

rownames(dispr.COMM.dist)<-NULL

dispr.COMM.dist<- dispr.COMM.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.COMM.dist$Jaccard_distance)

dispr.COMM.dist$Type_habitat_broad<-  factor(dispr.COMM.dist$Type_habitat_broad,
                                             levels = c("Artificial", "Natural"))


# is dispersion lower n ports ? Is homogenization greater in ports?
t.test(Jaccard_distance~Type_habitat_broad, data=dispr.COMM.dist, paired=FALSE, alternative=c("less")) # 0.001226


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


ggsave(plot=comm_homogen.plot, "Final_files_pubblication/Figures/biotic_homogenization.COMM.pdf",
       dpi=300,width=7, height=6)


## NAT  ----

dispr.NAT.dist<-as.data.frame(dispr.NAT$distances)

colnames(dispr.NAT.dist)[1]<-"Jaccard_distance"

dispr.NAT.dist$sample_name<-rownames(dispr.NAT.dist)

metaNAT.sub<-metaNAT %>% select (Type_habitat_broad, sample_name) 

dispr.NAT.dist<-left_join(dispr.NAT.dist, metaNAT.sub, by="sample_name")

rownames(dispr.NAT.dist)<-NULL

dispr.NAT.dist<- dispr.NAT.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.NAT.dist$Jaccard_distance)

dispr.NAT.dist$Type_habitat_broad<-  factor(dispr.NAT.dist$Type_habitat_broad,
                                             levels = c("Artificial", "Natural"))

t.test(Jaccard_distance~Type_habitat_broad, data=dispr.NAT.dist, paired=FALSE, alternative=c("less")) # 0.04493

# Plot

NAT.homogen.plot<-ggplot(dispr.NAT.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                        fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  ylab("Jaccard distance")+
  xlab("")+
  ggtitle("Native")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))




## UNCERTAIN  ----

dispr.UNC.dist<-as.data.frame(dispr.UNC$distances)

colnames(dispr.UNC.dist)[1]<-"Jaccard_distance"

dispr.UNC.dist$sample_name<-rownames(dispr.UNC.dist)

metaUNCERT.sub<-metaUNCERT %>% select (Type_habitat_broad, sample_name) 

dispr.UNC.dist<-left_join(dispr.UNC.dist, metaUNCERT.sub, by="sample_name")

rownames(dispr.UNC.dist)<-NULL

dispr.UNC.dist<- dispr.UNC.dist %>% column_to_rownames(var="sample_name")

# T test - Artificial sites are less dissimilar

shapiro.test(dispr.UNC.dist$Jaccard_distance)

dispr.UNC.dist$Type_habitat_broad<-  factor(dispr.UNC.dist$Type_habitat_broad,
                                            levels = c("Artificial", "Natural"))

t.test(Jaccard_distance~Type_habitat_broad, data=dispr.UNC.dist, paired=FALSE, alternative=c("less")) # 0.1619

# Plot

UNC.homogen.plot<-ggplot(dispr.UNC.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                                             fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  ylab("Jaccard distance")+
  xlab("")+
  ggtitle("Uncertain")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "none",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))



#d.  Venn diagram : unique and shared species inside and outside ----

# Find all the Native that are present both in the artificial and natural habitats and those that are shared between the two.

venn.plot<-function(ps){
  
  tax<-as.data.frame(tax_table(ps))
  tax$OTU<-rownames(tax)
  
  otu<-as.data.frame(t(otu_table(ps)))
  otu$sample_name<-rownames(otu)
  
  # Convert the binary tableOTU# Convert the binary table to a two-column dataset
  long_otu<- otu %>%
    pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
    filter(Presence == 1) %>%
    dplyr::select(sample_name, OTU)
  
  metaVenn<-sample_data(ps)%>% data.frame()
  
  metaVenn<-metaVenn %>% dplyr::select(sample_name,Type_habitat_broad)
  
  meta_long_otu<-left_join(metaVenn,long_otu,by="sample_name")
  
  venn_data <- meta_long_otu %>%
    group_by(Type_habitat_broad) %>%
    summarise(OTUs = list(unique(OTU)))
  
  # Prepare the OTU lists for the Venn diagram
  # venn_list.COMM <- list(
  #   Outside = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Natural"][[1]],
  #   Inside = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Artificial"][[1]]
  # )
  
  venn_list <- list(
    Natural = venn_data$OTUs[venn_data$Type_habitat_broad == "Natural"][[1]],
    Artificial = venn_data$OTUs[venn_data$Type_habitat_broad == "Artificial"][[1]]
  )
  
  # Create and plot the Venn diagram with shared OTUs included
  library(VennDiagram)
  venn.plot<- venn.diagram(
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
    # cat.fontface = "normal",
    cat.fontfamily = "calibri"
  )
  
  
  grid.draw(venn.plot)
  
  
}

# Venn COMM - NAT - UNCERT ----
venn.plot(COMM)
venn.plot(NAT)
venn.plot(UNCERT)

##### d.1 Make table with shared and unique species found in each habitat ----

tax.COMM<-as.data.frame(tax_table(COMM))
tax.COMM$OTU<-rownames(tax.COMM)

otu.COMM<-as.data.frame(t(otu_table(COMM)))
otu.COMM$sample_name<-rownames(otu.COMM)

# Convert the binary tableOTU# Convert the binary table to a two-column dataset
long_otu.COMM<- otu.COMM %>%
  pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  dplyr::select(sample_name, OTU)

metaVenn.COMM<-sample_data(COMM)%>% data.frame()

metaVenn.COMM<-COMM_metaVenn %>% dplyr::select(sample_name,Type_habitat_broad)

meta_long_otu.COMM<-left_join(metaVenn.COMM,long_otu.COMM,by="sample_name")

venn_data.COMM <- meta_long_otu.COMM %>%
  group_by(Type_habitat_broad) %>%
  summarise(OTUs = list(unique(OTU)))


venn_list.COMM <- list(
  Natural = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Natural"][[1]],
  Artificial = venn_data.COMM$OTUs[venn_data.COMM$Type_habitat_broad == "Artificial"][[1]]
)


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
  select(Species,Venn, Phylum_Chordata_splitted,
         Class,
         Order,
         Family,
         Genus,
         Species,
         STATUS)

write.csv(COMM_unique.shared, file="Final_files_pubblication/Tables/Species_unique_shared_different_habitats.csv")

### NAT ----

# Find all the Native that are present both in the artificial and natural habitats and those that are shared between the two.

# Combine taxa table otu table and metadata
tax.NAT<-as.data.frame(tax_table(NAT))
tax.NAT$OTU<-rownames(tax.NAT)

otu.NAT<-as.data.frame(t(otu_table(NAT)))
otu.NAT$sample_name<-rownames(otu.NAT)

# Convert the binary tableOTU# Convert the binary table to a two-column dataset
long_otu.NAT <- otu.NAT %>%
  pivot_longer(cols = -sample_name, names_to = "OTU", values_to = "Presence") %>%
  filter(Presence == 1) %>%
  dplyr::select(sample_name, OTU)

metaVenn.NAT<-sample_data(NAT)%>% data.frame()

metaVenn<-metaVenn %>% dplyr::select(sample_name,Type_habitat_broad)

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













#e. Beta diversity components and Traits tests -----
### e.1. COMM - Result of beta diversity components: test between nestedness and turnover ----

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

ggsave(plot=COMM_beta.components, "Final_files_pubblication/Figures/beta.components.boxplot.COMM.pdf",
       dpi=300,width=7, height=6)

## e.2. COMM - Traits explain the differences between habitats ----

## Dissimilarities between communities inside and outside ports are accounted for by different species traits

COMM_tax<-as.data.frame(tax_table(COMM))

COMM_tax %>% group_by(Trait1) %>% summarise(freq=(n()*100)/536) 
COMM_tax %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/536)

COMM_tax.trait2<-subset(COMM_tax, !Class %in% c("Teleostei",
                                                               "Elasmobranchii"))

COMM_tax.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/386)# all mobile species without fishes (benthic and zooplankton)

# Vertical distribution of fish
COMM_tax.fish<-subset(COMM_tax, Class %in% c("Teleostei",
                                                            "Elasmobranchii"))

COMM_tax.fish %>% group_by(Fishbase_depth_broad) %>% summarise(freq=(n()*100)/150)

# # Perform dbRDA 
# 
# jaccard.COMM
# 
# COMM_otu<-as.data.frame(otu_table(COMM)) 
# 
# COMM.dbrda<-capscale(jaccard.COMM~Type_habitat_broad +
#                        Condition(Locality), # set locality as condition since it not the focus of the analysis 
#                      distance="jaccard", data=metaCOMM, 
#                      comm = as.data.frame(t(COMM_otu)))
# 
# 
# anova.cca(COMM.dbrda)
# anova.cca(COMM.dbrda, by ="terms",permutations = 999)
# RsquareAdj(COMM.dbrda) 


# Test traits and plot dbRDA

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
                                                      Phylum_Chordata_splitted,
                                                      Class,
                                                      STATUS,
                                                      Species)


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
kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad)   
dunn.test::dunn.test(x=cap.COMM.3$CAP1, g=cap.COMM.3$Fishbase_depth_broad, method="bh")

COMM_traits.cap1.res[3,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[3,"Trait"]<-"Distribution in water column (Fishes)"
COMM_traits.cap1.res[3,"p value"]<-kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad)$p.value
COMM_traits.cap1.res[3,"stat"]<-kruskal.test(cap.COMM.3$CAP1~cap.COMM.3$Fishbase_depth_broad)$statistic
colnames(COMM_traits.cap1.res)[5]<-"significant_pairwise_comparison"
COMM_traits.cap1.res[3,"significant_pairwise_comparison"]<-"reef-associated vs bethopelagic (p=0.04),reef-associated vs demersal (p=0.05), slightly pelagic vs benthopelagic p=0.06"


#  PHylum
cap.COMM$Phylum_Chordata_splitted<-as.factor(cap.COMM$Phylum_Chordata_splitted)
kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Chordata_splitted) # 0.0009871

COMM_res.phyl<-FSA::dunnTest(x=cap.COMM$CAP1,g=cap.COMM$Phylum_Chordata_splitted, method = "bh")
COMM_res.phyl<-COMM_res.phyl$res
COMM_res.phyl<-COMM_res.phyl$Comparison[COMM_res.phyl$P.adj<=0.05] # "Arthropoda - Tunicata"    "Echinodermata - Tunicata" "Arthropoda - Vertebrata" 

COMM_traits.cap1.res[4,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[4,"Trait"]<-"Phylum"
COMM_traits.cap1.res[4,"stat"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Chordata_splitted)$statistic
COMM_traits.cap1.res[4,"p value"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Phylum_Chordata_splitted)$p.value
COMM_traits.cap1.res[4,"significant_pairwise_comparison"]<-c("Arthropoda - Tunicata (p=0.01),Echinodermata - Tunicata (p=0.02),Arthropoda - Vertebrata (p=0.007)" )


#Class

cap.COMM$Class<-as.factor(cap.COMM$Class)
kruskal.test(cap.COMM$CAP1,cap.COMM$Class) #  0.01888

COMM_res.class<-FSA::dunnTest(x=cap.COMM$CAP1,g=cap.COMM$Class, method = "bh")
COMM_res.class<-COMM_res.class$res
COMM_res.class<-COMM_res.class$Comparison[COMM_res.class$P.adj<=0.05] # "Ascidiacea - Copepoda" "Copepoda - Teleostei" 

COMM_traits.cap1.res[5,"test"]<-"Kruskal-Wallis"
COMM_traits.cap1.res[5,"Trait"]<-"Class"
COMM_traits.cap1.res[5,"stat"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Class)$statistic
COMM_traits.cap1.res[5,"p value"]<-kruskal.test(cap.COMM$CAP1,cap.COMM$Class)$p.value
COMM_traits.cap1.res[5,"significant_pairwise_comparison"]<-c("Ascidiacea - Copepoda (p=0.005), Copepoda - Teleostei (p=0.01)" )


# fishing vulnerability 

cap.COMM.6<-cap.COMM %>% subset(!Fishbase_Fishing_vulnerability_numeric=="") 
cap.COMM.6$Fishbase_Fishing_vulnerability_numeric<-as.numeric(cap.COMM.6$Fishbase_Fishing_vulnerability_numeric)
COMM_cor.tets<-cor.test(cap.COMM.6$CAP1,cap.COMM.6$Fishbase_Fishing_vulnerability_numeric, method = "spearman") # 0.06

COMM_traits.cap1.res[6,"test"]<-"Spearman_correlation"
COMM_traits.cap1.res[6,"Trait"]<-"Fishing_vulnerability"
COMM_traits.cap1.res[6,"stat"]<-COMM_cor.tets$statistic
COMM_traits.cap1.res[6,"p value"]<-COMM_cor.tets$p.value

# Status

cap.COMM$STATUS<-as.character(cap.COMM$STATUS)
cap.COMM$STATUS[cap.COMM$STATUS==""]<-"Uncertain" # resident species
cap.COMM$STATUS[cap.COMM$STATUS %in% c("NIS-EU","NIS-MED","Cryptogenic")]<-"NIS" 

table(cap.COMM$STATUS)

cap.COMM$STATUS<-as.factor(cap.COMM$STATUS)

shapiro.test(cap.COMM$CAP1)
kruskal.test(CAP1~STATUS, data=cap.COMM)

COMM_res.STATUS<-FSA::dunnTest(x=cap.COMM$CAP1,g=cap.COMM$STATUS, method = "bh")

COMM_traits.cap1.res[7,"test"]<-"Kruskal"
COMM_traits.cap1.res[7,"Trait"]<-"STATUS"
COMM_traits.cap1.res[7,"stat"]<-kruskal.test(cap.COMM$CAP1~cap.COMM$STATUS)$statistic
COMM_traits.cap1.res[7,"p value"]<-kruskal.test(cap.COMM$CAP1~cap.COMM$STATUS)$p.value
COMM_traits.cap1.res[7,"significant_pairwise_comparison"]<-c("NIS - Native (p= 0.0008), NIS - Uncertain (p=0.0013)")

write.csv(COMM_traits.cap1.res, file="Final_files_pubblication/Tables/Traits_tests_outputs_COMM.csv")

#### Plot dbRDA 

### e.2. COMM - PLOT Generate Plot and plot with indicators and traits ----- 
# 
# # Reporte the importance of the axes 
# summ.comm<-summary(COMM.dbrda)
# summ.comm$cont$importance[2, "CAP1"] # importanc eof axes
# summ.comm$cont$importance[2, "MDS1"] # importanc eof axes
# 
# # Plot dbRDA simple by habitats
# # Plot site only dbRDA for NIS 
# # Reporte the importance of the axes 
# 
# 
# # Plot sites and taxa together and colo by mobility
# sites.COMM<-as.data.frame(scores(COMM.dbrda, display = "sites"))
# sites.COMM$sample_name<-rownames(sites.COMM)
# sites.COMM<-left_join(sites.COMM, df[,c("sample_name","Type_habitat_broad")], by ="sample_name")
# 
# dbrda.COMM<-ggplot() +
#   geom_point(data=sites.COMM, 
#              aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
#   # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
#   ggConvexHull:: geom_convexhull(data=sites.COMM,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
#                                                     fill=Type_habitat_broad),alpha=0.0)+
#   
#   scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
#   xlab(paste0("CAP 1 (", round(summ.comm$cont$importance[2,"CAP1"]*100, 1),"%)"))+
#   ylab(paste0("PCoA 1 (", round(summ.comm$cont$importance[2,"MDS1"]*100, 1),"%)"))+
#   theme_minimal()+
#   theme(legend.position="top")+
#   theme(axis.text.x = element_text(size=15),
#         axis.text.y = element_text(size=15),
#         axis.title.x = element_text(size=17),
#         axis.title.y = element_text(size=17),
#         legend.text = element_text(size=13),
#         legend.title = element_text(size=14))
# 
# 
# 
# ggsave(plot=dbrda.COMM, "Final_files_pubblication/Figures/dbRDA_COMM_simple.pdf",
#        dpi=300, width=7,height=6)

# Use the indicator species obtained above at point C.b1.1 

# ## Add indicator to the PCOA plot 

COMM_indic.axis.top<-indicspecs.COMM.sign %>%
  arrange(desc(stat))


COMM_indic.axis.topN <- subset(COMM_indic.axis.top, s.Natural%in%1)
COMM_indic.axis.topA<-subset(COMM_indic.axis.top, s.Artificial%in%1)
COMM_indic.axis.topA<-COMM_indic.axis.topA[1:20,]
COMM_indic.axis.topNA<-rbind(COMM_indic.axis.topA,COMM_indic.axis.topN)


COMM_indic.axis.dbRDA<-cap.COMM %>% subset(Species %in% COMM_indic.axis.top$Species)
COMM_indic.axis.dbRDA.top<-COMM_indic.axis.dbRDA %>% subset(Species %in% COMM_indic.axis.topNA$Species)

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
                               label=Species),color="black",fontface="italic", size=3)+
  
  xlab("CAP 1 (7.6%)")+ylab("PCoA 1 (5.3%)")+
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17))

dbRDA.COMM.species

ggsave(plot=dbRDA.COMM.species, "Final_files_pubblication/Figures/dbRDA_COMM_indicator_species.pdf",
       dpi=300, width= 8, height=9)

# Generate bubblu plot with funcitonal traits 
cap.COMM.summary.Trait1<-cap.COMM.1 %>% group_by(Trait1)%>% summarise(mean_CAP1=mean(CAP1),
                                                                      mean_MDS1=mean(MDS1))
cap.COMM.summary.Trait1$Trait1<-as.character(cap.COMM.summary.Trait1$Trait1)
cap.COMM.summary.Trait1$Trait1[cap.COMM.summary.Trait1$Trait1=="Mobile"] <- "mobile"
cap.COMM.summary.Trait1$Trait1[cap.COMM.summary.Trait1$Trait1=="Sessile"] <- "sessile"

cap.COMM.summary.Trait2<-cap.COMM.2 %>% group_by(Trait2_Worms)%>% 
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

cap.COMM.summary.Fish_depth<-cap.COMM.3 %>%
  group_by(Fishbase_depth_broad)%>% summarise(mean_CAP1=mean(CAP1),
                                              mean_MDS1=mean(MDS1))

cap.COMM.summary.Phylum<-cap.COMM %>% group_by(Phylum_Chordata_splitted)%>%
  subset(Phylum_Chordata_splitted%in% c("Vertebrata","Arthropoda","Echinodermata", "Tunicata"))%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

cap.COMM.summary.Class<-cap.COMM %>% group_by(Class)%>%
  subset(Class%in% c("Ascidiacea","Teleostei","Copepoda"))%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))

cap.COMM.summary.nis<-cap.COMM %>% group_by(STATUS)%>%
  summarise(mean_CAP1=mean(CAP1),
            mean_MDS1=mean(MDS1))


# By each dataset: buble plot
bubble_plot<-ggplot()+
  
  geom_hline(yintercept = 0, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Trait1[!is.na(cap.COMM.summary.Trait1$Trait1),],
             aes(x = mean_CAP1,
                 y = 0),
             color = "darkgoldenrod3",
             alpha=5,
             size = 3)+ 
  # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Trait1[!is.na(cap.COMM.summary.Trait1$Trait1),],aes(
    x = mean_CAP1,
    y = 0,
    label = Trait1),
    size = 6, segment.color = "darkgoldenrod3",
    color="darkgoldenrod3",
    angle=45) +
  
  # Depth non fish 
  geom_hline(yintercept = -0.03, color = "grey80", alpha=0.3) +
  
  geom_point(data=cap.COMM.summary.Trait2[!is.na(cap.COMM.summary.Trait2$Trait2_Worms),], 
             aes(x = mean_CAP1,
                 y = -0.03),
             color = "steelblue",
             alpha=5,
             size = 3)+ # Labels (repelled so they don't overlap)
  ggrepel::geom_text_repel(data=cap.COMM.summary.Trait2[!is.na(cap.COMM.summary.Trait2$Trait2_Worms),],
                           aes(x = mean_CAP1,
                               y = -0.03,
                               label = Trait2_Worms),
                           size = 6, segment.color = "steelblue",
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
                           size = 6, segment.color = "#DC143C",
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
                               label = Phylum_Chordata_splitted),
                           size = 6, segment.color = "grey5",
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
                               label = Class),
                           size = 6, segment.color = "grey35",
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
                               label = STATUS),
                           size = 6, segment.color = "#117733",
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



ggsave(plot=bubble_plot, "Final_files_pubblication/Figures/bubble_plot_traits_COMM.pdf",
       width=14,height = 10, dpi = 300)



## e.3. NAT - Result of beta diversity components: test between nestedness and turnover ----

# Generate a table with the mean jaccard distance, nestedness and turnover for each artificial sample vs each natural sample in each location

final.combine.NAT<-data.frame()

for (i in unique_localities) {
  
  loc<-subset_samples(NAT, Locality %in% i)
  meta.loc<-subset(metaNAT, Locality %in% i)

    
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
  
  final.combine.NAT<-rbind(combine,final.combine.NAT)
  
}


theme_set(theme_minimal())

NAT_beta.components<-ggplot(final.combine.NAT , aes(Group, value, fill=Group))+
  geom_boxplot()+
  ylab("Distance")+xlab("")+
  ggtitle("Beta diversity components",
          "Artificial vs Natural sites - NAT")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))

ggsave(plot=NAT_beta.components, "Final_files_pubblication/Figures/beta.components.boxplot.NAT.pdf",
       dpi=300,width=7, height=6)
ggarrange(NAT_beta.components,COMM_beta.components)



## e.4. NAT - Traits explain the differences between habitats ----
## Dissimilarities between communities inside and outside ports are accounted for by different species traits

tax.NAT<-as.data.frame(tax_table(NAT))

tax.NAT %>% group_by(Trait1) %>% summarise(freq=(n()*100)/354 )
tax.NAT %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/354 )


tax.NAT.trait2<-subset(tax.NAT, !Class %in% c("Teleostei",
                                                    "Elasmobranchii"))

tax.NAT.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/237)

## what percent of mobile species is benthic or pelagic?

tax.NAT.trait2.mob<-subset(tax.NAT.trait2, Trait1 %in% "Mobile")
tax.NAT.trait2.mob %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/143  )


# Vertical distribution of fish
tax.NAT.fish<-subset(tax.NAT, Class %in% c("Teleostei",
                                                              "Elasmobranchii"))

tax.NAT.fish %>% group_by(Fishbase_depth_broad) %>% summarise(freq=(n()*100)/117)


# Test if the CAP1 is influenced by the traits 
cap.NAT<-as.data.frame(scores(NAT.dbrda, display = "species"))
cap.NAT$OTU_seed<-rownames(cap.NAT)

# Create a table where to store the results:

traits.NAT.cap1.res<-as.data.frame(matrix(0,nrow=1,ncol=5))
colnames(traits.NAT.cap1.res)<-c("Trait","test","stat","p value")

# Combine with the trait taxa table

trait.NAT<-as.data.frame(tax_table(NAT)) %>% select(OTU_seed, Trait1,
                                                        Trait2_Worms,
                                                        Fishbase_depth_broad,
                                                        Fishbase_Fishing_vulnerability_numeric,
                                                        Phylum_Chordata_splitted,
                                                        Class,
                                                        Species)


cap.NAT<-left_join(cap.NAT,trait.NAT, by="OTU_seed")%>% column_to_rownames(var="OTU_seed")
cap.NAT$OTU_seed<-rownames(cap.NAT)

# Motility trait 

cap.NAT$Trait1<-as.factor(cap.NAT$Trait1)
cap.NAT.1<-subset(cap.NAT, !Trait1%in%"")
  
shapiro.test(cap.NAT.1$CAP1)

wilcox.test(cap.NAT.1$CAP1~cap.NAT.1$Trait1) # 0.07

traits.NAT.cap1.res[1,"test"]<-"wilcox"
traits.NAT.cap1.res[1,"Trait"]<-"Mobility"
traits.NAT.cap1.res[1,"stat"]<-wilcox.test(cap.NAT.1$CAP1~cap.NAT.1$Trait1)$statistic
traits.NAT.cap1.res[1,"p value"]<-wilcox.test(cap.NAT.1$CAP1~cap.NAT.1$Trait1)$p.value

# Ecology 
cap.NAT$Trait2_Worms<-as.factor(cap.NAT$Trait2_Worms)
cap.NAT.2<-subset(cap.NAT, !Trait2_Worms%in%"")
wilcox.test(cap.NAT.2$CAP1~cap.NAT.2$Trait2_Worms) 

traits.NAT.cap1.res[2,"test"]<-"wilcox"
traits.NAT.cap1.res[2,"Trait"]<-"Distribution in water column"
traits.NAT.cap1.res[2,"stat"]<-wilcox.test(cap.NAT.2$CAP1~cap.NAT.2$Trait2_Worms)$statistic
traits.NAT.cap1.res[2,"p value"]<-wilcox.test(cap.NAT.2$CAP1~cap.NAT.2$Trait2_Worms)$p.value

# Ecology Fish 
cap.NAT$Fishbase_depth_broad<-as.factor(cap.NAT$Fishbase_depth_broad)

cap.NAT.3<-cap.NAT %>% subset(!Fishbase_depth_broad=="") 
kruskal.test(cap.NAT.3$CAP1~cap.NAT.3$Fishbase_depth_broad) # 0.1498 

traits.NAT.cap1.res[3,"test"]<-"Kruskal-Wallis"
traits.NAT.cap1.res[3,"Trait"]<-"Distribution in water column (Fishes)"
traits.NAT.cap1.res[3,"p value"]<-kruskal.test(cap.NAT.3$CAP1~cap.NAT.3$Fishbase_depth_broad)$p.value
traits.NAT.cap1.res[3,"stat"]<-kruskal.test(cap.NAT.3$CAP1~cap.NAT.3$Fishbase_depth_broad)$statistic


# Phylum
cap.NAT$Phylum_Chordata_splitted<-as.factor(cap.NAT$Phylum_Chordata_splitted)
kruskal.test(cap.NAT$CAP1,cap.NAT$Phylum_Chordata_splitted) # 0.0007298

res.NAT.phyl<-FSA::dunnTest(x=cap.NAT$CAP1,g=cap.NAT$Phylum_Chordata_splitted, method = "bh")
res.NAT.phyl<-res.NAT.phyl$res
res.NAT.phyl<-res.NAT.phyl$Comparison[res.NAT.phyl$P.adj<=0.05] # Arthropoda vs Vertebrata

traits.NAT.cap1.res[4,"test"]<-"Kruskal-Wallis"
traits.NAT.cap1.res[4,"Trait"]<-"Phylum"
traits.NAT.cap1.res[4,"stat"]<-kruskal.test(cap.NAT$CAP1,cap.NAT$Phylum_Chordata_splitted)$statistic
traits.NAT.cap1.res[4,"p value"]<-kruskal.test(cap.NAT$CAP1,cap.NAT$Phylum_Chordata_splitted)$p.value
traits.NAT.cap1.res[4,"significant_pairwise_comparison"]<-"Arthropoda - Vertebrata (p=0.0003)"


# Class
cap.NAT$Class<-as.factor(cap.NAT$Class)
kruskal.test(cap.NAT$CAP1,cap.NAT$Class) # 0.0108

res.NAT.class<-FSA::dunnTest(x=cap.NAT$CAP1,g=cap.NAT$Class, method = "bh")
res.NAT.class<-res.NAT.class$res
res.NAT.class<-res.NAT.class$Comparison[res.class$P.adj<=0.05] # Copepoda vs Teleostei Ascidiacea vs Copepoda

traits.NAT.cap1.res[5,"test"]<-"Kruskal-Wallis"
traits.NAT.cap1.res[5,"Trait"]<-"Class"
traits.NAT.cap1.res[5,"stat"]<-kruskal.test(cap.NAT$CAP1,cap.NAT$Class)$statistic
traits.NAT.cap1.res[5,"p value"]<-kruskal.test(cap.NAT$CAP1,cap.NAT$Class)$p.value
traits.NAT.cap1.res[5,"significant_pairwise_comparison"]<-"Copepoda - Teleostei (p=0.0002), Ascidiacea - Copepoda (p=0.05)"


# Fishing vulnerability

cap.NAT.6<-cap.NAT %>% subset(!Fishbase_Fishing_vulnerability_numeric=="") 
cap.NAT.6$Fishbase_Fishing_vulnerability_numeric<-as.numeric(cap.NAT.6$Fishbase_Fishing_vulnerability_numeric)
cor.tets<-cor.test(cap.NAT.6$CAP1,cap.NAT.6$Fishbase_Fishing_vulnerability_numeric, method = "spearman")

traits.NAT.cap1.res[6,"test"]<-"Spearman_correlation"
traits.NAT.cap1.res[6,"Trait"]<-"Fishing_vulnerability"
traits.NAT.cap1.res[6,"stat"]<-cor.tets$statistic
traits.NAT.cap1.res[6,"p value"]<-cor.tets$p.value

write.csv(traits.NAT.cap1.res, file="Final_files_pubblication/Tables/Traits_tests_outputs_NAT.csv")


## e.4.Plot -NAT -----
# Reporte the importance of the axes 
summ.nat<-summary(NAT.dbrda)
summ.nat$cont$importance[2, "CAP1"] # importanc eof axes
summ.nat$cont$importance[2, "MDS1"] # importanc eof axes

# Plot dbRDA simple by habitats
# Plot site only dbRDA for NIS 
# Reporte the importance of the axes 


# Plot sites and taxa together and colo by mobility
sites.NAT<-as.data.frame(scores(NAT.dbrda, display = "sites"))
sites.NAT$sample_name<-rownames(sites.NAT)
sites.NAT<-left_join(sites.NAT, metaNAT[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.NAT<-ggplot() +
  geom_point(data=sites.NAT, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.NAT,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                     fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
  xlab(paste0("CAP 1 (", round(summ.nat$cont$importance[2,"CAP1"]*100, 1),"%)"))+
  ylab(paste0("PCoA 1 (", round(summ.nat$cont$importance[2,"MDS1"]*100, 1),"%)"))+
  ggtitle("NAT")+
  theme_minimal()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))


## e.5. UNCERTAIN - Result of beta diversity components: test between nestedness and turnover ----

# Generate a table with the mean jaccard distance, nestedness and turnover for each artificial sample vs each natural sample in each location

final.combine.UNCERT<-data.frame()

for (i in unique_localities) {
  
  loc<-subset_samples(UNCERT, Locality %in% i)
  meta.loc<-subset(metaUNCERT, Locality %in% i)
  
  
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
  
  final.combine.UNCERT<-rbind(combine,final.combine.NAT)
  
}


theme_set(theme_minimal())

UNC_beta.components<-ggplot(final.combine.UNCERT , aes(Group, value, fill=Group))+
  geom_boxplot()+
  ylab("Distance")+xlab("")+
  ggtitle("Beta diversity components",
          "Artificial vs Natural sites - UNC")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))


## e.6. UNCERTAIN - Traits explains differences between habitats ####

## Dissimilarities between communities inside and outside ports are accounted for by different species traits

tax.UNC<-as.data.frame(tax_table(UNCERT))

tax.UNC %>% group_by(Trait1) %>% summarise(freq=(n()*100)/132 )
tax.UNC %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/132 )


tax.UNC.trait2<-subset(tax.UNC, !Class %in% c("Teleostei",
                                              "Elasmobranchii"))

tax.UNC.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/101)

## what percent of mobile species is benthic or pelagic?

tax.UNC.trait2.mob<-subset(tax.UNC.trait2, Trait1 %in% "Mobile")
tax.UNC.trait2.mob %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/39  )


# Vertical distribution of fish
tax.UNC.fish<-subset(tax.UNC, Class %in% c("Teleostei",
                                       "Elasmobranchii"))

tax.UNC.fish %>% group_by(Fishbase_depth_broad) %>% summarise(freq=(n()*100)/31)


# Test if the CAP1 is influenced by the traits 
cap.UNC<-as.data.frame(scores(UNCERT.dbrda, display = "species"))
cap.UNC$OTU_seed<-rownames(cap.UNC)

# Create a table where to store the results:

traits.UNC.cap1.res<-as.data.frame(matrix(0,nrow=1,ncol=5))
colnames(traits.UNC.cap1.res)<-c("Trait","test","stat","p value")

# Combine with the trait taxa table

trait.UNC<-as.data.frame(tax_table(UNCERT)) %>% select(OTU_seed, Trait1,
                                                    Trait2_Worms,
                                                    Fishbase_depth_broad,
                                                    Fishbase_Fishing_vulnerability_numeric,
                                                    Phylum_Chordata_splitted,
                                                    Class,
                                                    Species)


cap.UNC<-left_join(cap.UNC,trait.UNC, by="OTU_seed")%>% column_to_rownames(var="OTU_seed")
cap.UNC
shapiro.test(cap.UNC$CAP1)

# Motility trait 

cap.UNC$Trait1<-as.factor(cap.UNC$Trait1)
cap.UNC.1<-subset(cap.UNC, !Trait1%in%"")

shapiro.test(cap.UNC.1$CAP1)

wilcox.test(cap.UNC.1$CAP1~cap.UNC.1$Trait1) # 0.3133

traits.UNC.cap1.res[1,"test"]<-"wilcox"
traits.UNC.cap1.res[1,"Trait"]<-"Mobility"
traits.UNC.cap1.res[1,"stat"]<-wilcox.test(cap.UNC.1$CAP1~cap.UNC.1$Trait1)$statistic
traits.UNC.cap1.res[1,"p value"]<-wilcox.test(cap.UNC.1$CAP1~cap.UNC.1$Trait1)$p.value

# Ecology 
cap.UNC$Trait2_Worms<-as.factor(cap.UNC$Trait2_Worms)
cap.UNC.2<-subset(cap.UNC, !Trait2_Worms%in%"")
wilcox.test(cap.UNC.2$CAP1~cap.UNC.2$Trait2_Worms) 

traits.UNC.cap1.res[2,"test"]<-"wilcox"
traits.UNC.cap1.res[2,"Trait"]<-"Distribution in water column"
traits.UNC.cap1.res[2,"stat"]<-wilcox.test(cap.UNC.2$CAP1~cap.UNC.2$Trait2_Worms)$statistic
traits.UNC.cap1.res[2,"p value"]<-wilcox.test(cap.UNC.2$CAP1~cap.UNC.2$Trait2_Worms)$p.value

# Ecology Fish 
cap.UNC$Fishbase_depth_broad<-as.factor(cap.UNC$Fishbase_depth_broad)

cap.UNC.3<-cap.UNC %>% subset(!Fishbase_depth_broad=="") 
kruskal.test(cap.UNC.3$CAP1~cap.UNC.3$Fishbase_depth_broad) # 0.5418 

traits.UNC.cap1.res[3,"test"]<-"Kruskal-Wallis"
traits.UNC.cap1.res[3,"Trait"]<-"Distribution in water column (Fishes)"
traits.UNC.cap1.res[3,"p value"]<-kruskal.test(cap.UNC.3$CAP1~cap.UNC.3$Fishbase_depth_broad)$p.value
traits.UNC.cap1.res[3,"stat"]<-kruskal.test(cap.UNC.3$CAP1~cap.UNC.3$Fishbase_depth_broad)$statistic


# Phylum
cap.UNC$Phylum_Chordata_splitted<-as.factor(cap.UNC$Phylum_Chordata_splitted)
kruskal.test(cap.UNC$CAP1,cap.UNC$Phylum_Chordata_splitted) # 0.3134

traits.UNC.cap1.res[4,"test"]<-"Kruskal-Wallis"
traits.UNC.cap1.res[4,"Trait"]<-"Phylum"
traits.UNC.cap1.res[4,"stat"]<-kruskal.test(cap.UNC$CAP1,cap.UNC$Phylum_Chordata_splitted)$statistic
traits.UNC.cap1.res[4,"p value"]<-kruskal.test(cap.UNC$CAP1,cap.UNC$Phylum_Chordata_splitted)$p.value


# Class
cap.UNC$Class<-as.factor(cap.UNC$Class)
kruskal.test(cap.UNC$CAP1,cap.UNC$Class) # 0.2302

traits.UNC.cap1.res[5,"test"]<-"Kruskal-Wallis"
traits.UNC.cap1.res[5,"Trait"]<-"Class"
traits.UNC.cap1.res[5,"stat"]<-kruskal.test(cap.UNC$CAP1,cap.UNC$Class)$statistic
traits.UNC.cap1.res[5,"p value"]<-kruskal.test(cap.UNC$CAP1,cap.UNC$Class)$p.value


# Fishing vulnerability

cap.UNC.6<-cap.UNC %>% subset(!Fishbase_Fishing_vulnerability_numeric=="") 
cap.UNC.6$Fishbase_Fishing_vulnerability_numeric<-as.numeric(cap.UNC.6$Fishbase_Fishing_vulnerability_numeric)
cor.tets<-cor.test(cap.UNC.6$CAP1,cap.UNC.6$Fishbase_Fishing_vulnerability_numeric, method = "spearman")

traits.UNC.cap1.res[6,"test"]<-"Spearman_correlation"
traits.UNC.cap1.res[6,"Trait"]<-"Fishing_vulnerability"
traits.UNC.cap1.res[6,"stat"]<-cor.tets$statistic
traits.UNC.cap1.res[6,"p value"]<-cor.tets$p.value

write.csv(traits.UNC.cap1.res, file="Final_files_pubblication/Tables/Traits_tests_outputs_UNCERTAIN.csv")

## e.6.Plot -UNC -----
# Reporte the importance of the axes 
summ.unc<-summary(UNCERT.dbrda)
summ.unc$cont$importance[2, "CAP1"] # importanc eof axes
summ.unc$cont$importance[2, "MDS1"] # importanc eof axes

# Plot dbRDA simple by habitats
# Plot site only dbRDA for NIS 
# Reporte the importance of the axes 


# Plot sites and taxa together and colo by mobility
sites.UNC<-as.data.frame(scores(UNCERT.dbrda, display = "sites"))
sites.UNC$sample_name<-rownames(sites.UNC)
sites.UNC<-left_join(sites.UNC, metaUNCERT[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.UNC<-ggplot() +
  geom_point(data=sites.UNC, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.UNC,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                    fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"
  xlab(paste0("CAP 1 (", round(summ.unc$cont$importance[2,"CAP1"]*100, 1),"%)"))+
  ylab(paste0("PCoA 1 (", round(summ.unc$cont$importance[2,"MDS1"]*100, 1),"%)"))+
  theme_minimal()+ggtitle("Uncertain")+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))

ggarrange(dbrda.NAT,dbrda.NIS,dbrda.UNC, nrow=1, common.legend = TRUE)

# g. Generalities of NIS -----

taxNIS<-as.data.frame(tax_table(NIS))

taxNIS %>% group_by(STATUS)%>% summarise(total=n()) 

## h. Proportion of NIS in the community ----

### h.1. barplot with proportion NIS vs NAT inside and outside by locality ----

counts<-as.data.frame(otu_table(combined_taxa_ps_summer.sp))

counts$OTU_seed<-rownames(counts)

tax.all<-as.data.frame(tax_table(combined_taxa_ps_summer.sp))

tax.count<-left_join(counts,tax.all,by="OTU_seed")

rownames(tax.count)<-tax.count$OTU_seed

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

tax.countNIS.all<-subset(tax.count, STATUS %in% c("NIS-EU","NIS-MED","Cryptogenic")) # here I include all the NIS both from and not the Mediterranean Sea

total.NIS.all<-get_total(tax.countNIS.all,"total_number_NIS")

proportion.table<-list(total,total.NIS.all) %>% reduce(left_join, by="sample_name")

proportion.table<-proportion.table %>%  
  mutate(percentage_NIS_all=
           (proportion.table$total_number_NIS/proportion.table$total_number_species)*100)

# Mean proportion of NIS in the samples

mean(proportion.table$percentage_NIS_all) # 8.677012
sd(proportion.table$percentage_NIS_all) #5.216064

# SE
sd(proportion.table$percentage_NIS_all) / sqrt(length(proportion.table$percentage_NIS_all))

# Bind the sample data to this table
df.NIS<-sample_data(combined_taxa_ps_summer.sp)%>% data.frame()

meta.with.proportions<-left_join(df.NIS,proportion.table,by="sample_name")

rownames(meta.with.proportions)<-meta.with.proportions$sample_name

NIS.prop.in.Locs<-meta.with.proportions %>% group_by(Locality) %>% 
  summarise(mean.NIS<-mean(percentage_NIS_all),
            SE.NIS<-sd(percentage_NIS_all) / sqrt(length(percentage_NIS_all))
) 


# write.csv(NIS.prop.in.Locs,file="Final_files_pubblication/Tables/NIS.proportion.in.localities.csv")

# Save also table with number of NIS per sample and proportion 

write.csv(meta.with.proportions [,c("sample_name","Locality",
                                 "Type_habitat_broad","total_number_species",
                                 "total_number_NIS",
                                 "percentage_NIS_all")],
          file="Final_files_pubblication/Tables/Proportion_NIS_across_all_samples.csv")

# Number of NIS per locality
meta.with.proportions %>% group_by()


# PLot the prcentage of species (Non-Indigenous) by habitat type and locality 

View(meta.with.proportions %>%
  group_by(Type_habitat_broad,Locality) %>%
  summarise(NIS_richness=mean(total_number_NIS),
            SE.NIS.richness=sd(total_number_NIS) / sqrt(length(total_number_NIS))))

rownames(summary_NIS.bysample)<-summary_NIS.bysample$sample_name


summary_NIS_proportion<-meta.with.proportions %>% 
  group_by(Type_habitat_broad,Locality) %>% 
  
  summarise(
    mean_NIS.ALL=mean(percentage_NIS_all),
            sd_NIS.ALL=sd(percentage_NIS_all),
            SE.NIS=sd(percentage_NIS_all) / sqrt(length(percentage_NIS_all))) %>% 
  as.data.frame()

meta.with.proportions %>% group_by(Type_habitat_broad) %>% 
  summarise(mean_NIS.ALL=mean(percentage_NIS_all),
            sd_NIS.ALL=sd(percentage_NIS_all),
            SE.NIS=sd(percentage_NIS_all) / sqrt(length(percentage_NIS_all))) %>% 
  as.data.frame()

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

ggsave(plot=NIS_all,"Final_files_pubblication/Figures/proportion_NIS_by_habitat.pdf",
       dpi=300,width=10,height = 9)

### h.2. test on the proportions ----

shapiro.test(meta.with.proportions$percentage_NIS_all)

meta.with.proportions$Type_habitat_broad<-factor(meta.with.proportions$Type_habitat_broad, 
                                                 levels=c("Artificial","Natural"))

test<- t.test(meta.with.proportions$percentage_NIS_all~meta.with.proportions$Type_habitat_broad, alternative="greater") # t = 6.5877, df = 54.294, p-value = 9.402e-09 
p_val<-test$p.value
stat<-test$statistic
test<-"t.test"
res.prop<-data.frame(p_val,stat,test)
rownames(res.prop)<-"All_localities"


## GENERALIZED LINEAR MIXED MODEL  ----
meta.with.proportions$percentage_NIS_all_01<-(meta.with.proportions$percentage_NIS_all)/100
shapiro.test(meta.with.proportions$total_number_NIS / meta.with.proportions$total_number_species)

mod.NIS.prop <- glmmTMB(
  total_number_NIS / total_number_species ~
    Type_habitat_broad+ (1| Locality),
  weights = total_number_species,
  family = binomial(),
  data = meta.with.proportions
)


summary(mod.NIS.prop)

confint(mod.NIS.prop)

res<-DHARMa::simulateResiduals(mod.NIS.prop)
DHARMa::testDispersion(res)


# By locality - Artificial -- t test and wilcosx test

unique_localities<-unique(sample_data(NIS)$Locality)

result.data.finAll<-data.frame()

for (i in unique_localities) {
  
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
    test<-"t.test"
    result.data<-data.frame(p_val,stat,test,i)
    rownames(result.data)<-i
    result.data.finAll<-rbind(result.data,result.data.finAll)
  } else {
    
    test<-  wilcox.test(loc$percentage_NIS_all~loc$Type_habitat_broad, alternative="greater")
    p_val<-test$p.value
    stat<-test$statistic
    test<-"wilcox"
    result.data<-data.frame(p_val,stat,test,i)
    rownames(result.data)<-i
    result.data.finAll<-rbind(result.data,result.data.finAll)
  }
}


result.data.finAll$adusted.pval.holm<-p.adjust(result.data.finAll$p,method="holm")
result.data.finAll$adusted.pval.BH<-p.adjust(result.data.finAll$p,method="BH")

res.prop$i<-"All_localities"

result.data.finAll<-rbind(res.prop,result.data.finAll)

write.csv(result.data.finAll, file="Final_files_pubblication/Tables/proportion_test_NIS_localities.csv")

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
ggsave(plot=boxplot.proportions, "Final_files_pubblication/Figures/NIS_proportion_across_localities_no_habitat_distinction.pdf",
       width=8,height=6,dpi=300)

# Not significant difference in prop of NIS across localities, neither for artificial or natural sites

## i. Greater richness of NIS in ports than outside -----

### i.1. NIS alpha diversity is greater inside than outside + boxplot ----

richNIS<-estimate_richness(NIS, measures="Observed")

richNIS$sample_name<-rownames(richNIS)

metaNIS<-sample_data(NIS) %>% data.frame()

meta.rich.NIS<-left_join(metaNIS,richNIS,by="sample_name")

#### NIS.a.1. Alpha diversity (richness) all localities together ----

# Test different richness in ports and in natural habitats

# Perform one tail test - higher richness in ports
NIS.alpha.test.sp.ot<-alpha_tests.new(NIS, levels = c("Natural","Artificial"))   # invert the order to check wheter the diversity of NIS is greater in ports
rownames(NIS.alpha.test.sp.ot)<-"All locations"
print(NIS.alpha.test.sp.ot)
NIS.alpha.test.sp.ot$test_type<-"one_tail_higher_Port"

#### NIS.a.2. Alpha diversity by locality ----

# Higher richness in ports habitats 
unique_localities<-unique(metaNIS$Locality)

result.locality.NIS.higher.Artificial<-data.frame()

for (locality in unique_localities[c(1,3:12)]){
  
  ps.sub<-subset_samples(NIS, Locality == locality)
  
  result.test<-alpha_tests.new(ps.sub, levels=c("Natural","Artificial"))
  rownames(result.test)<-locality
  
  result.locality.NIS.higher.Artificial<-rbind(result.locality.NIS.higher.Artificial,result.test)
  
}

print(result.locality.NIS.higher.Artificial)

result.locality.NIS.higher.Artificial$adjusted.pvalue.holm<-p.adjust(result.locality.NIS.higher.Artificial$pvalue,method = "holm")
result.locality.NIS.higher.Artificial$adjusted.pvalue.BH<-p.adjust(result.locality.NIS.higher.Artificial$pvalue,method = "BH")

# Bastia repeat it separately 

bastia<-subset_samples(NIS, Locality == "Bastia")

rich<-estimate_richness(bastia, measures="Observed")

rich$sample_name<-rownames(rich)

meta<-sample_data(bastia) %>% data.frame()

meta.rich<-left_join(meta,rich,by="sample_name")

# Run wilcox

meta.rich$Type_habitat_broad<-  factor(meta.rich$Type_habitat_broad,
                                       levels = c("Natural","Artificial"))


bastia.result<- wilcox.test(Observed~Type_habitat_broad, data=meta.rich, paired = FALSE, exact=FALSE, alternative=c("less"))

result.locality.NIS.higher.Artificial[11,"pvalue"]<-bastia.result$p.value
result.locality.NIS.higher.Artificial[11,"stat"]<-bastia.result$stat
result.locality.NIS.higher.Artificial[11,"test"]<-"wilcox"
result.locality.NIS.higher.Artificial[11,"variable.tested"]<-"Habitat"
rownames(result.locality.NIS.higher.Artificial)[11]<-"Bastia"

result.locality.NIS.higher.Artificial$test_type<-"one_tail_higher_Port"

combined_NIS_div<-rbind(NIS.alpha.test.sp.ot,result.locality.NIS.higher.Artificial)
combined_NIS_div$Phlyum<-"All_Phyla"



#### NIS.a.3. Alpha diversity by phylum  ----
unique_phylum<-unique(as.data.frame(tax_table(NIS))$Phylum_Chordata_splitted)

# One tail

NIS.alpha.test.phyla.ot<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(NIS, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  result.test<-alpha_tests.new(ps.sub.phyl, levels = c("Natural","Artificial"))
  
  result.test$Phlyum<-phylum
  
  NIS.alpha.test.phyla.ot<-rbind(NIS.alpha.test.phyla.ot,result.test)
}

NIS.alpha.test.phyla.ot
NIS.alpha.test.phyla.ot$test_type<-"one_tail_higher_Port"

combined_NIS_div<-rbind(combined_NIS_div, NIS.alpha.test.phyla.ot)

# Save table
write.csv(combined_NIS_div,file="Final_files_pubblication/Tables/Alpha_diversity_tests_outputs_NIS.csv")


#### NIS.a.4. Plot richness across all dataset ---- 

NIS.richness.allLocs<-ggplot(meta.rich.NIS,aes(x = Type_habitat_broad, y=Observed, fill=Type_habitat_broad,
                                               color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point()+
  xlab("")+ylab("Richness of species")+
  # labs(title = "Richness of Species by habitat type")+
  # stat_compare_means(method = "t.test",label = "p.signif",label.x = 1.5,
  #                    label.y = 100)+   # show only *, **, etc.
  # #hide.ns = TRUE) +     # hide 'ns' when not significant
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold1"))+
  #facet_grid(~Locality, scale="free_x")+
  # scale_x_discrete(labels=c("Artificial"="Inside port", "Natural"="Outside port"))+
  # annotate("text", x=1.05, y=105,size=5, label= "Wilcoxon test, P value = 0.009", )+
  fig+theme(axis.text.x = element_blank(),
            legend.position = "top",
            legend.text=element_text(size=10),
            legend.title=element_text(size=11))

NIS.richness.allLocs

# Combine together teh boxplots of richness

alpha.div<-ggarrange(COMM.richness.allLocs+ggtitle("COMM"),
          NAT.richness.allLocs+ggtitle("NAT"),
          NIS.richness.allLocs+ggtitle("NIS"),
          common.legend = TRUE,
          nrow=1)

ggsave(plot=alpha.div, "Final_files_pubblication/Figures/alpha_diversity_boxplot_COOM.NAT.NIS.pdf",
       dpi=300,width=15, height =6)

#### GENERALIZED LINEAR MIXE DMODEL - NIS -----

# Perform generalized linear mixed models for the whole dataset and repeated for each phyla
meta.rich.NIS

shapiro.test(meta.rich.NIS$Observed)

# Set the levels to have as reference the Artificial - Port one
meta.rich.NIS$Type_habitat_broad<-  factor(meta.rich.NIS$Type_habitat_broad,
                                            levels = c("Natural", "Artificial"))


mod.NIS<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
             family=nbinom2,
             data=meta.rich.NIS)


# mod2<-glmmTMB(Observed~Type_habitat_broad*Locality, 
#              family=nbinom2,
#              data=COMM.meta.rich) # this overfit so better to not use it even if the AIC is th elowest. Too little replicates per each habotat to rely on it.
# 
# 
# mod3<-glmmTMB(Observed~Type_habitat_broad+(Type_habitat_broad|Locality), 
#               family=nbinom2,
#               data=COMM.meta.rich)


summary.mod.NIS<-summary(mod.NIS)
confint(mod.NIS)
# Wald intervals (fast)

# summarmod# summary(mod2)
# summary(mod3)

res<-DHARMa::simulateResiduals(mod.NIS)
DHARMa::testDispersion(res)

## get he coefficients

summary.mod.NIS<-as.data.frame(summary.mod.NIS$coefficients$cond)[2,]
summary.mod.NIS$Phylum<-"All"
summary.mod.NIS$ref.habitat<-"Port"

# By each phylum 
unique_phylum<-unique(as.data.frame(tax_table(NIS))$Phylum_Chordata_splitted)

glmm.phyl.NIS<-data.frame()

for (phylum in unique_phylum){
  
  ps.sub.phyl<-subset_taxa(NIS, Phylum_Chordata_splitted == phylum) #repeat the same loop by changing everytime the subset 
  
  rich<-estimate_richness(ps.sub.phyl, measures="Observed")
  
  rich$sample_name<-rownames(rich)
  
  meta<-sample_data(ps.sub.phyl) %>% data.frame()
  
  meta.rich<-left_join(meta,rich,by="sample_name")
  
  result.test<-glmmTMB(Observed~Type_habitat_broad + (1|Locality), 
                       family=nbinom2,
                       data=meta.rich)
  
  
  
  result.test<-as.data.frame(summary(result.test)$coefficients$cond)
  result.test<-result.test[2,]
  result.test$Phylum<-phylum
  
  
  glmm.phyl.NIS<-rbind(glmm.phyl.NIS,result.test)
}

glmm.phyl.NIS

glmm.phyl.NIS$ref.habitat<-"Port"

glmm.NIS<-rbind(summary.mod.NIS,glmm.phyl.NIS)

write.csv(glmm.NIS,file="Final_files_pubblication/Tables/glmmTMB_NIS.csv")




### i.1.NIS.b - Test NIS richness across localities ----
# Test and generate boxplot throughout localities splittniig inside and outside

## Only inside 
rich.NIS.in<-meta.rich.NIS[meta.rich.NIS$Type_habitat_broad %in% "Artificial",]
rich.NIS.in %>% group_by(Locality) %>% summarise(mean_NIS=mean(Observed),
                                                 sd_NIS=sd(Observed))

shapiro.test(rich.NIS.in$Observed)

kruskal.test(Observed~Locality, data=rich.NIS.in) # 0.01337

tbl.rni<-as.data.frame(dunn.test::dunn.test(x=rich.NIS.in$Observed,g=rich.NIS.in$Locality, method = "bh"))
tbl.rni$Habitat<-"Artificial"

## Only outside

rich.NIS.out<-meta.rich.NIS[meta.rich.NIS$Type_habitat_broad %in% "Natural",]

rich.NIS.out %>% group_by(Locality) %>% summarise(mean_NIS=mean(Observed),
                                                           sd_NIS=sd(Observed))


shapiro.test(rich.NIS.out$Observed)

kruskal.test(Observed~Locality, data=rich.NIS.out) # 0.01588

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

ggsave(plot=NIS.rich.plot,"Final_files_pubblication/Figures/NIS_richness_across_localities.pdf",
       width=11 ,height = 7,dpi = 300)


### l.1. Venn diagram  -----
venn.plot(NIS)

## m. Beta diversity of NIS ---- 

NIS_PCoA<-ordinate(NIS, method = "PCoA", distance="jaccard") 

theme_set(theme_minimal())

NIS_PCoA.plot<-plot_ordination(NIS, NIS_PCoA, axes=c(1,2),color= "Locality", 
                                shape="Type_habitat_broad")+
  geom_point(size=4, alpha=1, aes(color=Locality))+
  scale_color_manual(values = loc.cols)+ 
  labs(title = paste0(" PCoA - NIS"))+fig+
  theme(legend.position = "top")


ggsave(plot=NIS_PCoA.plot, "Final_files_pubblication/Figures/Unconstrained_MDS_loc_habitat_LEGEND.pdf",
      dpi=300,width=12,height = 6)

#### m.1. Beta diversity test on all localities together ----

jaccard.NIS<-distance(NIS, method = "jaccard", binary=TRUE)
metaNIS<- sample_data(NIS) %>% data.frame() # extract the dataframe
metaNIS[,"Type_habitat_broad"]<-as.factor(metaNIS[,"Type_habitat_broad"])

# check beta dispersion - betadisperser test for homogeneity of variance

dispr.NIS<-betadisper(jaccard.NIS,metaNIS$Type_habitat_broad) 
permutest(dispr.NIS) # P value 0.18   
boxplot(dispr.NIS)

#### NIS - dbRDA  ----

jaccard.NIS

NIS_otu<-as.data.frame(otu_table(NIS)) 

NIS.dbrda<-capscale(jaccard.NIS~Type_habitat_broad +
                      Condition(Locality), # set locality as condition since it not the focus of the analysis 
                    distance="jaccard", data=metaNIS, 
                    comm = as.data.frame(t(NIS_otu)))


anova.cca(NIS.dbrda)
anova.cca(NIS.dbrda, by ="terms",permutations = 999)
RsquareAdj(NIS.dbrda) 


# permutest result output

permutest.table <- rbind(permutest(dispr.COMM)$tab,
      permutest(dispr.NAT)$tab,
      permutest(dispr.NIS)$tab,
      permutest(dispr.UNC)$tab)

permutest.table$Category<-c("COMM","COMM",
                            "NAT","NAT",
                            "NIS","NIS",
                            "Uncertain","Uncertain")

dbRDA.table<-rbind(anova.cca(COMM.dbrda),
                   anova.cca(NAT.dbrda),
                   anova.cca(NIS.dbrda),
                   anova.cca(UNCERT.dbrda))

dbRDA.table$Category<-c("COMM","COMM",
                        "NAT","NAT",
                        "NIS","NIS",
                        "Uncertain","Uncertain")

write.csv(dbRDA.table, file="Final_files_pubblication/Tables/dbRDA_results_table.csv")
write.csv(permutest.table, file="Final_files_pubblication/Tables/permutes_results_table.csv")

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
          "Artificial vs Natural sites - NIS")+
  scale_fill_manual(values=c("grey15","grey50","grey90"))+fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18),
        title = element_text(size=18))


components.plots<-ggarrange(NAT_beta.components, NIS.plot.components, UNC_beta.components, nrow=1, common.legend = TRUE)

ggsave(plot=components.plots, "Final_files_pubblication/Figures/beta_components_plots_NAT.NIS.UNC.pdf",
       dpi=300, width=16, height=6)


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

wilcox.test(Jaccard_distance~Type_habitat_broad, data=dispr.NIS.dist, paired=FALSE, alternative=c("less")) #0.3837

# Plot

NIS_homogen.plot<-ggplot(dispr.NIS.dist, aes(x=Type_habitat_broad,y=Jaccard_distance, 
                                               fill=Type_habitat_broad, color=Type_habitat_broad))+
  geom_boxplot()+
  geom_point(size=2)+theme_set(theme_minimal())+
  scale_fill_manual(values=c("darkblue","yellow2"))+
  scale_color_manual(values=c("dodgerblue","gold3"))+
  ylab("Jaccard distance")+
  xlab("")+
  ggtitle("NIS")+
  # scale_x_discrete(labels = c("Artificial"="Inside port","Natural"="Outside port"))+
  fig+
  #annotate("text", x=1.05, y=0.75,size=5, label= "T test, P value = 0.001", )+
  theme(legend.position = "top",
        axis.text.x = element_blank(),
        axis.text.y=element_text(size=15),
        axis.title.y = element_text(size=18))

homogen.plots<-ggarrange(NAT.homogen.plot, NIS_homogen.plot, UNC.homogen.plot, nrow=1, common.legend = TRUE)

ggsave(plot=homogen.plots, "Final_files_pubblication/Figures/biotic_homogenization_plots_NAT.NIS.UNC.pdf",
       dpi=300, width=16, height=6)
## n.2. NIS Heatmap ----

# Step 1: Detemine the number of replicates per site in which the NIS is found

#obtain the samples names in each locality

site.list<-metaNIS%>% group_by(Locality) %>% summarise(sites=as.vector(sample_name))

replicates.final<-data.frame(matrix(NA, nrow = 50, ncol = 0))

for (locality in unique(site.list$Locality)) {
  
  Locality<-subset(site.list, Locality %in% locality)
  
  sub<-as.data.frame(tax.countNIS.all[,Locality$sites])
  
  replicates<-as.data.frame(rowSums(sub))
  
  colnames(replicates)<-locality
  
  replicates.final<-cbind(replicates,replicates.final)
  
}

replicates.final$OTU_seed<-rownames(replicates.final)

tax.sub<-taxNIS[,c("OTU_seed","Species","AphiaID","Phylum_Chordata_splitted",
                "collapsed_Markers",
                "Trait1",
                "STATUS")]

replicates.final<-left_join(replicates.final,tax.sub,by ="OTU_seed")

rownames(replicates.final)<-replicates.final$OTU_seed


# add column with "shared","unique artificial" and "unique natural" form the Venn diagram of NIS

NIS.unique.shared<-subset(COMM_unique.shared,STATUS %in% c("NIS-EU","NIS-MED","Cryptogenic"))
NIS.unique.shared$OTU_seed<-rownames(NIS.unique.shared)

replicates.final<-replicates.final %>% mutate(Habitat_Venn=case_when(OTU_seed%in% NIS.unique.shared$OTU_seed[NIS.unique.shared$Venn%in% "Shared"]~"Shared",
                                                                     OTU_seed%in% NIS.unique.shared$OTU_seed[NIS.unique.shared$Venn%in% "Unique_Artificial"]~"Port",
                                                                     OTU_seed%in% NIS.unique.shared$OTU_seed[NIS.unique.shared$Venn%in% "Unique_Natural"]~"Natural"))
                                                                     

write.csv2(replicates.final,file="Final_files_pubblication/Tables/All_NIS_replicates_localities.csv")

# Make heatmap with species and localities and numb of replicates

# load libraries
library(pheatmap)
library(grid)

# load matrix

mat.to.heatmap<-read_xlsx("Final_files_pubblication/Tables/matrix_NIS_med_sea_to_heatmap_.xlsx") %>%column_to_rownames(var="Species")

# Order all the species by Phylum

mat.to.heatmap<- mat.to.heatmap %>% arrange(Phylum_Chordata_splitted)
# italiacize the species names 

species_names_italics<-lapply(rownames(mat.to.heatmap), function(x)bquote(italic(.(x))))

# Define colors 
heat_color<-colorRampPalette(c("seashell","firebrick"))(100)


# Prepare annotation for row (species)
annotation_row <- data.frame(Habitat=mat.to.heatmap$Habitat,
                             Mobility=mat.to.heatmap$Trait1,
                             STATUS = mat.to.heatmap$STATUS,
                             Phylum= mat.to.heatmap$Phylum_Chordata_splitted)

rownames(annotation_row) <- rownames(mat.to.heatmap)

colnames(mat.to.heatmap)
mat.to.heatmap<-mat.to.heatmap[,-c(13:19)]

desired_order <- c("Port-la-Nouvelle", "Sete_1", "Sete_2", "Seyne-sur-mer", "Toulon", "Saint-Tropez",
                   "Nice", "Calvi", "Ajaccio", "Bonifacio", "Porto-Vecchio", "Bastia")

# Reorder columns of mat.to.heatmap

mat.to.heatmap <- mat.to.heatmap[, desired_order]

# Generate unique colors for Phylum (bright)
phylum_levels <- unique(annotation_row$Phylum)
phylum_colors.sub <- phylum_colors[names(phylum_colors)%in%phylum_levels]
  
# Black and white for Status
status_levels <- unique(annotation_row$STATUS)
status_colors <- setNames(c("black", "white", "grey50")[1:length(status_levels)], status_levels)

# Custom colors for Mobility
mobility_colors <- c("Mobile" = "deeppink", "Sessile" = "greenyellow")

# Custom colors for Habitat
# Blend yellow2 and blue for shared
shared_color <- colorRampPalette(c("yellow", "blue"))(3)[2]
habitat_colors <- c("Port" = "blue", "Natural" = "yellow", "Shared" = shared_color)

# Define annotation colors list
annotation_colors <- list(
  Phylum = phylum_colors.sub,
  STATUS = status_colors,
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


## Dissimilarities between communities inside and outside ports are accounted for by different species traits

NIS_tax<-as.data.frame(tax_table(NIS))

NIS_tax %>% group_by(Trait1) %>% summarise(freq=(n()*100)/50) 
NIS_tax %>% group_by(Trait1_3rd_lev_fish) %>% summarise(freq=(n()*100)/50)

NIS_tax.trait2<-subset(NIS_tax, !Class %in% c("Teleostei",
                                                "Elasmobranchii"))

NIS_tax.trait2 %>% group_by(Trait2_Worms) %>% summarise(freq=(n()*100)/48)# all mobile species without fishes (benthic and zooplankton)

# Vertical distribution of fish - not possible to do becasue only 2 fishes are present among NIS 
# NIS_tax.fish<-subset(NIS_tax, Class %in% c("Teleostei",
#                                              "Elasmobranchii"))


# Test traits and plot dbRDA

# Test if the CAP1 is influenced by the traits 
cap.NIS<-as.data.frame(scores(NIS.dbrda, display = "species"))
cap.NIS$OTU_seed<-rownames(cap.NIS)

# Create a table where to store the results:

NIS_traits.cap1.res<-as.data.frame(matrix(0,nrow=1,ncol=5))
colnames(NIS_traits.cap1.res)<-c("Trait","test","stat","p value")

# Combine with the trait taxa table

trait.NIS<-as.data.frame(tax_table(NIS)) %>% select(OTU_seed, Trait1,
                                                      Trait2_Worms,
                                                      Fishbase_depth_broad,
                                                      Fishbase_Fishing_vulnerability_numeric,
                                                      Phylum_Chordata_splitted,
                                                      Class,
                                                      STATUS,
                                                      Species)


cap.NIS<-left_join(cap.NIS,trait.COMM, by="OTU_seed")%>% column_to_rownames(var="OTU_seed")
cap.NIS$OTU_seed<-rownames(cap.NIS)

# Motility trait 
cap.NIS$Trait1<-as.factor(cap.NIS$Trait1)
cap.NIS.1<-subset(cap.NIS, !Trait1%in%"")

shapiro.test(cap.NIS.1$CAP1)

wilcox.test(cap.NIS.1$CAP1~cap.NIS.1$Trait1) 

NIS_traits.cap1.res[1,"test"]<-"wilcox"
NIS_traits.cap1.res[1,"Trait"]<-"Mobility"
NIS_traits.cap1.res[1,"stat"]<-wilcox.test(cap.NIS.1$CAP1~cap.NIS.1$Trait1)$statistic
NIS_traits.cap1.res[1,"p value"]<-wilcox.test(cap.NIS.1$CAP1~cap.NIS.1$Trait1)$p.value

# Ecology Trait

cap.NIS$Trait2_Worms<-as.factor(cap.NIS$Trait2_Worms)
cap.NIS.2<-subset(cap.NIS, !Trait2_Worms%in%"")
wilcox.test(cap.NIS.2$CAP1~cap.NIS.2$Trait2_Worms) # 5.628e-06

NIS_traits.cap1.res[2,"test"]<-"wilcox"
NIS_traits.cap1.res[2,"Trait"]<-"Distribution in water column"
NIS_traits.cap1.res[2,"stat"]<-wilcox.test(cap.NIS.2$CAP1~cap.NIS.2$Trait2_Worms)$statistic
NIS_traits.cap1.res[2,"p value"]<-wilcox.test(cap.NIS.2$CAP1~cap.NIS.2$Trait2_Worms)$p.value



NIS_traits.cap1.res[3,"test"]<-"Not possible to perform tests, only 2 fishes"
NIS_traits.cap1.res[3,"Trait"]<-"Distribution in water column - Fishes"


#  PHylum

cap.NIS$Phylum_Chordata_splitted<-as.factor(cap.NIS$Phylum_Chordata_splitted)
kruskal.test(cap.NIS$CAP1,cap.NIS$Phylum_Chordata_splitted) # 0.0009871

NIS_res.phyl<-FSA::dunnTest(x=cap.NIS$CAP1,g=cap.NIS$Phylum_Chordata_splitted, method = "bh")
NIS_res.phyl<-NIS_res.phyl$res
NIS_res.phyl<-NIS_res.phyl$Comparison[NIS_res.phyl$P.adj<=0.05] # None

NIS_traits.cap1.res[4,"test"]<-"Kruskal-Wallis"
NIS_traits.cap1.res[4,"Trait"]<-"Phylum"
NIS_traits.cap1.res[4,"stat"]<-kruskal.test(cap.NIS$CAP1,cap.NIS$Phylum_Chordata_splitted)$statistic
NIS_traits.cap1.res[4,"p value"]<-kruskal.test(cap.NIS$CAP1,cap.NIS$Phylum_Chordata_splitted)$p.value


#Class

cap.NIS$Class<-as.factor(cap.NIS$Class)
kruskal.test(cap.NIS$CAP1,cap.NIS$Class) #  0.1357

NIS_traits.cap1.res[5,"test"]<-"Kruskal-Wallis"
NIS_traits.cap1.res[5,"Trait"]<-"Class"
NIS_traits.cap1.res[5,"stat"]<-kruskal.test(cap.NIS$CAP1,cap.NIS$Class)$statistic
NIS_traits.cap1.res[5,"p value"]<-kruskal.test(cap.NIS$CAP1,cap.NIS$Class)$p.value


# fishing vulnerability 


NIS_traits.cap1.res[6,"test"]<-"Not posisble to perform test, only 2 fishes"
NIS_traits.cap1.res[6,"Trait"]<-"Fishing_vulnerability"

write.csv(NIS_traits.cap1.res, file="Final_files_pubblication/Tables/Traits_tests_outputs_NIS.csv")

#### Plot dbRDA 

# Plot site only dbRDA for NIS 

summ.nis<-summary(NIS.dbrda)
summ.nis$cont$importance[2,"CAP1"]
summ.nis$cont$importance[2,"MDS1"]


# Plot sites and taxa together and colo by mobility
sites.NIS<-as.data.frame(scores(NIS.dbrda, display = "sites"))
sites.NIS$sample_name<-rownames(sites.NIS)
sites.NIS<-left_join(sites.NIS, metaNIS[,c("sample_name","Type_habitat_broad")], by ="sample_name")

dbrda.NIS<-ggplot() +
  geom_point(data=sites.NIS, 
             aes(x = CAP1, y = MDS1, color = Type_habitat_broad),size = 3) +
  # geom_point(data=cap.COMM.1,aes(x = CAP1, y = MDS1, color=Trait1),alpha=0.4,size = 1)+
  ggConvexHull:: geom_convexhull(data=sites.NIS,aes(x = CAP1, y = MDS1, color=Type_habitat_broad,
                                                     fill=Type_habitat_broad),alpha=0.0)+
  
  scale_color_manual(values = c( "blue","yellow2"))+ #,"pink","grey","red","green"+
  xlab(paste0("CAP 1 (", round(summ.nis$cont$importance[2,"CAP1"]*100, 1),"%)"))+
  ylab(paste0("PCoA 1 (", round(summ.nis$cont$importance[2,"MDS1"]*100, 1),"%)"))+
  theme_minimal()+
  theme(legend.position="top")+
  theme(axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17),
        legend.text = element_text(size=13),
        legend.title = element_text(size=14))


dbRDA.plots<-ggarrange(dbrda.NAT,dbrda.NIS,dbrda.UNC, nrow=1,common.legend = TRUE)

ggsave(plot=dbRDA.plots, "Final_files_pubblication/Figures/dbRDA_plots_NAT.NIS.UNC.pdf",
       dpi=300, width=16, height=6)

##### END 