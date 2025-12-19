Files and codes used to generate the results reported in the manuscript "_Expanding on the portuarization syndrome from an ecological perspective using environmental DNA_" 

**Preparation_combined_table.R**: R script used to: i) filter the OTU tables form each marker (12S,18S,COI) from OTUs not assigned to marine Metazoan, ii) rarefy to the minimum number of reads per sample, iii) combine the OTUs tables of the three markers (12S,18S,COI) into one non-redundant taxonomic dataset, iv) remove all OTUs not assigned down to species level.

**Diversity_analyses.R**: R script reporting all the statistical analyses performed in the study, as well as the script to generate plots reported in the main text and in the supplementary material

**Metadata.csv**: Metadata used in the study, sampling location, date and type of habitat (Port and Natural) for each sample.

**OTU_tables_three_markers.csv**: OTUs tables with full FASTA sequences generated for the three markers used in the study (12S, 18S and COI) before being combined into the taxonomic dataset. This table reports the abundances of each OTU in read counts. 

**Species_occurrence_table.csv**: Table including the 538 species retrieved by combining the three OTUs tables ("OTU_tables_three_markers.csv") and keeping only the OTUs assigned down to the Species level. In this table, the OTU abundances have been transformed into presence/absence. Each species reported in the table is characterised by its status as NIS,  Native or "Not-Identified-as-NIS"(NINIS in the table) - i.e. species not recorded as NIS in the Mediterranean, but not Native of this region. The last two classifications were combined in the study into the "RESS" group. The database used to assign the NIS status is reported in the ANISE_DB column. The functional traits classification (motility and ecology) of each species is also reported in this table. 

**ps_538_species_summer2021.R** : phyloseq object used in the analyses, including taxonomic table, species table with presence and absence and metadata; it is possible to download this and run it in the "Diversity_analyses.R" script.
