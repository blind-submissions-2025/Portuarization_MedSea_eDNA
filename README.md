Files and codes used to generate the results reported in the manuscript "_An ecological perspective on the portuarization syndrome: eDNA reveals species-rich communities, non-indigenous species hotspots, and biotic homogenization in ports_" 

**Preparation_combined_table.R**: R script used to: i) filter the OTU tables form each marker (12S,18S,COI) from OTUs not assigned to marine Metazoan, ii) rarefy to the minimum number of reads per sample, iii) combine the OTUs tables of the three markers (12S,18S,COI) into one non-redundant taxonomic dataset, iv) remove all OTUs not assigned down to species level.

**Diversity_analyses.R**: R script reporting all the statistical analyses performed in the study, as well as the script to generate plots reported in the main text and in the supplementary material

**Metadata.csv**: Metadata used in the study, sampling location, date and type of habitat (Port and Natural) for each sample.

**OTU_tables_three_markers.csv**: OTUs tables with full FASTA sequences generated for the three markers used in the study (12S, 18S and COI) before being combined into the taxonomic dataset. This table reports the abundances of each OTU in read counts. 

**ps_536_species.R** : phyloseq object used in the analyses, including taxonomic table, species table with presence and absence and metadata; it is possible to download this and run it in the "Diversity_analyses.R" script. The taxonomy table and the species occurrence tables were used to generate the Species_occurrence_tables reported in the Appendix S.2 Table S1.
