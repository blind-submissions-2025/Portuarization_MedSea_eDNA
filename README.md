Files and codes used to generate the results reported in the manuscript "_Expanding on the portuarization syndrome from an ecological perspective: eDNA reveals rich diversity, non-indigenous hotspots, and biotic homogenization in ports_" 

**Diversity_analyses.R** : R script reporting all the statistical analyses performed in the study, as well as the script to generate plots reported in the main text and in the supplementary material

**Metadata.csv** : Metadata used in the study, sampling location, date and type of habitat (Port and Natural) for each sample.

**OTU_tables_three_markers.csv**: OTUs tables with full FASTA sequences generated for the three markers used in the study (12S, 18S and COI) before being combined into the final Species occurrence table. The abundances are not transformed into presence/absence at this stage. 

**Species_occurrence_table.csv**: Table including the 538 species retrieved by combining the trhee OTUs tables and keeping only OTUs assigned down to the Species level. Each species is characterised for its functionality (motility and habitat) and its status as NIS or Resident in the Mediterranean Sea.

**ps_538_species_summer2021.R** : phyloseq object used in the analyses, including taxonomic table, species table with presence and absence and metadata; it is possible to download this and run it in the "Diversity_analyses.R" script.
