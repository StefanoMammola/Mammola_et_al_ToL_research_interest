Data and R script to reproduce the analyses in:

Mammola S., et al. (2023) Drivers of species knowledge across the Tree of Life. eLife
 
- Folder "Scripts":

01.GetRandomSpecies: Extract random species from GBIF backbone
02.GetGBIFData.R: Extract data from GBIF
03.GetIUCNData.R: Extract data from IUCN
04.GetWoSdata.R: Extract data from the Web of Science
05.Analysis&Figures.R Reproduce analyses and figures

[Note that script 01 - 04 are provided to illustrate the procedure of sampling data. All analyses and figures can be ran simply using the "05.Analysis&Figures.R" and associated data]

- Folder "Functions":

Functions.R: this R file is sourced during analyses providing custom fuctions and plot parameters.

- Folder "Data":

SampleTREE.csv + SampleTREE.xlsx: Initially sampled species
SampleTREE_TraitCompiled.csv: Sample species + compiled traits (this is the database needed to reproduce analyses)
sel_rareas.rds: needed to run script "04.GetWoSdata.R"

- Folder "Phylopics":

Provides the .png of animal, plant, and fungi silhouettes used in figures. Silhouettes were taken from PhyloPics (http://phylopic.org/) - all with open licence.

- Folder "Figures" and "Tables":

Stores figures and tables generated during analyses.