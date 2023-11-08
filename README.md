# Code in support of Article in Conservation Biology

**Title: Effects of protection and temperature variation on temporal stability in a marine reserve network**

Srednick & Swearer 2023

doi: [10.1111/cobi.14220](https://conbio.onlinelibrary.wiley.com/doi/10.1111/cobi.14220)

## Included scripts

- 0_code_pathway_sourced.R --> runs all following scripts via source
- 1_data_curation.R --> curates and formats PISCO subtidal survey data
- 2_distance_method.R --> calculates over-water distance between sites; prep
- 3_geo_drivers_synchrony.R --> calculates SST metrics for each community; metacommunity
- 4_synchrony_stability.R --> calculates synchrony and stability for each community; metacommunity
- 5_species_interactions.R --> uses MARSS modelling to estimate interaction strength at community; metacommunity scale
- 6_hab_similarity.R --> calculates similarity indices for abiotic habitat structure at each community; metacommunity.
- 7_sem_network.R --> casual structural equation model (SEM) analyses
- 8_figures.R --> main figure generator
- 9_SEM_path_plots.R --> SEM figure generator


![script_flowchart](./Manuscript_scripts/script_flowchart.pdf)
