### Metatranscriptomic analysis of wastewater samples

- This repository contains the Nextflow script used to analyze wastewaster samples to idetify Antibiotic resistant genes.

#### Requirements 
- Nextflow
  - Nextflow installation - https://www.nextflow.io/docs/latest/getstarted.html
  -  Requirement (Java version 17+)   
- Docker or Conda
  - Docker or conda used to run each process in Nextflow.
  - Docker installation - https://docs.docker.com/engine/install/ubuntu/
  - Conda download and installation - https://www.anaconda.com/
 - Data files
   - Human genome retrieved from NCBI. Human host contaminating reads removed from the wastewater samples.
   - Ribosomal RNA reads from the SILVA database for multiple different species to leave behind messenger RNA for analysis.       

#### Steps 
We had a total of 24 samples for analysis collected across four sites.  RNA was extracted from individual samples, libraries were prepared and sequencing caried out using a Miseq (Ilumina).  We concatanated samples present in similar sites and ended up with four samples for analysis. 




  
