### Metatranscriptomic analysis of wastewater samples

- This repository contains the Nextflow script used to analyze wastewaster samples to idetify Antibiotic resistant genes.

#### Requirements 
- Nextflow
  - Nextflow installation - https://www.nextflow.io/docs/latest/getstarted.html
  -  Requirement (Java version 17+)   
- Docker or Conda
  - Docker or conda used to run each process in Nextflow.
  - Docker installation (https://docs.docker.com/engine/install/ubuntu/)
 - Data files
   - We retrieved the human genome from NCBI. Wastewater samples contain human reads that need to be remove host contaminating reads
   - We retrieved ribosomal RNA reads from the SILVA database      

#### Steps 
We had a total of 24 samples for analysis collected across four sites.  RNA was extracted from individual samples, libraries were prepared and sequencing caried out using a Miseq (Ilumina). 


  
