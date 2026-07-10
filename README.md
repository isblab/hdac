
[Add PubMed link]: [![PubMed]

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11056108.svg)](https://zenodo.org/doi/10.5281/zenodo.11056108)

# Multi-complex Integrative Structure Determination of the HDAC1/2 Interactome

This repository is for the integrative models of the HDAC1 corepressor complexes - NuRD, Sin3A, coREST, based on  data from chemical crosslinking, cryo-EM maps, X-ray crystallography, homology modeling using [Modeller](https://salilab.org/modeller/), and structure prediction from [Alphafold](https://www.alphafold.ebi.ac.uk/). It contains input data, scripts for data preprocessing, modeling and results including bead models and localization probability density maps. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*). 

These integrative structures are deposited in the [PDB-IHM](https://pdb-ihm.org/) . The accession codes are 9A8O: NuRD, 9A8P: CoREST, PA8Q: Sin3A. 

![main_fig](main_figure.png)


## Directory structure
1. [input](input/) : contains the subdirectories for the input data used for modeling all the corepressor complexes.
2. [scripts](scripts/) : contains all the scripts used for pre-processing, modeling and analysis of the models.
3. [results](results/) : contains the models and the localization probability densities of the top cluster of the corepressor complexes.
4. [test](test/) : scripts for testing the sampling.


## Protocol
### Preprocessing 
1. In case of structures predicted by AlphaFold2, only regions of high confidence (>70 pLDDT and <5 PAE) were used. Following scripts extracts regions of high confidence:  
    ```
    python get_high_confidence_region_from_AF2.py af2_struct.cif af2_struct.json
    ```

2. For the presence of multiple paralogs of a protein, XLs from all paralogs were mapped to the paralog with the highest number of XLs.  
   Following script generates the mapped XLs:
    ```
    python paralog_alignment.py
    ```

3. EM maps, where available, were converted to Gaussian Mixture Models (GMM) which were used as input for the modeling.
   Following script generates GMMs for the input EM map:
   ```
   ./create_gmm.sh EM_map.mrc threshold
   ```
   Threshold can be obtained from the Validation section on [EMDB](https://www.ebi.ac.uk/emdb/) for the specific EM map.
   The minimum number of Gaussians which yield a cross-correlation of >0.95 with the original EM map is used.


### Sampling
To run the sampling, run modeling scripts like this   
```
./run_modeling.sh
```


### Analysis
To run the analysis, run the end-t-end analysis script like   
```
python end_to_end_analysis.py
```

The above script does the following -
#### 1. Getting the good-scoring models
  Good-scoring models were selected using `pmi_analysis` (Please refer to [pmi_analysis tutorial](https://github.com/salilab/PMI_analysis) for more detailed explaination) along with our `variable_filter_v1.py` script.  
  
  Following are the scripts used:
  1. `run_analysis_trajectories.py`  
      
  2. `variable_filter_v1.py` on the major cluster if the number of models exceeds 30000.   

  3. The selected good-scoring models were then extracted using `run_extract_models.py`.:   
      
#### 2. Running the sampling exhaustiveness tests (Sampcon)
A `density_{}.txt` (Nurd, corest, Sin3a) file was created. This file contains the details of the domains to be split for visualizing the localisation probability densities. Finally, sampling exhaustiveness tests were performed using `imp-sampcon`. 

#### 3. Analysing the major cluster
1. Compute crosslink violations using `get_xl_viol_validation_set.py` script.   
      
2. Create contact maps for the component proteins in the complex using `contact_maps_all_pairs_surface.py` script. The proteins to be considered are specified as lists `protein1` and `protein2`.

3. Obtain domainwise precision using [PrISM](https://doi.org/10.1093/bioinformatics/btac400).

### Results

For the simulations, the [results](results/) directory consists of a subdirectory for each complex comprising of -
* `contact_maps` : Directory containing of the contact map for the component proteins of the complex.
* `models_and_densities` : Directory containing sampcon output for the largest cluster.
* `prism` : Directory containing the PrISM output.
* `xl_violations` : Directory containing the logs for crosslink violations.


## Information
**Author(s):** Jules Nde*, Kartik Majila*, Rosalyn C. Zimmermann, Cassandra Kempf, Ying Zhang, Joseph Cesare, Janet L. Thornton, Jerry L. Workman, Laurence Florens, Shruthi Viswanath, Michael P. Washburn  
**Date**: September 12th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** Not tested  
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:** Multicomplex Integrative Structural Modeling of a Human Histone Deacetylase Interactome, in revisions 

