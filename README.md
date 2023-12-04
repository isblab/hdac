
[Add PubMed link]: [![PubMed]

[Add Zenodo link]: [![DOI]

# Multi-complex Integrative Structure Determination of the HDAC1/2 Interactome

This repository is for the integrative models of the HDAC1 corepressor complexes - NuRD, Sin3A, coREST, based on  data from chemical crosslinking, cryoEM maps, X-ray crystallography, homology modeling using [Modeller](https://salilab.org/modeller/), and structure prediction from [Alphafold](https://www.alphafold.ebi.ac.uk/). It contains input data, scripts for data preprocessing, modeling and results including bead models and localization probability density maps. The modeling was performed using [IMP](https://integrativemodeling.org) (*Integrative Modeling Platform*). 

These integrative structures will be deposited in the PDB-Dev database with accession codes ***AddPDBdev Ids***

![main_fig](F1.png)


## Directory structure
1. [input](input/) : contains the subdirectories for the input data used for modeling all the corepressor complexes.
2. [scripts](scripts/) : contains all the scripts used for pre-processing, modeling and analysis of the models.
3. [results](results/) : contains the models and the localization probability densities of the top cluster of the corepressor complexes.
4. [test](test/) : scripts for testing the sampling.


## Protocol
### Preprocessing the crosslinks
1. For crosslinks in sheetA of sheetA of `data/xlinks/original_suppmat_DataS3.xlsx`:  
    ```
    python get_protein_uniprot_mapping.py -x /home/shreyas/Dropbox/washburn_wdr_spin/xls_sheet1.csv
    ```
    Make a file `proteins_of_interest.txt` and use it to run:
    ```
    python get_protein_uniprot_mapping.py -x /home/shreyas/Dropbox/washburn_wdr_spin/xls_sheet1.csv -m mapping -p proteins_of_interest.txt
    ```
    Finally, to generate the input file for modeling, do:
    ```
    python xl_preprocessing.py ~/Dropbox/washburn_wdr_spin/xls_sheet1.csv uniprot_mapping.yaml
    ```
2. For crosslinks in sheetA of sheetA of `data/xlinks/original_suppmat_DataS3.xlsx`:  
    Run the following command to generate the crosslinks input file for modeling:
    ```
    python xl_change_protnames_from_ncbi2name.py
    ```


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
#### 1. Getting the good scoring models
  Good scoring models were selected using `pmi_analysis` (Please refer to [pmi_analysis tutorial](https://github.com/salilab/PMI_analysis) for more detailed explaination) along with our `variable_filter_v1.py` script.  
  Following are the scripts used:
  1. `run_analysis_trajectories.py`  
      
  2. `variable_filter_v1.py` on the major cluster if the no. of models exceed 30000.   

  3. The selected good scoring models were then extracted using `run_extract_models.py`.:   
      
#### 2. Running the sampling exhaustiveness tests (Sampcon)
A `density_{}.txt` (Nurd, corest, Sin3a) file was created. This file contains the details of the domains to be split for plotting the localisation probability densities. Finally, sampling exhaustiveness tests were performed using `imp-sampcon`. 

#### 3. Analysing the major cluster
1. Compute crosslink violations using `get_xl_viol_validation_set.py` script.   
      
2. Create contact maps for the component proteins in the complex using `contact_maps_all_pairs_surface.py` script. The proteins to be considered are specified as lists `protein1` and `protein2`.

3. Obtain domainwise precision using PrISM [PrISM](https://doi.org/10.1093/bioinformatics/btac400).


### Results

For the simulations, the [results](results/) directory consists of a subdirectory for each complex comprising of -
* `Contact_maps` : Directory containing of the contact map for the component proteins of the complex.
* `cluster.0` : Directory containing sampcon output for the largest cluster.
* `prism` : Directory containing the PrISM output.
* `xl_violations` : Directory containing the logs for crosslink violations.


## Information
**Author(s):** Rosalyn C. Zimmermann*, Kartik Majila*, Cassandra Kempf*, Charles A.S. Banks, Mark K. Adams, Sayem Miah, Janet L. Thornton, Mihaela Sardiu, Laurence Florens, Shruthi Viswanath, Michael P. Washburn  
**Date**: September 12th, 2023  
**License:** [CC BY-SA 4.0](https://creativecommons.org/licenses/by-sa/4.0/)
This work is licensed under the Creative Commons Attribution-ShareAlike 4.0
International License.  
**Last known good IMP version:** Not tested  
**Testable:** Yes  
**Parallelizeable:** Yes  
**Publications:**  

