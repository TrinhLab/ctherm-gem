# Files in this directory
## Models
- iCBI676.xml, corresponds to the model as downloaded from KBase in 10/3/2018
- The remaning models are generated from the notebook steps, which are executed in the order indicated by their file name. 
## Notebooks 
These are executed sequentially to obtain the final model.


# Consesus model
The consensus model, iCBI676.xml, was constructed from iSG676 with the following changes:
0. ACS is not blocked
1. Remove amino acid (ATP consuming) uptake reactions
2. Link trna cycling to biomass reactions
3. Remove the selenate pathway
4. Remove nad and nadp dependent hydrogenases. Keep ferredoxin hydrogenases with and without external protons
5. Add fumarate reductase
6. update all reactions based on ModelSEED metabolite charges

iCBI676.xml was downloaded from KBase


