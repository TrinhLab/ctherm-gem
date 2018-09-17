# Consesus model
The consensus model, iCBI676, was constructed from iSG676 with the following changes:
- ACS is not blocked
1. Remove amino acid (ATP consuming) uptake reactions
2. Link trna cycling to biomass reactions
3. Remove the selenate pathway
4. Remove nad and nadp dependent hydrogenases. Keep ferredoxin hydrogenases with and without external protons
5. Add fumarate reductase
6. update all reactions based on ModelSEED metabolite charges

It was downloaded from KBase

# Compound mapping
- compund.tsv is a reduced version (less columns) of the modelSEED/KBase compound database. Used to map the modelSEED ids to BiGG IDs
- bigg2ms.csv was is a manually extended version of the output of create_id_map.ipynb (bigg2ms_raw.csv)
