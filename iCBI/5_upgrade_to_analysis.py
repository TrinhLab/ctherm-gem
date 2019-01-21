"""
This script generates the final version of the model and does two things:

1. Update the model according to the essentiality analysis findings (see /analysis/essentiality/known_essentiality.ipynb)
2.
"""
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import cobra as cb
import settings


# 1
model = cb.io.load_json_model(os.path.join(settings.PROJECT_ROOT,'iCBI','iCBI655bigg_cellb_batch.json'))

model.reactions.DRPA.bounds = (0,0)
model.reactions.DRPA.notes['curation_notes'] = 'The bounds are set to zero to correctly predict hydg-ech-pfl lethal mutant '


# 2
isg = cb.io.load_json_model(os.path.join(settings.PROJECT_ROOT,'iSG676','iSG676_cb.json'))
for met in model.metabolites:
    if met in isg.metabolites:
        met.formula = isg.metabolites.get_by_id(met.id).formula
cb.io.save_json_model(model, os.path.join(settings.PROJECT_ROOT,'iCBI','iCBI655bigg_cellb_batch_v2.json'))
cb.io.write_sbml_model(model, os.path.join(settings.PROJECT_ROOT,'iCBI','iCBI655bigg_cellb_batch_v2.sbml'))
cb.io.save_matlab_model(model, os.path.join(settings.PROJECT_ROOT,'iCBI','iCBI655bigg_cellb_batch_v2.mat'))
