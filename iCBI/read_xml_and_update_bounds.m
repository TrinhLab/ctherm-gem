%%This code uses functions from the COBRA toolbox
%% Required input files:1. iCBI676.xml 2.fba_results.xls 
model=readCbModel('iCBI676.xml');

%change reaction bounds based on the bounds defined in the flux file which
%was generated using KBase

rxn_bnds = xlsread('fba_result.xls','ModelReactions');
met_bnds = xlsread('fba_result.xls','ModelCompounds');

model.lb(1:799)=rxn_bnds(:,3); % update lower bounds
model.ub(1:799)=1000* ones(799,1); % update upper bounds
model.lb(800:853)=-met_bnds(:,7); % update uptake lower bound
model.ub(800:853)=-met_bnds(:,5); % update uptake upper bound

% update objective function to Biomass_cellobiose
model.c(799)=0;
model.c(797)=1;

writeCbModel(model,'sbml','iCBI665_updated.xml');
