% This scripts validates the configuration of the matlab model by ensuring
% tha that hydG-ech-pfl growth prediction is correct. 

lin =  load('/home/sergio/Dropbox/cthermgem-dev/iSG676/iSG676_cb.mat');
model = lin.iSG;

s = optimizeCbModel(model);
wt_gr = s.f;

delete_id = {'BIF','H2ASE_syn','ECH','PFL'};
model_ko = changeRxnBounds(model, delete_id, 0, 'b');
s = optimizeCbModel(model_ko);
ko_gr = s.f;
max_expected_growth_frac = 0.21;
assert(ko_gr/wt_gr <= max_expected_growth_frac)
fprintf('Fraction of hydg-ech-pfl growth %f\n', ko_gr/wt_gr)