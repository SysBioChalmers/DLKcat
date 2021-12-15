%% updateDLkcat to in vitro measurement
function enzymedata = updatekcats(enzymedata,enzymedata_classic)
idx = find(enzymedata_classic.kcat_conf == 4);
enzyme_list = enzymedata_classic.enzyme(idx);
kcat_list = enzymedata_classic.kcat(idx);
rxn_list = enzymedata_classic.rxn_list(idx);
conf_list = enzymedata_classic.kcat_conf(idx);
for i = 1:length(enzyme_list)
    enzymename = enzyme_list(i);
    kcat_tmp = kcat_list(i);
    conf_tmp = conf_list(i);
    rxn_tmp = rxn_list(i);
    idx = ismember(enzymedata.enzyme,enzymename)&ismember(enzymedata.rxn_list,rxn_tmp);
    enzymedata.kcat(idx) = kcat_tmp;
    enzymedata.kcat_conf(idx) = conf_tmp;
end
