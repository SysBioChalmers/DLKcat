% this function is to sort MNX prop and get the mapping info for all MNX
% IDs and deprecated IDs

missingID = {'Starch;MNXM725906'
'galactitol;MNXM1233'
'L-xylo-3-hexulose;MNXM8830'
'inulin;MNXM737727'
'D-erythrulose;MNXM3572'
'D-erythrulose 4-phosphate;MNXM48433'
'erythritol;MNXM734797'
'erythritol;MNXM734797'
'D-gluconate;MNXM341'
'nitrate;MNXM732398'
'Nitrite;MNXM107'
'propane-1,2-diol;MNXM1118'
'propane-1,2-diol;MNXM1118'
'beta-cellobiose;MNXM726890'
'beta-cellobiose;MNXM726890'
'lactose;MNXM734836'
'lactose;MNXM734836'
'2-dehydro-D-gluconate;MNXM582'
'2-dehydro-D-gluconate;MNXM582'
'D-ribulose;MNXM735687'
'D-glucose 6-phosphate;MNXM729991'
'D-arabinose;MNXM544'
'D-fructose;MNXM729613'
'D-galactose;MNXM735264'
'glycogen;MNXM738130'
'D-mannose;MNXM731377'
'D-mannose 6-phosphate;MNXM735948'
% 'L-sorbose;MNXM729623' this can be found in pubchem
'D-ribose;MNXM732917'};
missingID = split(missingID,';');

fid2 = fopen('chem_prop.tsv'); % downloaded from https://www.metanetx.org/mnxdoc/mnxref.html
format = repmat('%s ',1,9);
format = strtrim(format);
met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXprop(:,i) = met_temp{i};
end
commentLines = startsWith(MNXprop(:,1),'#');
MNXprop(commentLines,:) = []; % MNX_ID name REFERENCE formula charge mass standardInChI standardInChIKey SMILES 
fclose(fid2);


fid2 = fopen('chem_depr.tsv'); % downloaded from https://www.metanetx.org/mnxdoc/mnxref.html
format = repmat('%s ',1,3);
format = strtrim(format);
met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXdepr(:,i) = met_temp{i};
end
commentLines = startsWith(MNXdepr(:,1),'#');
MNXdepr(commentLines,:) = []; % MNX_ID_old MNX_ID_new
fclose(fid2);

fid2 = fopen('MNX_met_yeast8.3.4.txt'); % downloaded from https://www.metanetx.org/cgi-bin/mnxget/user_model/yeastnet_v8_3_4.chemicals.tsv
format = repmat('%s ',1,8);
format = strtrim(format);
met_temp = textscan(fid2,format,'Delimiter','\t','HeaderLines',0);
for i = 1:length(met_temp)
    MNXYeast(:,i) = met_temp{i};
end
commentLines = startsWith(MNXYeast(:,1),'#');
MNXYeast(commentLines,:) = []; % MNX_ID name REFERENCE formula charge mass standardInChI standardInChIKey SMILES 
fclose(fid2);

MNXYeastList = [];
for i = 1:length(MNXYeast(:,1))
    ID_tmp = split(MNXYeast(i,3),';');
    MNXYeastList = [MNXYeastList;repmat(MNXYeast(i,1),length(ID_tmp),1),ID_tmp];
end
MNXYeast = MNXYeastList;
MNXYeast(startsWith(MNXYeast(:,1),'UNK:'),:) = [];
clearvars -except MNXprop MNXdepr MNXYeast missingID
save('MNXprop.mat')