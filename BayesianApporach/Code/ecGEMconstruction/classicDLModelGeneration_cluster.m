function classicDLModelGeneration_cluster(a,b)

initcluster
modelpath = '../../Results/ssGEMs';
dbpath = '../../Results/Proteinfasta';
kcatpredictionPath = '../../Results/PredcitedKcat343species';
load('StrianData.mat')
species = StrianData.strains;
species = species(a:b);
classicDLModelGeneration(species,modelpath,dbpath,kcatpredictionPath);
end