function classicDLModelGeneration_cluster(a,b)

initcluster
modelpath = '../../Results/model_build_files/splitedmodel';
dbpath = '../../Results/Proteinfasta';
kcatpredictionPath = '../../Results/PredcitedKcat343species';
load('strains.mat')
species = strains(a:b);
classicDLModelGeneration(species,modelpath,dbpath,kcatpredictionPath);
end