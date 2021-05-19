function autoDLModelGeneration_cluster(a,b)
addpath(genpath('/proj/nobackup/snic2021-22-16/cplex1210/cplex/matlab'))
addpath(genpath('/home/f/feiranl/tools'))
addpath(genpath('/proj/nobackup/snic2021-22-16/MLKcat/'))
%addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/home/f/feiranl/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

modelpath = '/proj/nobackup/snic2021-22-16/models';
dbpath = '/proj/nobackup/snic2021-22-16/newpep';
kcatpredictionPath = '/proj/nobackup/snic2021-22-16/MLKcat/ComplementaryScripts/ForKcatPrediction/KcatPredictionResult';
load('StrianData.mat')
species = StrianData.strains;
species = species(a:b);
autoDLModelGeneration(species,modelpath,dbpath,kcatpredictionPath);
end