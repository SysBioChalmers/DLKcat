function WriteFile_pre_cluster(a,b)

% init cluster function %% please modify this part based on the position of
% those tool box
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
workdir = pwd;
cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)

load('strains.mat')
inputpath = '/cephyr/users/feiranl/Hebbe/Yeast-Species-GEMs/Reconstruction_script/ModelFiles/mat';
species = strains;
species = species(a:b);
writeFileForKcatPrediction(species,inputpath)

end