addpath(genpath('/proj/nobackup/snic2021-22-16/cplex1210/cplex/matlab'))
addpath(genpath('/home/f/feiranl/tools'))
addpath(genpath('/proj/nobackup/snic2021-22-16/DLKcat/'))
%addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/home/f/feiranl/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)
sched = parcluster('local');
sched.JobStorageLocation = getenv('TMPDIR');
parpool(sched, 18)

% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/libSBML-5.15.0-matlab'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/MATLAB-git'))
% addpath(genpath('/cephyr/users/feiranl/Hebbe/DLKcat'))
% workdir = pwd;
% cd '/cephyr/users/feiranl/Hebbe/tools/cobratoolbox'
% initCobraToolbox
% savepath '~/pathdef.m'
% cd(workdir)
% sched = parcluster('local');
% sched.JobStorageLocation = getenv('TMPDIR');
% parpool(sched, 18)