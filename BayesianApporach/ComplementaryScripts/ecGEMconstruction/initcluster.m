addpath(genpath('/proj/nobackup/snic2021-22-16/cplex1210/cplex/matlab'))
addpath(genpath('/home/f/feiranl/tools'))
addpath(genpath('/proj/nobackup/snic2021-22-16/MLKcat/'))
%addpath(genpath('/cephyr/users/feiranl/Hebbe/tools/RAVEN'))
workdir = pwd;
cd '/home/f/feiranl/tools/cobratoolbox'
initCobraToolbox
savepath '~/pathdef.m'
cd(workdir)
sched = parcluster('local');
sched.JobStorageLocation = getenv('TMPDIR');
parpool(sched, 18)