function writeClusterFile

%% get model_auto and model_dl
num = 343;
k = 1:10:num;
for i = 1:length(k)
    subfileName = ['sub',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    %fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:10
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
        end
        for m = 1:10
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "GrowthAnalysis_cluster($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    else
        for m = 1:num-k(i)+1
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(1*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(1*m-1) '\n']);
        end
        for m = 1:num-k(i)+1
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "GrowthAnalysis_cluster($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end


%% get model_bayesian
num = 343;
k = 1:1:num;
for i = 171:length(k)
    subfileName = ['sub',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-10:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    fprintf(fptr,'b=$i\n');

    fprintf(fptr,'matlab -nodesktop -singleCompThread -r "BayesianModelGeneration_cluster($i,$b,false)"');
    
    fclose(fptr);
end