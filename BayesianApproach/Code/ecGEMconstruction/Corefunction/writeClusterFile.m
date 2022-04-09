function writeClusterFile

%% get model_auto and model_dl
num = 343;
k = 1:100:num;
for i = 1:length(k)
    subfileName = ['sub',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    %fprintf(fptr,'#SBATCH -A C3SE2021-1-16\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-5:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    %fprintf(fptr,'module load MATLAB intel/2018b GMP\n');
        fprintf(fptr,'module load MATLAB/2021a GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    if i ~= length(k)
        for m = 1:5
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(20*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(20*m-1) '\n']);
        end
        for m = 1:5
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "classicDLModelGeneration_cluster($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
    else
        for m = 1:4
            fprintf(fptr,['a',num2str(m),'=$i+' num2str(10*(m-1)) '\n']);
            fprintf(fptr,['b',num2str(m),'=$i+' num2str(10*m-1) '\n']);
        end
        fprintf(fptr,['a',num2str(5),'=$i+' num2str(10*4) '\n']);
        fprintf(fptr,['b',num2str(5),'=' num2str(num) '\n']);
        for m = 1:4
            fprintf(fptr,['matlab -nodesktop -singleCompThread -r "classicDLModelGeneration_cluster($a',num2str(m),',$b',num2str(m),')" &\n']);
        end
        
       fprintf(fptr,['matlab -nodesktop -singleCompThread -r "classicDLModelGeneration_cluster($a',num2str(5),',$b',num2str(5),')" &\n']);  
    end
    fprintf(fptr,'wait;\n');
    fclose(fptr);
end


%% get model_bayesian
num = 343;
k = 171:1:num;
for i = 1:length(k)
    subfileName = ['sub_k',num2str(k(i)),'.sh'];
    fptr = fopen(subfileName,'w');
    fprintf(fptr,'#!/bin/bash\n');
    %     fprintf(fptr,'#SBATCH -A SNIC2020-7-19\n');
    fprintf(fptr,'#SBATCH -A snic2021-22-21\n');
    fprintf(fptr,'#SBATCH -n 20\n');
    fprintf(fptr,'#SBATCH -o out.txt\n');
    fprintf(fptr,'#SBATCH -t 0-10:00:00\n');
    fprintf(fptr,'#SBATCH --mail-user=feiranl@chalmers.se\n');
    fprintf(fptr,'#SBATCH --mail-type=end\n');
    fprintf(fptr,'module load MATLAB/2021a GMP\n');
    fprintf(fptr,['i=',num2str(k(i)),'\n']);
    fprintf(fptr,'b=$i\n');

    fprintf(fptr,'matlab -nodesktop -singleCompThread -r "BayesianModelGeneration_cluster($i,$b,false)"');
    
    fclose(fptr);
end