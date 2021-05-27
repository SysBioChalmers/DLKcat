function PCAsampledKcatsOneSpecies(species)

current_path = pwd;
cd ../../Results/model_bayesian/
cd(species)
nfound = length(dir('kcat_genra*.txt'));
for i = 1:nfound
    display([num2str(i) '/' num2str(nfound)]);
    tmp = readmatrix(['kcat_genra',num2str(i),'.txt'],'FileType','text','Delimiter',',');
    theta_100 = tmp(end,:); % is the rmse error
    kcat_100 = tmp(1:end-2,:);
    all(i*100-99:i*100,:) = kcat_100';
    theta_all(i*100-99:i*100,1) = theta_100';
end

% all(theta_all > 10,:) = [];
% theta_all(theta_all > 10) = [];

    mean_tmp = mean(all,1);
    std_tmp = std(all,1);
    all_norm = (all - mean_tmp) ./ std_tmp;
    
    [~, score, ~, ~, explained, ~] = pca(all_norm,'NumComponents',2);
  save(['res_ForKcatPCA',species,'.mat'],'score','explained','theta_all')
end