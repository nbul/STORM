clc 
clear variables 

[fname, pname] = uigetfile('*.csv');
FileName = fullfile(pname, fname);
Data = csvread(FileName);

GMModels = cell(2,1); % Preallocation
options = statset('MaxIter',1000);
rng(1); % For reproducibility

Data(:,(Data(1,:) == 0)) = [];
numComponents = zeros(1, length(Data(1,:)));
Model = cell(length(Data(1,:)),1); % Preallocation
Clustering = zeros(2,length(Data(1,:)));
CleanData = zeros(2,length(Data(1,:)));

for i = 1:length(Data(1,:))
    Temp = Data(:,i);
    Temp(Temp == 0) = [];
    
    AIC = zeros(2,1);
    for j = 1:2
        GMModels{j} = fitgmdist(Temp,j,'Options',options, 'CovarianceType','diagonal','RegularizationValue',0.01);
        AIC(j)= GMModels{j}.BIC;
    end
    [minAIC,numComponents(i)] = min(AIC);
    
    Model{i} = GMModels{numComponents(i)};
    Clustering(1:length(Temp),i) = cluster(Model{i}, Temp);
    
    if numComponents(i) == 2
        if Model{i}.mu(1) < Model{i}.mu(2)
            Temp(Clustering(1:length(Temp),i) == 2) = [];
        else
            Temp(Clustering(1:length(Temp),i) == 1) = [];
        end
    end
    CleanData(1:length(Temp),i) = Temp;
    
end


