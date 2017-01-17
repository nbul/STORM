clc
clear variables
close all


%% Defining extension
% Default cutoff
cutoff =0.80;
Column_intensity = 18;
Ring = struct([]);
coordinates = struct([]);
Cluster = struct([]);
Centroid = struct([]);
coordinates2 = struct([]);
coordinates3 = struct([]);
IntensityCluster = struct([]);


usedefault = questdlg(strcat('Use default settings: (Intensity column = ',...
    num2str(Column_intensity),'?)'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No');
    parameters = inputdlg({'Enter intensity column:'},'Parameters',1, {num2str(Column_intensity)});
    % Redefine extension
    Column_intensity = str2double(parameters{1});
else
    parameters{1} = num2str(Column_intensity);
end

%% Determening paths
currdir = pwd;
addpath(pwd);
filedir = uigetdir();
files = dir(strcat(filedir,'/*.csv'));
cd(filedir);
mkdir(filedir,'/results');
resultdir = [filedir, '/results'];
mkdir(resultdir,'/Clustering');
clustdir = [resultdir,'/Clustering'];
mkdir(resultdir,'/Lines');
linedir = [resultdir,'/Lines'];


Index = zeros(numel(files),1);
Distance = zeros(numel(files),1);

%% Reading files and all coordinates
for i=1:numel(files)
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4)];
end

Sigma = 'full'; % possible to change to 'diagonal'
nSigma = numel(Sigma);
SharedCovariance = false; % possible to change to 'true'
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000);

%% Clustering
for l=1:numel(files)
    
    k=2;
    d = 500;
    x1 = linspace(min(coordinates{l}(:,1)) - 2,max(coordinates{l}(:,1)) + 2,d);
    x2 = linspace(min(coordinates{l}(:,2)) - 2,max(coordinates{l}(:,2)) + 2,d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];
    i = 0;
    while i == 0;
        image1 = figure;
        gmfit = fitgmdist(coordinates{l},k,'CovarianceType',Sigma,...
            'SharedCovariance',SharedCovariance,'Options',options);
        clusterX = cluster(gmfit,coordinates{l});
        mahalDist = mahal(gmfit,X0);
        coordinates2{l} = [coordinates{l} clusterX];
        h1 = gscatter(coordinates{l}(:,1),coordinates{l}(:,2),clusterX);
        hold on;
        for m = 1:k;
            idx = mahalDist(:,m)<=threshold;
            Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
            h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
            uistack(h2,'bottom');
        end
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10)
        
        usedefault = questdlg(strcat('Are you happy about clustering)'),'Settings','Yes','No','Yes');
        if strcmp(usedefault, 'Yes');
            i=1;
            cd(clustdir);
            print(image1, [num2str(l),'_cluster.tif'], '-dtiff', '-r150');
            cd(currdir);
        end
        close all;
    end
end

%% Identifying a brighter cluster for fit
for i=1:numel(files)
    IntensityCluster{i}(1) = sum(Ring{i}(coordinates2{i}(:,3)==1,Column_intensity));
    IntensityCluster{i}(2) = sum(Ring{i}(coordinates2{i}(:,3)==2,Column_intensity));
    if IntensityCluster{i}(1)>IntensityCluster{i}(2)
        Index(i) = 1;
    else
        Index(i) = 2;
    end
end

%% Splitting proportianally to the signal
for i=1:numel(files)
    
    k = 1;
    for o=1:length(Ring{i})
        m=ceil(Ring{i}(o,18)/100);
        for f=1:m
            coordinates3{i}(k,:) = coordinates2{i}(o,:);
            k=k+1;
        end
    end
    
end

%% Fitting clusters with lines and finding distance
for i=1:numel(files)
    Line1 = polyfit(coordinates3{i}(coordinates3{i}(:,3)==Index(i),1),coordinates3{i}(coordinates3{i}(:,3)==Index(i),2), 1);
    image2 = figure;
    plot(coordinates2{i}(coordinates2{i}(:,3)==1,1),coordinates2{i}(coordinates2{i}(:,3)==1,2), 'r*');
    hold on
    plot(coordinates2{i}(coordinates2{i}(:,3)==2,1),coordinates2{i}(coordinates2{i}(:,3)==2,2), 'g*');
    hold on
    x1 = min(coordinates2{i}(coordinates2{i}(:,3)==Index(i),1)):0.005:max(coordinates2{i}(coordinates2{i}(:,3)==Index(i),1));
    y1 = Line1(1,1)*x1 + Line1(1,2);
    plot(x1, y1 , '-b', 'LineWidth',3);
    hold on;
    ft = fittype('a*x + b','coefficients','b','independent','x','problem','a');
    Line2 = fit(coordinates3{i}(coordinates3{i}(:,3)~=Index(i),1),coordinates3{i}(coordinates3{i}(:,3)~=Index(i),2), ft, 'problem',Line1(1,1));
    x2 = min(coordinates2{i}(coordinates2{i}(:,3)~=Index(i),1)):0.005:max(coordinates2{i}(coordinates2{i}(:,3)~=Index(i),1));
    y2 = Line2.a*x2 + Line2.b;
    plot(x2, y2 , '-b', 'LineWidth',3);
    cd(linedir);
    print(image2, [num2str(i),'_lines.tif'], '-dtiff', '-r150');
    cd(currdir); 
    close all;
    
    Distance(i) = abs(Line1(1,2)-Line2.b)/sqrt(Line2.a*Line2.a + 1);
    
end
Filenames = 1:1:numel(files);
Distance2 = [Filenames' Distance];
cd(resultdir);
Otput_summary = 'Distances.csv';
headers2 = {'Ring', 'Distance'};
csvwrite_with_headers(Otput_summary,Distance2, headers2);
cd(currdir); 

cd(pwd);