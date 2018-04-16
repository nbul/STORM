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
if strcmp(usedefault, 'No')
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

if exist([filedir,'/results'],'dir') == 0
    mkdir(filedir,'/results');
end
resultdir = [filedir, '/results'];

if exist([filedir,'/Clustering'],'dir') == 0
    mkdir(filedir,'/Clustering');
end
clustdir = [resultdir,'/Clustering'];

if exist([filedir,'/Lines'],'dir') == 0
    mkdir(filedir,'/Lines');
end
linedir = [resultdir,'/Lines'];

if exist([filedir,'/Distributions'],'dir') == 0
    mkdir(filedir,'/Distributions');
end
distdir = [resultdir,'/Distributions'];



Index = zeros(numel(files),1);
Distance = zeros(numel(files),1);
Line1 = zeros(numel(files),2);
Line2 = zeros(numel(files),2);

%% Reading files and all coordinates
for i=1:numel(files)
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4) Ring{i}(:,Column_intensity)];
end

Sigma = 'full'; % possible to change to 'diagonal'
nSigma = numel(Sigma);
SharedCovariance = false; % possible to change to 'true'
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000);
k=2;
d = 500;
%% Clustering
for l=1:numel(files)
    x1 = linspace(min(coordinates{l}(:,1)) - 2,max(coordinates{l}(:,1)) + 2,d);
    x2 = linspace(min(coordinates{l}(:,2)) - 2,max(coordinates{l}(:,2)) + 2,d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];
    i = 0;
    while i == 0;
        image1 = figure;
        gmfit = fitgmdist(coordinates{l}(:,1:2),k,'CovarianceType',Sigma,...
            'SharedCovariance',SharedCovariance,'Options',options);
        clusterX = cluster(gmfit,coordinates{l}(:,1:2));
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
    IntensityCluster{i}(1) = sum(coordinates2{i}(coordinates2{i}(:,4)==1,3));
    IntensityCluster{i}(2) = sum(coordinates2{i}(coordinates2{i}(:,4)==2,3));
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
    Line1(i,:) = polyfit(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1),coordinates3{i}(coordinates3{i}(:,4)==Index(i),2), 1);
    image2 = figure;
    plot(coordinates2{i}(coordinates2{i}(:,4)==1,1),coordinates2{i}(coordinates2{i}(:,4)==1,2), 'r*');
    hold on
    plot(coordinates2{i}(coordinates2{i}(:,4)==2,1),coordinates2{i}(coordinates2{i}(:,4)==2,2), 'g*');
    hold on
    x1 = min(coordinates2{i}(coordinates2{i}(:,4)==Index(i),1)):0.005:max(coordinates2{i}(coordinates2{i}(:,4)==Index(i),1));
    y1 = Line1(i,1)*x1 + Line1(i,2);
    plot(x1, y1 , '-b', 'LineWidth',3);
    hold on;
    ft = fittype('a*x + b','coefficients','b','independent','x','problem','a');
    Line2fit = fit(coordinates3{i}(coordinates3{i}(:,4)~=Index(i),1),coordinates3{i}(coordinates3{i}(:,4)~=Index(i),2),...
        ft, 'problem',Line1(i,1),'StartPoint', Line1(i,2));
    Line2(i,1) = Line1(i,1);
    Line2(i,2) = coeffvalues(Line2fit);
    x2 = min(coordinates2{i}(coordinates2{i}(:,4)~=Index(i),1)):0.005:max(coordinates2{i}(coordinates2{i}(:,4)~=Index(i),1));
    y2 = Line2(i,1)*x2 + Line2(i,2);
    plot(x2, y2 , '-b', 'LineWidth',3);
    cd(linedir);
    print(image2, [num2str(i),'_lines.tif'], '-dtiff', '-r150');
    cd(currdir); 
    close all;
    
    Distance(i) = abs(Line1(i,2)-Line2(i,2))/sqrt(Line2(i,1)*Line2(i,1) + 1);
    
end
Filenames = 1:1:numel(files);
Distance2 = [Filenames' Distance];
cd(resultdir);
Otput_summary = 'Distances.csv';
headers2 = {'Ring', 'Distance'};
csvwrite_with_headers(Otput_summary,Distance2, headers2);
cd(currdir); 

%% Moving, rotating and summing up
coordinates_rot1 = struct([]);
coordinates_all = zeros(1,3);
rotation = zeros(numel(files),1);
Center = zeros(numel(files),2);
for i=1:numel(files)
    rotation(i) = atan2d(Line1(i,1),1);
    Center(i,1) = mean(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1)) - ...
        sind(rotation(i))*(Line2(i,2)-Line1(i,2))/sqrt(Line2(i,1)*Line2(i,1) + 1)/2;
    Center(i,2) = Line2(i,1)*Center(i,1) + (Line1(i,2)+Line2(i,2))/2;
    
    R = [cosd(-rotation(i)) -sind(-rotation(i)); sind(-rotation(i)) cosd(-rotation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(coordinates2{i}))';
    
    coordinates_rot1{i} = coordinates2{i}(:,1:2) - Center_mat_border;
    coordinates_rot1{i} = (R* coordinates_rot1{i}')';
    coordinates_rot1{i}(:,3) = coordinates2{i}(:,3);
    
    coordinates_all = [coordinates_all;coordinates_rot1{i}(:,1:3)];
end
coordinates_all(coordinates_all(:,3) ==0,:) =[];

bin_size = 2;
dist_length = max(max(abs(coordinates_all(:,1:2))));
binrange = 0 : bin_size : (dist_length*2);
bincenter=binrange(1:(end-1)) + bin_size/2;
final_lines = zeros(length(bincenter)+2,length(bincenter)+2);
for n=1:length(coordinates_all)
   final_lines(ceil((coordinates_all(n,2)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_all(n,1)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1) = ...
       final_lines(ceil((coordinates_all(n,2)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_all(n,1)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1) + coordinates_all(n,3);
end

final = mtimes(sum(final_lines,2),sum(final_lines,1));
cd(resultdir);
image1 = figure;
imshow(final, [0, max(max(final))]);
print(image1, 'averaged_image.tif', '-dtiff', '-r150');
image2 = figure;
imshow(final_lines, [0, max(max(final_lines))]);
print(image2, 'summed_image.tif', '-dtiff', '-r150');
image3 = figure;
[fit_line, gof_line] = fit((1:bin_size:size(final_lines,2)*2)',sum(final_lines,2), 'gauss2');
plot((1:bin_size:size(final_lines,2)*2)', sum(final_lines,2), 'o',...
    (1:bin_size:size(final_lines,2)*2)', fit_line((1:bin_size:size(final_lines,2)*2)'),'r');
print(image3, 'fitted_image.tif', '-dtiff', '-r150');
cd(currdir);
clear variables
close all