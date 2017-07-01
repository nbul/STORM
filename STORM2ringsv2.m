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


usedefault = questdlg(strcat('Use default settings: (Cut-off value = ', num2str(cutoff),'?);(Intensity column = ',...
    num2str(Column_intensity),'?)'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No');
    parameters = inputdlg({'Enter cut-off value:','Enter intensity column:'},'Parameters',1,...
        {num2str(cutoff), num2str(Column_intensity)});
    % Redefine extension
    cutoff = str2double(parameters{1});
    Column_intensity = str2double(parameters{2});
else
    parameters{1} = num2str(cutoff);
    parameters{2} = num2str(Column_intensity);
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

%% Determining rotation angle and shift
Center = zeros(numel(files),2);
rotation = zeros(numel(files),1);
for i=1:numel(files)
    mu=mean(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1:2),1);
    X_minus_mu=coordinates3{i}(coordinates3{i}(:,4)==Index(i),1:2)-...
        repmat(mu, size(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1:2),1), 1);
    Sigma=(X_minus_mu'*X_minus_mu)/size(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1:2),1);
    [V, D]=eig(Sigma);
    rotation(i) = atan2d(V(2,2), V(2,1));
    Center(i,:) = (mean(coordinates3{i}(coordinates3{i}(:,4)==Index(i),1:2),1) +...
        mean(coordinates3{i}(coordinates3{i}(:,4)~=Index(i),1:2),1))/2;
end


%% Moving, rotating and summing up
coordinates_rot1 = struct([]);
coordinates_all = zeros(1,3);


for i=1:numel(files)
    
    R = [cosd(-rotation(i)) -sind(-rotation(i)); sind(-rotation(i)) cosd(-rotation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(coordinates2{i}))';
    
    coordinates_rot1{i} = coordinates2{i}(:,1:2) - Center_mat_border;
    coordinates_rot1{i} = (R* coordinates_rot1{i}')';
    coordinates_rot1{i}(:,3) = coordinates2{i}(:,3);
    
    coordinates_all = [coordinates_all;coordinates_rot1{i}(:,1:3)];
end
coordinates_all(coordinates_all(:,3) ==0,:) =[];

%% Individual distributions and distances
bin_size = 2;
dist_length = max(max(abs(coordinates_all(:,1:2))));
binrange = 0 : bin_size : (dist_length*2);
bincenter=binrange(1:(end-1)) + bin_size/2;
Dist = zeros(length(bincenter),numel(files)+1);
Dist2 = zeros(length(bincenter),numel(files)+1);
curve = struct([]);
gof = struct([]);
curve1 = struct([]);
gof1 = struct([]);
curve2 = struct([]);
gof2 = struct([]);
image1 = figure;
Distance = zeros(numel(files),7);
options = fitoptions('gauss2','Lower', [0 0 0 0 dist_length 0],...
        'Upper', [Inf dist_length Inf Inf dist_length*2 Inf]);
for i=1:numel(files)
    Dist(:,i)=histcounts(coordinates_rot1{i}(:,2),binrange-dist_length)';
    Dist(:,i)=  Dist(:,i) / sum( Dist(:,i));
    Dist2(:,i)=histcounts(coordinates_rot1{i}(:,1),binrange-dist_length)';
    Dist2(:,i)=  Dist2(:,i) / sum( Dist2(:,i));
    [curve1{i},gof1{i}] = fit(bincenter',Dist(:,i),'gauss2',options);
    [curve2{i},gof2{i}] = fit(bincenter',Dist(end:-1:1,i),'gauss2',options);
    if gof1{i}.rsquare>gof2{i}.rsquare
        curve{i} = curve1{i};
        gof{i} = gof1{i};
    else
        curve{i} = curve2{i};
        gof{i} = gof2{i};
        Dist(:,i) = Dist(end:-1:1,i);
    end
    subplot(4, ceil(numel(files)/4),i);
    plot(bincenter', Dist(:,i), 'o',...
        bincenter', curve{i}(bincenter'),'r');
    title(num2str(gof{i}.rsquare));
    if gof{i}.rsquare>cutoff
        Distance(i,7) = 1;
    end 
    
    Distance(i,1) = curve{i}.b1;
    Distance(i,2) = curve{i}.c1;
    Distance(i,3) = curve{i}.b2;
    Distance(i,4) = curve{i}.c2;
    Distance(i,5) = gof{i}.rsquare;
    Distance(i,6) = abs(curve{i}.b1-curve{i}.b2);
end

Filenames = 1:1:numel(files);
Distance2 = [Filenames' Distance];
cd(resultdir);


%% Building final image and summary distribution

final_lines = zeros(length(bincenter)+2,length(bincenter)+2);
for n=1:length(coordinates_all)
   final_lines(ceil((coordinates_all(n,2)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_all(n,1)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1) = ...
       final_lines(ceil((coordinates_all(n,2)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_all(n,1)+max(max(abs(coordinates_all(:,1:2)))))/bin_size)+1) + coordinates_all(n,3);
end
image2 = figure;
imshow(final_lines, [0, max(max(final_lines))]);
print(image2, 'summed_image.tif', '-dtiff', '-r150');


for i=1:numel(files)
    if Distance(i,7) == 1
        Dist(:,numel(files)+1) = Dist(:,numel(files)+1) + Dist(:,i);
        Dist2(:,numel(files)+1) = Dist2(:,numel(files)+1) + Dist2(:,i);
    end
end
Dist(:,numel(files)+1)=  Dist(:,numel(files)+1) / sum( Dist(:,numel(files)+1));
Dist2(:,numel(files)+1)=  Dist2(:,numel(files)+1) / sum( Dist2(:,numel(files)+1));

final = mtimes(Dist(:,numel(files)+1),Dist2(:,numel(files)+1)');
cd(resultdir);
print(image1, 'Distributions.tif', '-dtiff', '-r150');
image3 = figure;
imshow(final, [0, max(max(final))]);
print(image3, 'averaged_image.tif', '-dtiff', '-r150');

image4 = figure;
[curve1{numel(files)+1},gof1{numel(files)+1}] = fit(bincenter',Dist(:,numel(files)+1),'gauss2',options);
[curve2{numel(files)+1},gof2{numel(files)+1}] = fit(bincenter',Dist(end:-1:1,numel(files)+1),'gauss2',options);
if gof1{numel(files)+1}.rsquare>=gof2{numel(files)+1}.rsquare
    curve{numel(files)+1} = curve1{numel(files)+1};
    gof{numel(files)+1} = gof1{numel(files)+1};
else
    curve{numel(files)+1} = curve2{numel(files)+1};
    gof{numel(files)+1} = gof2{numel(files)+1};
    Dist(:,numel(files)+1) = Dist(end:-1:1,numel(files)+1);
end
plot(bincenter', Dist(:,numel(files)+1), 'o',...
        bincenter', curve{numel(files)+1}(bincenter'),'r');
title(num2str(gof{numel(files)+1}.rsquare));
print(image4, 'average_distribution.tif', '-dtiff', '-r150');

Distance2(numel(files)+1,2) = curve{numel(files)+1}.b1;
Distance2(numel(files)+1,3) = curve{numel(files)+1}.c1;
Distance2(numel(files)+1,4) = curve{numel(files)+1}.b2;
Distance2(numel(files)+1,5) = curve{numel(files)+1}.c2;
Distance2(numel(files)+1,6) = gof{numel(files)+1}.rsquare;
Distance2(numel(files)+1,7) = abs(curve{numel(files)+1}.b1-curve{numel(files)+1}.b2);
Otput_summary = 'Distances.csv';
headers2 = {'Ring', 'Center1', 'Width1', 'Center2', 'Width2', 'gof','Distance','Usage'};
csvwrite_with_headers(Otput_summary,Distance2, headers2);

cd(currdir);
clear variables
close all


