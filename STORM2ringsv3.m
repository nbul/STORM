clc
clear variables
close all


%% Defining extension
% Default cutoff
cutoff =0.60;
cutoff2 =0.30;
Column_intensity = 18;
Ring = struct([]);
coordinates = struct([]);
Cluster = struct([]);
Centroid = struct([]);
coordinates2 = struct([]);
coordinates3 = struct([]);
coordinates4 = struct([]);
IntensityCluster = struct([]);


usedefault = questdlg(strcat('Use default settings: (Cut-off value = ', num2str(cutoff),'?);Cut-off value2 = ',...
    num2str(cutoff2),'?);(Intensity column = ',...
    num2str(Column_intensity),'?)'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No')
    parameters = inputdlg({'Enter cut-off value:','Enter cut-off value2:','Enter intensity column:'},'Parameters',1,...
        {num2str(cutoff),num2str(cutoff2), num2str(Column_intensity)});
    % Redefine extension
    cutoff = str2double(parameters{1});
    cutoff2 = str2double(parameters{2});
    Column_intensity = str2double(parameters{3});
else
    parameters{1} = num2str(cutoff);
    parameters{2} = num2str(cutoff2);
    parameters{3} = num2str(Column_intensity);
end

usedefault2 = questdlg(strcat('Do you want to rotate and crop the data?'),'Settings','Yes','No','No');
if strcmp(usedefault2, 'No')
    choice = 1;
else
    choice = 0;
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
clustdir = [filedir,'/Clustering'];

if exist([filedir,'/Lines'],'dir') == 0
    mkdir(filedir,'/Lines');
end
linedir = [filedir,'/Lines'];

if exist([filedir,'/Distributions'],'dir') == 0
    mkdir(filedir,'/Distributions');
end
distdir = [filedir,'/Distributions'];


Index = zeros(numel(files),1);
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
Center = zeros(numel(files),2);
rotation = zeros(numel(files),1);
%% Clustering
for l=1:numel(files)
    x1 = linspace(min(coordinates{l}(:,1)) - 2,max(coordinates{l}(:,1)) + 2,d);
    x2 = linspace(min(coordinates{l}(:,2)) - 2,max(coordinates{l}(:,2)) + 2,d);
    [x1grid,x2grid] = meshgrid(x1,x2);
    X0 = [x1grid(:) x2grid(:)];
    i = 0;
    t = 1;
    %% Splitting proportianally to the signal
    for o=1:length(Ring{l})
        m=ceil(Ring{l}(o,18)/100);
        for f=1:m
            coordinates4{l}(t,:) = coordinates{l}(o,:);
            t=t+1;
        end
    end
    while i == 0
        image1 = figure;
        axis equal; 
        gmfit = fitgmdist(coordinates4{l}(:,1:2),k,'CovarianceType',Sigma,...
            'SharedCovariance',SharedCovariance,'Options',options);
        clusterX = cluster(gmfit,coordinates4{l}(:,1:2));
        mahalDist = mahal(gmfit,X0);
        coordinates3{l} = [coordinates4{l} clusterX];
        h1 = gscatter(coordinates4{l}(:,1),coordinates4{l}(:,2),clusterX);
        hold on;
        for m = 1:k
            idx = mahalDist(:,m)<=threshold;
            Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
            h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
            uistack(h2,'bottom');
        end
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10);
        legend(h1);
        title(num2str(l));
        usedefault = questdlg(strcat('Are you happy about clustering)'),'Settings','Yes','No','Yes');
        if strcmp(usedefault, 'Yes')
            i=1;
            cd(clustdir);
            print(image1, [num2str(l),'_cluster.tif'], '-dtiff', '-r150');
            cd(currdir);
        end
        close all;
    end
    %% Identifying a brighter cluster for fit
    IntensityCluster{l}(1) = length(coordinates3{l}(coordinates3{l}(:,4)==1,3));
    IntensityCluster{l}(2) = length(coordinates3{l}(coordinates3{l}(:,4)==2,3));
    if IntensityCluster{l}(1)>IntensityCluster{l}(2)
        Index(l) = 1;
    else
        Index(l) = 2;
    end

    %% Determining rotation angle and shift
    mu=mean(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1:2),1);
    mu2=mean(coordinates3{l}(coordinates3{l}(:,4)~=Index(l),1:2),1);
    X_minus_mu=coordinates3{l}(coordinates3{l}(:,4)==Index(l),1:2)-...
        repmat(mu, size(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1:2),1), 1);
    Sigma2=(X_minus_mu'*X_minus_mu)/size(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1:2),1);
    [V, D]=eig(Sigma2);
    rotation(l) = atan2d(V(2,2), V(2,1));
    Center(l,:) = (mean(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1:2),1) +...
        mean(coordinates3{l}(coordinates3{l}(:,4)~=Index(l),1:2),1))/2;
    image12 = figure;
    
    h2 = gscatter(coordinates3{l}(:,1),coordinates3{l}(:,2),clusterX);
    hold on;
    x3 = min(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1)):0.005:max(coordinates3{l}(coordinates3{l}(:,4)==Index(l),1));
    y3 = tand(rotation(l))*x3 + (mu(2)-tand(rotation(l))*mu(1));
    h3 = plot(x3, y3 , '-k', 'LineWidth',3);
    hold on;
    x4 = min(coordinates3{l}(coordinates3{l}(:,4)~=Index(l),1)):0.005:max(coordinates3{l}(coordinates3{l}(:,4)~=Index(l),1));
    y4 = tand(rotation(l))*x4 + (mu2(2)-tand(rotation(l))*mu2(1));
    h4 = plot(x4, y4 , '-k', 'LineWidth',3);
    hold on;
    h5 = plot(Center(l,1),Center(l,2),'kx','LineWidth',2,'MarkerSize',10);
    axis equal;
    xlim([min(coordinates3{l}(:,1)) max(coordinates3{l}(:,1))]);
    ylim([min(coordinates3{l}(:,2)) max(coordinates3{l}(:,2))]);
    legend(h2);
    cd(clustdir);
    print(image12, [num2str(l),'_cluster_lines.tif'], '-dtiff', '-r150');
    cd(currdir);
    close all;
end



%% Moving, rotating and summing up
coordinates_rot1 = struct([]);
coordinates_all = zeros(1,3);


for i=1:numel(files)
    
    R = [cosd(-rotation(i)) -sind(-rotation(i)); sind(-rotation(i)) cosd(-rotation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(coordinates3{i}))';
    
    coordinates_rot1{i} = coordinates3{i}(:,1:2) - Center_mat_border;
    coordinates_rot1{i} = (R* coordinates_rot1{i}')';
    coordinates_rot1{i}(:,3) = coordinates3{i}(:,3);
    
    coordinates_all = [coordinates_all;coordinates_rot1{i}(:,1:3)]; %#ok<AGROW>
end
coordinates_all(coordinates_all(:,3) ==0,:) =[];

if choice == 0
    rotatecrop;
end

%% Individual distributions and distances
bin_size = 2;
dist_length = max(max(abs(coordinates_all(:,1:2))));
binrange = -dist_length : bin_size : dist_length;
bincenter=binrange(1:(end-1)) + bin_size/2;
Dist = zeros(length(bincenter),numel(files)+1);
Dist3 = zeros(length(bincenter),numel(files)+1);
curve = struct([]);
gof = struct([]);
Distance = zeros(numel(files),7);
cd(distdir);
warning('off','all');
for i=1:numel(files)   
    Dist(:,i)=histcounts(coordinates_rot1{i}(:,2),binrange)';
    Dist(:,i)=  Dist(:,i) / sum( Dist(:,i));
    Dist3(:,i)=histcounts(coordinates_rot1{i}(:,1),binrange)';
    Dist3(:,i)=  Dist3(:,i) / sum( Dist3(:,i));
    Y = Dist(:,i);
    X = bincenter';
    Dist2 = [X,Y];
    Dist2 = Dist2(find(Dist2(:,2),1,'first'):find(Dist2(:,2),1,'last'),:);
    options2 = fitoptions('gauss2','Lower', [max(Dist(:,i))/5 min(Dist2(:,1)) 15 max(Dist(:,i))/5 0 15],...
        'Upper', [max(Dist(:,i)) 0 max(Dist2(:,1))/2 max(Dist(:,i)) max(Dist2(:,1)) max(Dist2(:,1))/2],...
        'Robust','LAR');
    [curve{i},gof{i}] = fit(Dist2(:,1),Dist2(:,2),'gauss2',options2);
    image1 = figure;
    plot(bincenter', Dist(:,i), 'o',...
        bincenter', curve{i}(bincenter'),'r');
    title(num2str(gof{i}.rsquare));
    if gof{i}.rsquare>cutoff
        Distance(i,7) = 1;
    end 

    print(image1, [num2str(i),'_distribution.tif'], '-dtiff', '-r150');    
    Distance(i,1) = curve{i}.b1;
    Distance(i,2) = curve{i}.c1;
    Distance(i,3) = curve{i}.b2;
    Distance(i,4) = curve{i}.c2;
    Distance(i,5) = gof{i}.rsquare;
    Distance(i,6) = abs(curve{i}.b1-curve{i}.b2);
end
close all;
cd(currdir);
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
final_lines = imadjust(double(final_lines/max(final_lines(:))));
final_lines = imgaussfilt(final_lines,2);
final_lines = imadjust(final_lines);
image2 = figure;
imshow(final_lines, [0, max(max(final_lines))]);
print(image2, 'summed_image.tif', '-dtiff', '-r150');


for i=1:numel(files)
    if Distance(i,7) == 1
        Dist(:,numel(files)+1) = Dist(:,numel(files)+1) + Dist(:,i);
        Dist3(:,numel(files)+1) = Dist3(:,numel(files)+1) + Dist3(:,i);
    end
end
Dist(:,numel(files)+1)=  Dist(:,numel(files)+1) / sum( Dist(:,numel(files)+1));
Dist3(:,numel(files)+1)=  Dist3(:,numel(files)+1) / sum( Dist3(:,numel(files)+1));

final = mtimes(Dist(:,numel(files)+1),Dist3(:,numel(files)+1)');
cd(resultdir);
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


