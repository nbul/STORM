%% Lead script to analyse rings
% The data should be in a single folder in .csv format 
% The output data:

% aligned_distribution.tif - summarised circular distribution of signal
% from center normalised to signal intensity and shifted to start with
% highest peak

% angles_clusters.csv - distributions of pair-wise angles between clusters
% identified with k-means after and before removing 2 weakest clusters

% distribution_circular.csv - individual circular distributions from center
% normalised to signal intensity and shifted to start with
% highest peak (first column is the bin, last column - summarised values)

% distribution_radial - individual radial distributions from center
% normalised to signal intensity (first column is the bin, last column -
% summarised values) 

% peaks_distribution.csv - destribution of all peaks found in individual
% circular destributions.

% peaks_distribution.tif - visualisation of peaks_distribution.csv.

% radiuses_curve.csv - data about individual distributions fitted with the
% equation from the paper (last raw is the summarised ring)

% radiuses.csv - data about individual distributions fitted with gaussian
% (last raw is the summarised ring)

% summarized_distributions_modified.tif - visualisation of
% radiuses_curve.csv 

% summarized_distributions.tif - visualisation of radiuses.csv

% summarized_ring_fitted.tif - visualisation of summarised ring fitted with
% the equation from the paper

% summarized_ring.tif - visualisation of non-fitted summarised ring



clc
clear variables
close all

bin_size = 2;
binAnlges = 10;
min_size = 50;


%% Defining extension
% Default cutoff
cutoff =0.80;
Column_intensity = 18;

usedefault = questdlg(strcat('Use default settings: (Cut-off value = ', num2str(cutoff),'?);(Intensity column = ', num2str(Column_intensity),'?)'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No');
    parameters = inputdlg({'Enter cut-off value:', 'Intensity column:'},'Parameters',1,{num2str(cutoff), num2str(Column_intensity)});
    % Redefine extension 
    cutoff = str2double(parameters{1});
    Column_intensity = str2double(parameters{2});
else
    parameters{1} = num2str(cutoff);
    parameters{2} = num2str(Column_intensity);
end

%% Determening paths
currdir = pwd;
filedir = uigetdir();
files = dir(strcat(filedir,'/*.csv'));
cd(filedir);
mkdir(filedir,'/results');
resultdir = [filedir, '/results'];
sum_ring = zeros(512,512);

%% Finding the center
cd(currdir);
STORMcsvcenter;

%% Individual distributions and gaussian

STORMgauss;

%% Curve fitting

STORMcurve;


%% Make rings and get radiuses

STORMrings;

%% Angle distributions

STORMcsvangle;

clc
clear variables
close all

%% Testing clustering 
%$First define number of ring by inputing it into variable i in command
%%window. e.g. i=5. Then uncomment the following part of code, copy to
%%command window and run by pressing enter.

% figure;
% l=0;
% for n=1:ClusterNumber(i)
%     if (n~=RemoveIdx1(i,1) && n~=RemoveIdx1(i,2))
%         plot(coordinates2{i}(coordinates2{i}(:,1)==n,2),coordinates2{i}(coordinates2{i}(:,1)==n,3),'*')
%         hold on
%         l=l+1;
%         plot(Centroid{i}(Centroid{i}(:,3)==n,1),Centroid{i}(Centroid{i}(:,3)==n,2),'k+','MarkerSize', 15)
%         hold on
%     end
% end
% rectangle('position',[centerX(i)-Rfit(i),centerY(i)-Rfit(i),Rfit(i)*2,Rfit(i)*2],...
%     'curvature',[1,1],'linestyle','-','edgecolor','r');
% hold off
% figure;
% for n=1:10
%     plot(coordinates{i}(Cluster{i}==n,1),coordinates{i}(Cluster{i}==n,2),'*')
%     hold on
% end
% hold off
