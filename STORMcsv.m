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
cutoff2 = 0.98;
Column_intensity = 18;

usedefault = questdlg(strcat('Use default settings: (Cut-off value = ', num2str(cutoff),'?);(Cut-off value for circle= ',...
    num2str(cutoff2),'?);(Intensity column = ', num2str(Column_intensity),'?)'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No');
    parameters = inputdlg({'Enter cut-off value:','Enter cut-off value for circle', 'Intensity column:'},...
        'Parameters',1,{num2str(cutoff), num2str(cutoff2), num2str(Column_intensity)});
    % Redefine extension
    cutoff = str2double(parameters{1});
    cutoff2 = str2double(parameters{2});
    Column_intensity = str2double(parameters{3});
else
    parameters{1} = num2str(cutoff);
    parameters{2} = num2str(cutoff2);
    parameters{3} = num2str(Column_intensity);
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
sum_ring = zeros(512,512);
usage = zeros(numel(files)+1,1);

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

