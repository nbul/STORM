%% Lead script to analyse rings
% The data should be in a single folder in .csv format 

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


%% Save means and SD

peaks2 = zeros(numel(files)+1,6);
Number2 = 1:(numel(files)+1);
ClusterNumber = [ClusterNumber'; 0];
peaks = [Number2' radius' sigma1' pvalue' ClusterNumber usage];
for i=1:(numel(files)+1)
peaks2(i,:) = [i curve2{i}.R curve2{i}.sigma gof2{i}.rsquare ClusterNumber(i) usage(i)];
end

cd(resultdir);
headers = {'ring', 'radius', 'width', 'p-value', 'Number of clusters','Used in sum?'};
csvwrite_with_headers('radiuses.csv', peaks, headers);
csvwrite_with_headers('radiuses_curve.csv', peaks2, headers);
cd(currdir);


%% Making summarized ring
% creating empty image and specifying its center
final_ring = zeros(length(bincenter)*2,length(bincenter)*2);
final_ring_fitted = zeros(length(bincenter)*2,length(bincenter)*2);
[columnsInImage, rowsInImage] = meshgrid(1:(length(bincenter)*2), 1:(length(bincenter)*2));
centerX2 = length(bincenter);
centerY2 = length(bincenter);

%Drowing the ring
sum_intensity = [sum_intensity; zeros(ceil(sqrt(2*length(bincenter)*length(bincenter))-length(sum_intensity))+1, 1)];
for l=1:(length(bincenter)*2)
    for n = 1:(length(bincenter)*2)
            final_ring(l,n) = abs(sum_intensity(ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)/max(sum_intensity));
            final_ring_fitted(l,n) = abs(curve2{numel(files)+1}((ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)*bin_size)...
                /max(curve2{numel(files)+1}(binsnew'*bin_size)));
    end
end

%Saving image
image3 = figure;
imshow(final_ring);
image4 = figure;
imshow(final_ring_fitted);
cd(resultdir);
print(image3, 'summarized_ring.tif', '-dtiff', '-r150');
print(image4, 'summarized_ring_fitted.tif', '-dtiff', '-r150');
cd(currdir);

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
