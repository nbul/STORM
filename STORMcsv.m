clear all
close all

bin_size = 2;
binAnlges = 8;
min_size = 50;
cutoff =0.90;
Column_intensity = 18;

cd('STORMcsv/');
files = dir('*.csv');
cd('../');
sum_ring = zeros(512,512);

%% Finding the center

STORMcsvcenter;

%% Individual distributions and gaussian

STORMgauss;

%% Curve fitting

STORMcurve;


%% Save means and SD
Number2 = [1:(numel(files)+1)];
peaks = [Number2' radius' sigma1' pvalue' usage];
for i=1:(numel(files)+1)
peaks2(i,:) = [i curve2{i}.R curve2{i}.sigma gof2{i}.rsquare usage(i)];
end
cd('STORMcsv/result/');
csvwrite('radiuses.csv', peaks);
cd('../../');

cd('STORMcsv/result/');
csvwrite('radiuses_curve.csv', peaks2);
cd('../../');


%% Making summarized ring
final_ring = zeros(length(bincenter)*2,length(bincenter)*2);
final_ring_fitted = zeros(length(bincenter)*2,length(bincenter)*2);
[columnsInImage rowsInImage] = meshgrid(1:(length(bincenter)*2), 1:(length(bincenter)*2));

centerX2 = length(bincenter);
centerY2 = length(bincenter);
sum_intensity = [sum_intensity; zeros(ceil(sqrt(2*length(bincenter)*length(bincenter))-length(sum_intensity))+1, 1)];
for l=1:(length(bincenter)*2)
    for n = 1:(length(bincenter)*2)
            final_ring(l,n) = abs(sum_intensity(ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)/max(sum_intensity));
            final_ring_fitted(l,n) = abs(curve2{numel(files)+1}((ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)*bin_size)...
                /max(curve2{numel(files)+1}(binsnew'*bin_size)));
    end
end

image3 = figure;
imshow(final_ring);
image4 = figure;
imshow(final_ring_fitted);
cd('STORMcsv/result/');
print(image3, 'summarized_ring.tif', '-dtiff', '-r150');
print(image4, 'summarized_ring_fitted.tif', '-dtiff', '-r150');
cd('../../');

%STORMcsvangle;

% % Testing clustering
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
