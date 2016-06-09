clear all
close all

bin_size = 2;
min_size = 50;
cutoff =0.9;
Column_intensity = 18;

cd('STORMcsv/');
files = dir('*.csv');
cd('../');
sum_ring = zeros(512,512);
%% Finding the center
for i=1:numel(files)
    clear test_k silh
    cd('STORMcsv/');
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4)];
    for q=8:12
        test_k(q-7,:) = kmeans(coordinates{i},q);
        silh(q-7) = mean(silhouette(coordinates{i},test_k(q-7,:)));
    end
    [num1 idx1] = max(silh);
    ClusterNumber(i) = idx1+7;
    counter=zeros(ClusterNumber(i),numel(files));
    IntensityCluster{i}=zeros(ClusterNumber(i),1);
    [Cluster{i},Centroid{i}]= kmeans(coordinates{i},ClusterNumber(i));
    coordinates2{i} = [Cluster{i} coordinates{i}];
    for k=1:ClusterNumber(i)
        IntensityCluster{i}(k) = sum(Ring{i}(Cluster{i}(:,1)==k,Column_intensity));
        for l=1:length(Cluster{i})
            if Cluster{i}(l)==k
                counter(k,i)=counter(k,i)+1;
            end
        end
    end
    cd('../');
    ClusterNames = [1:1:ClusterNumber(i)];
    Centroid{i} = [Centroid{i} ClusterNames'];
    [RemoveNum1 RemoveIdx1(i,1)] = min(IntensityCluster{i});
    IntensityCluster2{i}=IntensityCluster{i};
    IntensityCluster2{i}(RemoveIdx1(i,1),:) = [];
    Centroid{i}(RemoveIdx1(i,1),:) = [];
    coordinates2{i}(any(coordinates2{i}==RemoveIdx1(i,1),2),:) = [];
    [RemoveNum1 Idxtemp] = min(IntensityCluster2{i});
    RemoveIdx1(i,2) = Centroid{i}(Idxtemp,3);
    Centroid{i}(Idxtemp,:) = [];
    coordinates2{i}(any(coordinates2{i}==RemoveIdx1(i,2),2),:) = [];
    
    [centerX(i),centerY(i),Rfit(i)] = circfit(Centroid{i}(:,1),Centroid{i}(:,2));
    distances{i} = sqrt((Ring{i}(:,3)-centerX(i)).*(Ring{i}(:,3)-centerX(i)) + (Ring{i}(:,4)-centerY(i)).*(Ring{i}(:,4)-centerY(i)));
end

dist_length = ceil(max(Rfit)*2);
binrange = [0 : bin_size : dist_length];
bincenter=binrange(1:(end-1)) + bin_size/2;

for i=1:numel(files)
    clear distances_indexed distances_added
    [N, bins] = histc(distances{i},binrange);
    distances_indexed(:,2) = Ring{i}(:,Column_intensity);
    distances_indexed(:,1) = bins;
    distances_added = zeros(length(bincenter),2);
    distances_added(:,1) = bincenter;
    for k=1:length(bincenter)
        for ii = 1:length(Ring{i}(:,Column_intensity))
            if distances_indexed(ii,1) == k
                distances_added(k,2) = distances_added(k,2) + distances_indexed(ii,2);
            end
        end
    end
    
    total = sum(distances_added);
    for k=1:length(bincenter)
        distances_added_norm(k,1) = (k-1)* bin_size + bin_size/2;
    end
    distances_added_norm(:,i+1) = distances_added(:,2)/total(2);
end

%% Getting the radius individual rings and selecting rings
image1 = figure;
counter = 0;
binsnew = [1:length(bincenter)];
usage = zeros(numel(files)+1,1);
for i=1:numel(files)
    % normalize y to be a probability (sum = 1)
    p{i} = distances_added_norm(:,i+1);
    % compute weighted mean and standard deviation
    m{i} = sum(binsnew' .* p{i});
    s{i} = sqrt(sum((binsnew' - m{i}) .^ 2 .* p{i}));
    % compute theoretical probabilities
    pth{i} = normpdf(binsnew, m{i}, s{i});
    resudials{i} = p{i} - pth{i}';
    chi_square(i) = sum( resudials{i}.*resudials{i})/max(pth{i});
    pvalue(i) = 1-chi2cdf(chi_square(i), 2);
        if pvalue(i) > cutoff
            counter = counter + 1; 
            good_ring(counter) = i;
            usage(i) = 1;
        end
    % plot data and theoretical distribution
    subplot(4, ceil(numel(files)/4),i);
    plot(binsnew'*bin_size, p{i}, 'o', binsnew'*bin_size, pth{i});
    title(num2str(pvalue(i)));
    hold on;
    %[bestfit,resid]=nlinfit(y, x, norm_func, initGuess);
end        
%% Getting summarized distribution
sum_intensity = zeros(length(bincenter),1);


for i=1:counter
    sum_intensity = sum_intensity + distances_added_norm(:,(good_ring(i)+1));
end
sum_intensity = sum_intensity / counter;
distances_added_norm(:,numel(files)+1) = sum_intensity;


p{numel(files)+1} = sum_intensity / sum(sum_intensity);

m{numel(files)+1} = sum(binsnew' .* p{numel(files)+1});
s{numel(files)+1} = sqrt(sum((binsnew' - m{numel(files)+1}) .^ 2 .* p{numel(files)+1}));
pth{numel(files)+1} = normpdf(binsnew, m{numel(files)+1}, s{numel(files)+1});
resudials{numel(files)+1} = p{numel(files)+1} - pth{numel(files)+1}';
chi_square(numel(files)+1) = sum( resudials{numel(files)+1}.*resudials{numel(files)+1})/max(pth{numel(files)+1});
pvalue(numel(files)+1) = 1-chi2cdf(chi_square(numel(files)+1), 2);
% plot data and theoretical distribution
subplot(4, ceil(numel(files)/4),numel(files)+1);
plot(binsnew'*bin_size, p{numel(files)+1}, 'o', binsnew'*bin_size, pth{numel(files)+1});
title(num2str(pvalue(numel(files)+1)));

cd('STORMcsv/');
mkdir('result');
cd('result/');
print(image1, 'summarized_distributions.tif', '-dtiff', '-r150');
csvwrite('distributions.csv', peaks)
cd('../../');

%% Save means and SD
Number2 = [1:(numel(files)+1)];
m2 = cell2mat(m);
s2 = cell2mat(s);
peaks = [Number2' m2'*bin_size s2'*bin_size pvalue' usage];

cd('STORMcsv/result/');
csvwrite('radiuses.csv', peaks);
cd('../../');

%% Making summarized ring
final_ring = zeros(length(bincenter)*2,length(bincenter)*2);
final_ring_fitted = zeros(length(bincenter)*2,length(bincenter)*2);
[columnsInImage rowsInImage] = meshgrid(1:(length(bincenter)*2), 1:(length(bincenter)*2));

for r=length(bincenter):-1:1
    centerX2 = length(bincenter);
    centerY2 = length(bincenter);
    circlePixels = (rowsInImage - centerY2).^2 + (columnsInImage - centerX2).^2 <= r.^2;
    for l=1:(length(bincenter)*2)
        for n = 1:(length(bincenter)*2)
            if circlePixels(l,n) == 1
                final_ring(l,n) = abs(sum_intensity(r)/max(sum_intensity));
                final_ring_fitted(l,n) = abs(pth{numel(files)+1}(r)/max(pth{numel(files)+1}));
            end
        end
    end
end

image2 = figure;
imshow(final_ring);
image3 = figure;
imshow(final_ring_fitted);
cd('STORMcsv/result/');
print(image2, 'summarized_ring.tif', '-dtiff', '-r150');
print(image3, 'summarized_ring_fitted.tif', '-dtiff', '-r150');
cd('../../');

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
