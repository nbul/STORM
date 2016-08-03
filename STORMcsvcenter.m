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

for i=1:numel(files)
    clear test_k silh
    cd('STORMcsv/');
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4)];
    for q=6:15
        test_k(q-5,:) = kmeans(coordinates{i},q);
        silh(q-5) = mean(silhouette(coordinates{i},test_k(q-5,:)));
    end
    [num1 idx1] = max(silh);
    ClusterNumber(i) = idx1+7;
    counter=zeros(ClusterNumber(i),numel(files));
    IntensityCluster{i}=zeros(ClusterNumber(i),1);
    [Cluster{i},Centroid{i}]= kmeans(coordinates{i},ClusterNumber(i));
    Centroid2{i}= Centroid{i};
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

%% Individual distributions
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
    p{i} = p{i}/max(p{i});
    [curve{i},gof{i}] = fit(binsnew'*bin_size,p{i},'gauss1');
    radius(i) = curve{i}.b1;
    sigma1(i) = sqrt(curve{i}.c1*curve{i}.c1/2);
    pvalue(i) = gof{i}.rsquare;
    % compute theoretical probabilities
        if gof{i}.rsquare > cutoff
            counter = counter + 1; 
            good_ring(counter) = i;
            usage(i) = 1;
        end
    % plot data and theoretical distribution
    subplot(4, ceil(numel(files)/4),i);
    plot(binsnew'*bin_size, p{i}, 'o', binsnew'*bin_size, curve{i}(binsnew'*bin_size));
    title(num2str(gof{i}.rsquare));
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
p{numel(files)+1} = p{numel(files)+1}/max(p{numel(files)+1});
[curve{numel(files)+1},gof{numel(files)+1}] = fit(binsnew'*bin_size,p{numel(files)+1},'gauss1');
radius(numel(files)+1) = curve{numel(files)+1}.b1;
sigma1(numel(files)+1) = sqrt(curve{numel(files)+1}.c1*curve{numel(files)+1}.c1/2);
pvalue(numel(files)+1) = gof{numel(files)+1}.rsquare;
usage(numel(files)+1) = 1;
% plot data and theoretical distribution
subplot(4, ceil(numel(files)/4),numel(files)+1);
plot(binsnew'*bin_size, p{numel(files)+1}, 'o', binsnew'*bin_size, curve{numel(files)+1}(binsnew'*bin_size));
title(num2str(gof{numel(files)+1}.rsquare));

%% Curve fitting
fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [curve{numel(files)+1}.a1 sigma1(numel(files)+1) curve{numel(files)+1}.b1]);
ft = fittype('2*pi*A*exp(-(x*x + R*R)/(2*sigma*sigma))*besseli(0, (R*x/(sigma*sigma)))',...
         'coefficients',{'A','sigma', 'R'},'independent','x','options',fo);
[curve2,gof2] = fit(binsnew'*bin_size,p{numel(files)+1},ft);
subplot(4, ceil(numel(files)/4),numel(files)+2);
plot(binsnew'*bin_size, p{numel(files)+1}, 'o', binsnew'*bin_size, curve2(binsnew'*bin_size));
title(num2str(gof2.rsquare));


cd('STORMcsv/');
mkdir('result');
cd('result/');
print(image1, 'summarized_distributions.tif', '-dtiff', '-r150');
csvwrite('distributions.csv', peaks)
cd('../../');

%% Save means and SD
Number2 = [1:(numel(files)+1)];
peaks = [Number2' radius' sigma1' pvalue' usage; numel(files)+2 curve2.R curve2.sigma gof2.rsquare 1];

cd('STORMcsv/result/');
csvwrite('radiuses.csv', peaks);
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
            final_ring_fitted(l,n) = abs(curve2((ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)*bin_size)/max(curve2(binsnew'*bin_size)));
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

%% Pairwise angles between clusters
Angles_all = zeros(1,1);
for i=1:numel(files)
    CounterAngles = 0;
    for k=1:(length(Centroid{i})-1)
        V1 = [Centroid{i}(k,1)-centerX(i) Centroid{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid{i})
            V2 = [Centroid{i}(m,1)-centerX(i) Centroid{i}(m,2)-centerY(i) 0];
            CounterAngles = CounterAngles +1;
            Angles{i}(CounterAngles,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
    Angles_all = [Angles_all; Angles{i}];
end

binrangeAngles = [0 : binAnlges : 180];
bincenterAngles=binrangeAngles(1:(end-1)) + bin_size/2;
[NAngles, binsAngles] = histc(Angles_all,binrangeAngles);
figure; bar(binrangeAngles, NAngles);

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
