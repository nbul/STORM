clear all
close all

bin_size = 2;
%dist_length = 250;
cutoff =0.95;
Column_intensity = 5;

cd('STORM csv/');
files = dir('*.csv');
cd('../');
sum_ring = zeros(512,512);

for i=1:numel(files)
    cd('STORM csv/');
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    cd('../');
    Ring_binned = hist3(Ring{i}(:,3:4), {min(Ring{i}(:,3)) : bin_size : max(Ring{i}(:,3)) min(Ring{i}(:,4)) : bin_size : max(Ring{i}(:,4))});
    Ring_binned_bw = im2bw(Ring_binned);
    o = regionprops (uint8(Ring_binned_bw), 'centroid');
    center = cat(1, o.Centroid);
    centerX = min(Ring{i}(:,3)) + center(1) * bin_size;
    centerY = min(Ring{i}(:,4)) + center(2) * bin_size;
    distances{i} = sqrt((Ring{i}(:,3)-centerX).*(Ring{i}(:,3)-centerX) + (Ring{i}(:,4)-centerY).*(Ring{i}(:,4)-centerY));
end
for i=1:numel(files)
    max_temp(i) = max(distances{i});
end

dist_length = ceil(max(max_temp)/2)+1;
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
    temp = zeros(dist_length-length(bincenter),1);
    for k=1:dist_length
        distances_added_norm(k,1) = (k-1)* bin_size + bin_size/2;
    end
    distances_added_norm(:,i+1) = [distances_added(:,2)/total(2); temp];
end
        
%% Getting the radius individual rings and selecting rings
image1 = figure;
counter = 0;
binsnew = [1:dist_length];
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
    plot(binsnew', p{i}, 'o', binsnew', pth{i});
    title(num2str(pvalue(i)));
    hold on;
    %[bestfit,resid]=nlinfit(y, x, norm_func, initGuess);
end        
%% Getting summarized distribution
sum_intensity = zeros(dist_length,1);


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
plot(binsnew', p{numel(files)+1}, 'o', binsnew', pth{numel(files)+1});
title(num2str(pvalue(numel(files)+1)));

cd('STORM csv/');
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

cd('STORM csv/result/');
csvwrite('radiuses.csv', peaks);
cd('../../');

%% Making summarized ring
final_ring = zeros(dist_length*2,dist_length*2);
final_ring_fitted = zeros(dist_length*2,dist_length*2);
[columnsInImage rowsInImage] = meshgrid(1:(dist_length*2), 1:(dist_length*2));

for r=dist_length:-1:1
    centerX2 = dist_length;
    centerY2 = dist_length;
    circlePixels = (rowsInImage - centerY2).^2 + (columnsInImage - centerX2).^2 <= r.^2;
    for l=1:(dist_length*2)
        for n = 1:(dist_length*2)
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
cd('STORM csv/result/');
print(image2, 'summarized_ring.tif', '-dtiff', '-r150');
print(image3, 'summarized_ring_fitted.tif', '-dtiff', '-r150');
cd('../../');