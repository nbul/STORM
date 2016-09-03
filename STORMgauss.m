%% Individual distributions and gaussian

%% Assign memory
dist_length = ceil(max(Rfit)*2);
binrange = 0 : bin_size : dist_length;
bincenter=binrange(1:(end-1)) + bin_size/2;
distances_added_norm = double(zeros(length(bincenter),numel(files)+1));

p = struct([]);
curve = struct([]);
gof = struct([]);

radius = zeros(1, numel(files));
sigma1 = zeros(1, numel(files));
pvalue = zeros(1, numel(files));
good_ring = zeros(1, numel(files));

%% Obtaining distributions of signal from the center
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
binsnew = 1:length(bincenter);
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

%% Writing output file with distributions
cd(resultdir);
print(image1, 'summarized_distributions.tif', '-dtiff', '-r150');
csvwrite('distributions_radial.csv', distances_added_norm);
cd(currdir);