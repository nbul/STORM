clear all
close all

pixelsize = 0.8;
fudgeFactor = 0.2; %rerun edge detector with fudgefactor
%area1 = 10000; %cutoff area to remove unneeded objects
Smoothing = 4; % Gaussian smouthing diameter
Noise = 11; % Noise reduction diameter
Dil_factor = 11; % Se values for Dialtion
Er_factor = 3; % Se value for Erosion
Cutof_Ecc = 0.75; %Cut off Eccentricity value
Cutof_Area = 30000; %Cut off Area value
Inc = 1.2; %Increase in best-fit ellipse diameter
cd('STORM tifs/small');
files = dir('*.tif');
cd('../../');

majorAxisLengths = zeros(numel(files),1);
center_x = zeros(numel(files),1);
center_y = zeros(numel(files),1);
im_x = zeros(numel(files),1);
im_y = zeros(numel(files),1);

%% Getting individual rings and placing them in center
for i=1:numel(files)
    cd('STORM tifs/small');
    Number1 = [num2str(i),'.tif'];
    Ring{i} = imread(Number1);
    imbackground = imopen(Ring{i},strel('disk',60));  
    imcorrected = Ring{i} - imbackground;
    cd('../../');
    Signal1 = imgaussfilt(imcorrected,Smoothing);
    Signal4 = wiener2(Signal1,[Noise Noise]);
    [junk threshold] = edge(Signal4, 'sobel'); %detect edges of i to estimate threshold
    
    BW = edge(Signal4,'sobel', threshold * fudgeFactor);
    %figure, imshow(BW);
    
    se90 = strel('line', Dil_factor, 90); %structuring elements the amount we are averaging over in x direction
    se0 = strel('line', Dil_factor, 0); %as above in y direction -make each element a 3x3 box alter size of box as necessary
    BWalldil = imdilate(BW, [se90 se0]);
    BWallfill = imfill(BWalldil, 'holes');
    BWallfill = imdilate(BWallfill, [se90 se0]);
    BWallfill = imfill(BWalldil, 'holes');
    %figure, imshow(BWallfill)
    [im_x(i) im_y(i)] = size(Ring{i});
    seD = strel('diamond',Er_factor);
    BWallclean = imerode(BWallfill,seD); %did erosion once
    BWallclean = imerode(BWallclean,seD);
    %figure, imshow(BWallclean), title('segmented image');
    o = regionprops (uint8( BWallclean), 'majorAxisLength', 'centroid');
    majorAxisLengths(i) = cat(1, o.MajorAxisLength);
    center = cat(1, o.Centroid);
    center_x(i) = center(1);
    center_y(i) = center(2);
end

%% Getting circular profiles of individual rings and their sum

for i=1:numel(files)
    intensity_total{i} = zeros(ceil(majorAxisLengths(i)/2),1);
    for k=1:360
        x = ceil(center_x(i)) + ceil(majorAxisLengths(i)*cosd(k)/2);
        y = ceil(center_y(i)) + ceil(majorAxisLengths(i)*sind(k)/2);
        if x<1 x2 = [1 ceil(center_x(i))];
        elseif x>(im_x(i)-1) x2 = [ceil(center_x(i)) im_x(i)];
        else x2 = [ceil(center_x(i)) x];
        end
        if y<1 y2 = [1 ceil(center_y(i))];
        elseif y>(im_y(i)-1) y2 = [ceil(center_y(i)) im_y(i)];
        else y2 = [ceil(center_y(i)) y];
        end
        x3 = [ceil(im_x(i)/2) ceil( im_x(i)/2) + ceil(majorAxisLengths(i)*cosd(k)/2)];
        y3 = [ceil(im_y(i)/2) ceil( im_y(i)/2) + ceil(majorAxisLengths(i)*sind(k)/2)];
        [cx{k},cy{k},intensity{k}] = improfile(Ring{i}, x2, y2); 
        intensity{k}(isnan(intensity{k})) = 0 ;
        [cx2{k},cy2{k},intensity2{k}] = improfile(Ring{i}, x3, y3);
        intensity{k} = [intensity{k}; zeros((length(intensity2{k})-length(intensity{k})),1)]; 
        intensity{k} = imresize(intensity{k}, [ceil(majorAxisLengths(i)/2) 1]);
        intensity_total{i} = intensity_total{i} + intensity{k};    
    end
    intensity_total{i} = intensity_total{i}/360;
    
    %figure, plot(1:1:length(intensity_total{i}), intensity_total{i});
   
end



%% Getting the radius individual rings and selecting rings
image3 = figure;
counter = 0;

usage = zeros(numel(files)+1,1);
for i=1:numel(files)
    bins = [1:ceil(majorAxisLengths(i)/2)];
    % normalize y to be a probability (sum = 1)
    p{i} = intensity_total{i}/ sum(intensity_total{i});
    % compute weighted mean and standard deviation
    m{i} = sum(bins' .* p{i});
    s{i} = sqrt(sum((bins' - m{i}) .^ 2 .* p{i}));
    % compute theoretical probabilities
    pth{i} = normpdf(bins, m{i}, s{i});
    resudials{i} = p{i} - pth{i}';
    chi_square(i) = sum( resudials{i}.*resudials{i})/max(pth{i});
    pvalue(i) = 1-chi2cdf(chi_square(i), 2);
        if pvalue(i) > 0.98
            counter = counter + 1; 
            good_ring(counter) = i;
            usage(i) = 1;
        end
    % plot data and theoretical distribution
    subplot(4, ceil(numel(files)/4),i);
    plot(bins', p{i}, 'o', bins', pth{i});
    title(num2str(pvalue(i)));
    hold on;
    %[bestfit,resid]=nlinfit(y, x, norm_func, initGuess);
end

%% Getting summarized distribution
sum_intensity = zeros(ceil(max(majorAxisLengths/2)),1);


for i=1:counter
    intensity_total{good_ring(i)} = [intensity_total{good_ring(i)}; zeros(ceil(max(majorAxisLengths/2))-length(intensity_total{good_ring(i)}),1)];
    sum_intensity = sum_intensity + intensity_total{good_ring(i)};
end
sum_intensity = sum_intensity / counter;

bins = [1:ceil(max(majorAxisLengths/2))];

p{numel(files)+1} = sum_intensity / sum(sum_intensity);

m{numel(files)+1} = sum(bins' .* p{numel(files)+1});
s{numel(files)+1} = sqrt(sum((bins' - m{numel(files)+1}) .^ 2 .* p{numel(files)+1}));
pth{numel(files)+1} = normpdf(bins, m{numel(files)+1}, s{numel(files)+1});
resudials{numel(files)+1} = p{numel(files)+1} - pth{numel(files)+1}';
chi_square(numel(files)+1) = sum( resudials{numel(files)+1}.*resudials{numel(files)+1})/max(pth{numel(files)+1});
pvalue(numel(files)+1) = 1-chi2cdf(chi_square(numel(files)+1), 2);

% plot data and theoretical distribution
subplot(4, ceil(numel(files)/4),numel(files)+1);
plot(bins', p{numel(files)+1}, 'o', bins', pth{numel(files)+1});
title(num2str(pvalue(numel(files)+1)));

cd('STORM tifs/');
print(image3, 'summarized_distributions.tif', '-dtiff', '-r150');
cd('../');

%Save means and SD
Number2 = [1:(numel(files)+1)];
m2 = cell2mat(m);
s2 = cell2mat(s);
distributionsall= cell2mat(intensity_total);
peaks = [Number2' m2'*pixelsize s2'*pixelsize pvalue' usage];

cd('STORM tifs/');
csvwrite('radiuses.csv', peaks);
csvwrite('distributions.csv', distributionsall);
cd('../');

%% Making summarized ring
final_ring = zeros(ceil(max(majorAxisLengths/2))*2,ceil(max(majorAxisLengths/2))*2);
[columnsInImage rowsInImage] = meshgrid(1:ceil(max(majorAxisLengths/2))*2, 1:ceil(max(majorAxisLengths/2))*2);

for r=ceil(max(majorAxisLengths/2)):-1:1
    centerX = ceil(max(majorAxisLengths/2));
    centerY = ceil(max(majorAxisLengths/2));
    circlePixels = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= r.^2;
    for l=1:(ceil(max(majorAxisLengths/2))*2)
        for n = 1:(ceil(max(majorAxisLengths/2))*2)
            if circlePixels(l,n) == 1
                final_ring(l,n) = abs(sum_intensity(r)/max(sum_intensity));
            end
        end
    end
end

image2 = figure;
imshow(final_ring);
cd('STORM tifs/');
print(image2, 'summarized_ring.tif', '-dtiff', '-r150');
cd('../');


%% Getting intensity at perifery for checking period
image4 = figure;

Thetas = 0:1:360;
for i=1:counter
    clear pr_temp maximum idx tempx1 tempx2 tempy1 tempy2
    ring2{i} = zeros(im_x(i), im_y(i));
    if center_x(i)<ceil(max(majorAxisLengths/2)+5)
        tempx1 = zeros(ceil(max(majorAxisLengths/2)+5)-ceil(center_x(i)), im_y(i));
        center_x(i) = center_x(i) + ceil(max(majorAxisLengths/2)+5)-ceil(center_x(i));
        ring2{i} = [tempx1; Ring{good_ring(i)}];
    
    elseif (ceil(center_x(i))+ceil(max(majorAxisLengths/2)+5))>im_x(i)
        tempx2 = zeros((ceil(center_x(i))+ceil(max(majorAxisLengths/2)+5))-im_x(i), im_y(i));
        ring2{i} = [Ring{good_ring(i)}; tempx2];
    else
        ring2{i} = Ring{good_ring(i)}; 
    end
    
    if ceil(center_y(i))<ceil(max(majorAxisLengths/2)+5)
        tempy1 = zeros(size(ring2{i},1), ceil(max(majorAxisLengths/2)+5)-ceil(center_y(i)));
        center_y(i) = center_y(i) + ceil(max(majorAxisLengths/2)+5)-ceil(center_y(i));
        ring2{i} = [tempy1, ring2{i}];
    end
    
    if (ceil(center_y(i))+ceil(max(majorAxisLengths/2)+5))>im_y(i)
        tempy2 = zeros(size(ring2{i},1), (ceil(center_y(i))+ceil(max(majorAxisLengths/2)+5))-im_y(i));
        ring2{i} = [ring2{i},tempy2];
    end
    
    
    for k=-3:1:3
        Circle_radius = m{good_ring(i)}+k;
        outline_x1 = center_x(i)+Circle_radius*cosd(Thetas);
        outline_y1 = center_y(i)+Circle_radius*sind(Thetas);
        pr_temp{k+4} = improfile(ring2{i}, outline_x1, outline_y1);  
    end
    pr{i} = zeros(360, 1);
    for k=1:7
        pr_temp{k} = imresize(pr_temp{k}, [360 1]);
        pr{i} = pr{i} + pr_temp{k};
    end
    pr{i} = pr{i}/k;
    [num idx] = max([pr{i}]);
    for s=1:360
        if s<idx
            pr_reshaped{i}(360+s+1-idx,1) = pr{i}(s,1);
        else
            pr_reshaped{i}(s+1-idx,1) = pr{i}(s,1);
        end
    end
    subplot(4, ceil(numel(files)/4),i);
    plot(pr_reshaped{i});   
end

profile_sum = zeros(360, 1);
for k=1:counter
    profile_sum = profile_sum + [pr_reshaped{k}];
end
profile_sum = profile_sum/numel(files);
subplot(4, ceil(numel(files)/4),numel(files)+1);
plot(profile_sum);
cd('STORM tifs/');
print(image4, 'summarized_intensity.tif', '-dtiff', '-r150');
cd('../');
