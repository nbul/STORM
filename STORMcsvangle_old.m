%% This script tests periodicity of clusters

%% Assign memory
Anglenumber1 = 0;
for i=1:numel(files)
    Anglenumber1 = Anglenumber1 + (ClusterNumber(i)-3)*(ClusterNumber(i)-2)/2;
end
Anglenumber2 = 0;
for i=1:numel(files)
    Anglenumber2 = Anglenumber2 + (ClusterNumber(i)-1)*(ClusterNumber(i))/2;
end
Angles_clean = zeros(Anglenumber1,1);
Angles_all = zeros(Anglenumber2,1);
binrangeAngles = 0 : binAnlges : 180;
bincenterAngles=binrangeAngles(1:(end-1)) + bin_size/2;

peakindex = zeros(1, numel(files));
shiftposition = zeros(1, numel(files));

angle = struct([]);
Anglesind = struct([]);
PeakPos = struct([]);
PeakMag = struct([]);
PeakPos_reshaped = struct([]);
PeakMag_reshaped = struct([]);

anglerange = 0 : bin_size : 360;
all_peaks_range = 0 : binAnlges : 360;

distributions = zeros(length(anglerange), numel(files)+1);
distributions_reshaped = zeros(length(anglerange)-1, numel(files)+2);
profile_all = zeros(180,1);


%% Pairwise angles between clusters (with removed 2 clusters from each ring)
counter = 0;
for i=1:numel(files)
    CounterAngles = 0;
    for k=1:(length(Centroid{i})-1)
        V1 = [Centroid{i}(k,1)-centerX(i) Centroid{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid{i})
            counter = counter + 1;
            V2 = [Centroid{i}(m,1)-centerX(i) Centroid{i}(m,2)-centerY(i) 0];
            Angles_clean(counter,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
    
end

%% Pairwise angles between clusters (all)
counter = 0;
for i=1:numel(files)
    CounterAngles2 = 0;
    for k=1:(length(Centroid2{i})-1)
        V1 = [Centroid2{i}(k,1)-centerX(i) Centroid2{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid2{i})
            counter = counter + 1;
            V2 = [Centroid2{i}(m,1)-centerX(i) Centroid2{i}(m,2)-centerY(i) 0];
            Angles_all(counter,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
end

%% Binning and histogram
[NAngles, binsAngles] = histc(Angles_clean,binrangeAngles);
NAngles2 = histc(Angles_all,binrangeAngles);

Data_angles = [binrangeAngles' NAngles NAngles2];

%% Alternative approach by making circular distribution

%Collecting angles of each event relative to center and shift it within 0-180°
%combining with intensity

for i=1:numel(files)
    for m=1:length(Ring{i}(:,3))
        if Ring{i}(m,3)>centerX(i)
            angle{i}(m,1)= atand((Ring{i}(m,3)-centerX(i))./(Ring{i}(m,4)-centerY(i)));
        else angle{i}(m,1)= 180 +atand((Ring{i}(m,3)-centerX(i))./(Ring{i}(m,4)-centerY(i)));
        end
    end
    
    for k=1:length(angle{i})
        if angle{i}(k,1)<0 
            angle{i}(k,1) = 270 - angle{i}(k,1);
        end
    end
angle{i}(:,2) = Ring{i}(:,Column_intensity); 
angle{i} = sortrows(angle{i},1);
[~, Anglesind{i}] = histc(angle{i}(:,1),anglerange);
angle{i}(:,3) = Anglesind{i};
end

%% Making distributions
distributions(:,1) = anglerange;
distributions_reshaped(:,1) = anglerange(1:(end-1));
for i=1:numel(files)
    for m = 1:length(anglerange)
        for k=1:length(angle{i}(:,2))
            if angle{i}(k,3) == m
                distributions(m,i+1) = distributions(m,i+1) + angle{i}(k,2);
            end
        end
    end
    distributions(:,i+1) = distributions(:,i+1)/sum(distributions(:,i+1));
end

%% Resaping distributions to start with maximum
for i=1:numel(files)
    [PeakPos{i}, PeakMag{i}] = peakfinder(distributions(:,i+1));
    [~, peakindex(i)]= max(PeakMag{i});
    shiftposition(i) = PeakPos{i}(peakindex(i));
    for s=0 : length(anglerange)-1
        if s<shiftposition(i)
            distributions_reshaped(180+s+1-shiftposition(i),i+1) = distributions(s+1,i+1);
        else
            distributions_reshaped(s+1-shiftposition(i),i+1) = distributions(s+1,i+1);
        end
    end
    %figure; bar(anglerange(1:(end-1)), distributions_reshaped(:,i+1));
    [PeakPos_reshaped{i}, PeakMag_reshaped{i}] = peakfinder(distributions_reshaped(:,i+1));
end

%% Summarized distribution and peaks distribution
distributions_reshaped(:,numel(files)+2) = sum(distributions_reshaped(:,2:(end-1)),2);
distributions_reshaped(:,numel(files)+2) = distributions_reshaped(:,numel(files)+2)/sum(distributions_reshaped(:,numel(files)+2));
%collecting all peaks
PeaksNumber = 0;
for i=1:numel(files)
    PeaksNumber = PeaksNumber + length(PeakPos_reshaped{i});
end
all_peaks = zeros(PeaksNumber,2);
counter1 = 0;
for i=1:numel(files)
    counter2 = 0;
    for k = 1:length(PeakPos_reshaped{i})
        counter2 = counter2 +1;
        counter1 = counter1 +1;
        all_peaks(counter1,1) = PeakPos_reshaped{i}(counter2,1)*bin_size;
        all_peaks(counter1,2) = PeakMag_reshaped{i}(counter2,1)*bin_size;
    end
end
all_peaks = sortrows(all_peaks,1);
%getting distribution of all peaks by intensity
angle_distribution = zeros(length(all_peaks_range)-1,2);
angle_distribution(:,1) = all_peaks_range(1:end-1);
for m = 2:length(all_peaks_range)
    for k=1:length(all_peaks)
        if all_peaks(k,1) <= m*binAnlges && all_peaks(k,1) > (m-1)*binAnlges
            angle_distribution(m,2) = angle_distribution(m,2) + all_peaks(k,2);
        end
    end
end
angle_distribution(:,2) = angle_distribution(:,2)/sum(angle_distribution(:,2));

%% Saving data
cd(resultdir);
headers = {'Angle', 'cluters with 2 removed', 'clusters all'};
csvwrite_with_headers('angles_clusters.csv', Data_angles, headers);
image5 = figure; 
bar(angle_distribution(:,1), angle_distribution(:,2));
print(image5, 'peaks_distribution.tif', '-dtiff', '-r150');
image6 = figure; 
bar(anglerange(1:(end-1)), distributions_reshaped(:,numel(files)+2));
print(image6, 'aligned_distribution.tif', '-dtiff', '-r150');
csvwrite('distributions_circular.csv',distributions_reshaped);
headers = {'Angle', 'Number of peaks'};
csvwrite_with_headers('peaks_distribution.csv', angle_distribution, headers);


%% getting all angles for each cluster from the center
angle_center = struct([]);
angle_cluster = zeros(1,1);
m=0;
for i=1:numel(files)
        angle_center{i}(:,1) = atan2d(Centroid{i}(:,2)-centerY(i), Centroid{i}(:,1)-centerX(i));
        angle_center{i}(:,2:3) = Centroid{i}(:,1:2);
        angle_center{i} = sortrows(angle_center{i},1);
        angle_center{i} = [angle_center{i}(end,:);angle_center{i}];
        for k=2:length(angle_center{i})
            m=m+1;
            u = [angle_center{i}(k,2)-centerX(i) angle_center{i}(k,3)-centerY(i) 0];
            v = [angle_center{i}(k-1,2)-centerX(i) angle_center{i}(k-1,3)-centerY(i) 0];
            angle_cluster(m) = atan2d(norm(cross(u,v)),dot(u,v));
        end
end

dist_angles=histcounts(angle_cluster,0:bin_size:max(angle_cluster))';
bincenter_cluster = (0:bin_size:max(angle_cluster)) + bin_size/2;
options_cluster = fitoptions('gauss2','Lower', [0 0 0 0 max(angle_cluster)/2 0],...
        'Upper', [Inf max(angle_cluster)/2 Inf Inf max(angle_cluster) Inf]);
[curve_cluster,gof_cluster] = fit(bincenter_cluster(1:end-1)',dist_angles,'gauss2', options_cluster);
image7 = figure;
plot(bincenter_cluster(1:end-1)', dist_angles', 'o',...
        bincenter_cluster(1:end-1)', curve_cluster(bincenter_cluster(1:end-1)'),'r');
title(num2str(gof_cluster.rsquare));
print(image6, 'aligned_distribution.tif', '-dtiff', '-r150');

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