%% This script finds the centers of clusters and rings

%% Assign memory

silh = zeros(1, numel(files));
ClusterNumber = zeros(1, numel(files));
centerX = zeros(1, numel(files));
centerY = zeros(1, numel(files));
Rfit = zeros(1, numel(files));
RemoveIdx1 = zeros(2, numel(files));

Ring = struct([]);
Cluster = struct([]);
IntensityCluster = struct([]);
Centroid = struct([]);
Centroid2 = struct([]);
coordinates = struct([]);
coordinates2 = struct([]);
distances = struct([]);

cd(filedir);
%% Reading files and all coordinates
for i=1:numel(files)
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4)];
end

%% Determining optimum number of clusters and clustering
for i=1:numel(files)
    clear test_k silh;
    for q=6:15
        test_k = kmeans(coordinates{i},q);
        silh(q-5) = mean(silhouette(coordinates{i},test_k));
    end
    [num1, idx1] = max(silh);
    ClusterNumber(i) = idx1+7;
    [Cluster{i}, Centroid{i}]= kmeans(coordinates{i},ClusterNumber(i));
    Centroid2{i}= Centroid{i};
    coordinates2{i} = [Cluster{i} coordinates{i}];
end

%% Identifying cluster intensities
for i=1:numel(files)
    IntensityCluster{i}=zeros(ClusterNumber(i),1);
    for k=1:ClusterNumber(i)
        IntensityCluster{i}(k) = sum(Ring{i}(Cluster{i}(:,1)==k,Column_intensity));
    end
end

cd(currdir);

%% Removing two clusters with least intensity
for i=1:numel(files)
    clear IntensityCluster2;
    ClusterNames = 1:1:ClusterNumber(i);
    Centroid{i} = [Centroid{i} ClusterNames'];
    [~, RemoveIdx1(1,i)] = min(IntensityCluster{i});
    IntensityCluster2=IntensityCluster{i};
    IntensityCluster2(RemoveIdx1(1,i),:) = [];
    Centroid{i}(RemoveIdx1(1,i),:) = [];
    coordinates2{i}(any(coordinates2{i}==RemoveIdx1(1,i),2),:) = [];
    [~, Idxtemp] = min(IntensityCluster2);
    RemoveIdx1(2,i) = Centroid{i}(Idxtemp,3);
    Centroid{i}(Idxtemp,:) = [];
    coordinates2{i}(any(coordinates2{i}==RemoveIdx1(2,i),2),:) = [];
end

%% Getting center and distances of all coordinates from the center
for i=1:numel(files)
    [centerX(i),centerY(i),Rfit(i)] = circfit(Centroid{i}(:,1),Centroid{i}(:,2));
    distances{i} = sqrt((Ring{i}(:,3)-centerX(i)).*(Ring{i}(:,3)-centerX(i)) + (Ring{i}(:,4)-centerY(i)).*(Ring{i}(:,4)-centerY(i)));
end