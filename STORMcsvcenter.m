%% This script finds the centers of clusters and rings

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