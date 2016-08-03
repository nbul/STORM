clear all
close all

bin_size = 2;
binAnlges = 4;
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

%% Pairwise angles between clusters (clean)
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

%% Pairwise angles between clusters (all)
Angles_all2 = zeros(1,1);
for i=1:numel(files)
    CounterAngles2 = 0;
    for k=1:(length(Centroid2{i})-1)
        V1 = [Centroid2{i}(k,1)-centerX(i) Centroid2{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid2{i})
            V2 = [Centroid2{i}(m,1)-centerX(i) Centroid2{i}(m,2)-centerY(i) 0];
            CounterAngles2 = CounterAngles2 +1;
            Angles2{i}(CounterAngles2,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
    Angles_all2 = [Angles_all2; Angles2{i}];
end

NAngles2 = histc(Angles_all2,binrangeAngles);

Data_angles = [binrangeAngles' NAngles NAngles2];
cd('STORMcsv/');
mkdir('result');
cd('result/');
csvwrite('angles.csv', Data_angles)
cd('../../');