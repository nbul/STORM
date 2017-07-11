%% This script finds the centers of clusters and rings

%% Assign memory

silh = zeros(1, numel(files));
ClusterNumber = zeros(1, numel(files));
centerX = zeros(1, numel(files));
centerY = zeros(1, numel(files));
Rfit = zeros(1, numel(files));
R_square_2 = zeros(1, numel(files));
R_square_all = zeros(1, numel(files));
RemoveIdx1 = zeros(2, numel(files));
usage = zeros(1,1);

Ring = struct([]);
Cluster = struct([]);
IntensityCluster = struct([]);
Centroid = struct([]);
Centroid2 = struct([]);
coordinates = struct([]);
coordinates2 = struct([]);
distances = struct([]);


Sigma = 'full'; % possible to change to 'diagonal'
nSigma = numel(Sigma);
SharedCovariance = true; % possible to change to 'true'
threshold = sqrt(chi2inv(0.99,2));
options = statset('MaxIter',1000);

cd(filedir);
d = 500;
for i=1:numel(files)
    %% Reading files and all coordinates
    Number1 = [num2str(i),'.csv'];
    Ring{i} = csvread(Number1, 1, 1);
   
    %% Determining optimum number of clusters and clustering
    l = 0;
    while l==0
        clear test_k silh;
        coordinates{i} = [Ring{i}(:,3) Ring{i}(:,4)];
        for q=6:15
            gmfit = fitgmdist(coordinates{i}(:,1:2),q,'CovarianceType',Sigma,...
            'SharedCovariance',SharedCovariance,'Options',options);
            AIC(q-5)= gmfit.AIC;
            %             test_k = kmeans(coordinates{i},q);
            %             silh(q-5) = mean(silhouette(coordinates{i},test_k));
        end
        [num1, idx1] = min(AIC);
        ClusterNumber(i) = idx1+7;
        x1 = linspace(min(coordinates{i}(:,1)) - 2,max(coordinates{i}(:,1)) + 2,d);
        x2 = linspace(min(coordinates{i}(:,2)) - 2,max(coordinates{i}(:,2)) + 2,d);
        [x1grid,x2grid] = meshgrid(x1,x2);
        X0 = [x1grid(:) x2grid(:)];
        
        clusterX = cluster(gmfit,coordinates{i}(:,1:2));
        mahalDist = mahal(gmfit,X0);
        coordinates2{i} = [coordinates{i} clusterX];
        h1 = gscatter(coordinates{i}(:,1),coordinates{i}(:,2),clusterX);
        hold on;
        for m = 1:ClusterNumber(i)
            idx = mahalDist(:,m)<=threshold;
            Color = h1(m).Color*0.75 + -0.5*(h1(m).Color - 1);
            h2 = plot(X0(idx,1),X0(idx,2),'.','Color',Color,'MarkerSize',1);
            uistack(h2,'bottom');
        end
        plot(gmfit.mu(:,1),gmfit.mu(:,2),'kx','LineWidth',2,'MarkerSize',10);
        
        
%         [Cluster{i}, Centroid{i}]= kmeans(coordinates{i},ClusterNumber(i));
        Centroid2{i}= Centroid{i};
        coordinates2{i} = [Cluster{i} coordinates{i}];
        coordinates{i} = [Cluster{i} coordinates{i}];
        
        
        %% Identifying cluster intensities
        
        IntensityCluster{i}=zeros(ClusterNumber(i),1);
        for k=1:ClusterNumber(i)
            IntensityCluster{i}(k) = sum(Ring{i}(Cluster{i}(:,1)==k,Column_intensity));
        end
        
        
        
        %% Removing two clusters with least intensity
        
        clear IntensityCluster2;
        ClusterNames = 1:1:ClusterNumber(i);
        Centroid{i} = [Centroid{i} ClusterNames'];
        Centroid2{i} = [Centroid2{i} ClusterNames'];
        [~, RemoveIdx1(1,i)] = min(IntensityCluster{i});
        IntensityCluster2=IntensityCluster{i};
        IntensityCluster2(RemoveIdx1(1,i),:) = [];
        Centroid{i}(RemoveIdx1(1,i),:) = [];
        coordinates2{i}(any(coordinates2{i}==RemoveIdx1(1,i),2),:) = [];
        [~, Idxtemp] = min(IntensityCluster2);
        RemoveIdx1(2,i) = Centroid{i}(Idxtemp,3);
        Centroid{i}(Idxtemp,:) = [];
        coordinates2{i}(any(coordinates2{i}==RemoveIdx1(2,i),2),:) = [];
        
        
        %% Getting center and distances of all coordinates from the center
        
        [centerX(i),centerY(i),Rfit(i)] = circfit(Centroid{i}(:,1),Centroid{i}(:,2));
        distances{i} = sqrt((Ring{i}(:,3)-centerX(i)).*(Ring{i}(:,3)-centerX(i)) + (Ring{i}(:,4)-centerY(i)).*(Ring{i}(:,4)-centerY(i)));
        
        
        %% Getting GOF (R^2) of circle fit
        SStot_2 = 0;
        SSres_2 = 0;
        SStot_all = 0;
        SSres_all = 0;
        for k=1:length(Centroid{i})
            SStot_2 = SStot_2 + (Centroid{i}(k,1) - centerX(i)) * (Centroid{i}(k,1) - centerX(i))...
                + (Centroid{i}(k,2) - centerY(i)) * (Centroid{i}(k,2) - centerY(i));
            SSres_2 = SSres_2 + (sqrt((Centroid{i}(k,1) - centerX(i)) * (Centroid{i}(k,1) - centerX(i))...
                + (Centroid{i}(k,2) - centerY(i)) * (Centroid{i}(k,2) - centerY(i)))-Rfit(i))*...
                (sqrt((Centroid{i}(k,1) - centerX(i)) * (Centroid{i}(k,1) - centerX(i))...
                + (Centroid{i}(k,2) - centerY(i)) * (Centroid{i}(k,2) - centerY(i)))-Rfit(i));
        end
        
        for k=1:length(Centroid2{i})
            SStot_all = SStot_all + (Centroid2{i}(k,1) - centerX(i)) * (Centroid2{i}(k,1) - centerX(i))...
                + (Centroid2{i}(k,2) - centerY(i)) * (Centroid2{i}(k,2) - centerY(i));
            SSres_all = SSres_all + (sqrt((Centroid2{i}(k,1) - centerX(i)) * (Centroid2{i}(k,1) - centerX(i))...
                + (Centroid2{i}(k,2) - centerY(i)) * (Centroid2{i}(k,2) - centerY(i)))-Rfit(i))*...
                (sqrt((Centroid2{i}(k,1) - centerX(i)) * (Centroid2{i}(k,1) - centerX(i))...
                + (Centroid2{i}(k,2) - centerY(i)) * (Centroid2{i}(k,2) - centerY(i)))-Rfit(i));
        end
        
        R_square_2(i) = 1-SSres_2/SStot_2;
        R_square_all(i) = 1-SSres_all/SStot_all;
        
        image1 = figure;
        subplot(1,2,1);
        for k=1:ClusterNumber(i)
            plot(coordinates{i}(coordinates{i}(:,1)==k,2),coordinates{i}(coordinates{i}(:,1)==k,3),'*')
            hold on
            plot(Centroid2{i}(Centroid2{i}(:,3)==k,1),Centroid2{i}(Centroid2{i}(:,3)==k,2),'k+','MarkerSize', 15)
            hold on
            
        end
        rectangle('position',[centerX(i)-Rfit(i),centerY(i)-Rfit(i),Rfit(i)*2,Rfit(i)*2],...
            'curvature',[1,1],'linestyle','-','edgecolor','r');
        hold on
        title(num2str(R_square_all(i)));
        hold off
        subplot(1,2,2);
        for k=1:ClusterNumber(i)
            if (k~=RemoveIdx1(1,i) && k~=RemoveIdx1(2,i))
                plot(coordinates{i}(coordinates{i}(:,1)==k,2),coordinates{i}(coordinates{i}(:,1)==k,3),'*')
                hold on
                plot(Centroid2{i}(Centroid2{i}(:,3)==k,1),Centroid2{i}(Centroid2{i}(:,3)==k,2),'k+','MarkerSize', 15)
                hold on
                rectangle('position',[centerX(i)-Rfit(i),centerY(i)-Rfit(i),Rfit(i)*2,Rfit(i)*2],...
                    'curvature',[1,1],'linestyle','-','edgecolor','r');
                hold on
                
            end
        end
        title(num2str(R_square_2(i)));
        hold off
        
        usedefault = questdlg(strcat('Are you happy about clustering)'),'Settings','Yes','No','Yes');
        if strcmp(usedefault, 'Yes');
            l=1;
            cd(clustdir);
            print(image1, [num2str(i),'_cluster.tif'], '-dtiff', '-r150');
            cd(filedir);
        end
        
        if R_square_2(i)>cutoff2
            usage(i) = 1;
        end
        close all;
    end
end

Filenames = 1:1:numel(files);
Rsquare = [Filenames' R_square_all' R_square_2'];
cd(clustdir);
Otput_summary = 'GOF_circle.csv';
headers2 = {'Ring', 'R_square all clusters', 'R_square - 2 clusters'};
csvwrite_with_headers(Otput_summary,Rsquare, headers2);
cd(currdir); 

