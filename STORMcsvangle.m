%% This script tests periodicity of clusters


%% getting all angles for each cluster from the center
angle_center = struct([]);
angle_cluster = zeros(1,1);
m=0;
for i=1:numel(files)
    if usage(i)>0
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
end

%bin_size2 = 5;
dist_angles=histcounts(angle_cluster,0:bin_size:max(angle_cluster))';
bincenter_cluster = (0:bin_size:max(angle_cluster)) + bin_size/2;
options_cluster = fitoptions('gauss2','Lower', [0 0 0 0 max(angle_cluster)/2 0],...
        'Upper', [Inf max(angle_cluster)/2 Inf Inf max(angle_cluster) Inf]);
[curve_cluster,gof_cluster] = fit(bincenter_cluster(1:end-1)',dist_angles,'gauss2',options_cluster);
[curve_cluster2,gof_cluster2] = fit(bincenter_cluster(1:end-1)',dist_angles,'gauss1');
conf_cluster = confint(curve_cluster);
conf_cluster2 = confint(curve_cluster2);

Result_angle = zeros(7,5);
Result_angle(1,1) = curve_cluster.b1;
Result_angle(2,1) = conf_cluster(1,2);
Result_angle(3,1) = conf_cluster(2,2);
Result_angle(5,1) = curve_cluster2.b1;
Result_angle(6,1) = conf_cluster2(1,2);
Result_angle(7,1) = conf_cluster2(2,2);
Result_angle(1,2) = curve_cluster.c1;
Result_angle(2,2) = conf_cluster(1,3);
Result_angle(3,2) = conf_cluster(2,3);
Result_angle(5,2) = curve_cluster2.c1;
Result_angle(6,2) = conf_cluster2(1,3);
Result_angle(7,2) = conf_cluster2(2,3);
Result_angle(1,3) = curve_cluster.b2;
Result_angle(2,3) = conf_cluster(1,5);
Result_angle(3,3) = conf_cluster(2,5);
Result_angle(1,4) = curve_cluster.c2;
Result_angle(2,4) = conf_cluster(1,6);
Result_angle(3,4) = conf_cluster(2,6);
Result_angle(1,5) = gof_cluster.rsquare;
Result_angle(5,5) = gof_cluster2.rsquare;

cd(resultdir);
image5 = figure;
subplot(1,2,1);
plot(bincenter_cluster(1:end-1)', dist_angles', 'o',...
        bincenter_cluster(1:end-1)', curve_cluster(bincenter_cluster(1:end-1)'),'r');
title(num2str(gof_cluster.rsquare));
hold on;
subplot(1,2,2);
plot(bincenter_cluster(1:end-1)', dist_angles', 'o',...
        bincenter_cluster(1:end-1)', curve_cluster2(bincenter_cluster(1:end-1)'),'r');
title(num2str(gof_cluster2.rsquare));
hold off;
print(image5, 'Angle_between_clusters.tif', '-dtiff', '-r150');

Otput_summary = 'Angle_between_clusters.csv';
headers2 = {'Center1', 'Width1', 'Center2', 'Width2', 'gof'};
csvwrite_with_headers(Otput_summary,Result_angle, headers2);
cd(currdir);