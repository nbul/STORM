clc;
clear variables;
close all;

currdir = pwd;
addpath(pwd);

res_dir = '/Volumes/DataGurdon/Khodor Hazime/STORMmodel';
im_dir = '/Volumes/DataGurdon/Khodor Hazime/STORMmodel/images';

radius = [92, 183];

angles = 0:40:320;

spread = 5:1:50;

counter = 0;
result = zeros(length(spread),9);
cd(im_dir);
for l=5:1:50
    counter = counter + 1;
    line = -radius(2)-50:0.01:radius(2)+50;
    norm1 = zeros(1,length(line));
    norm2 = zeros(1,length(line));
    for i=0:0.1:39.9
        
        x1=radius(1)*cosd(angles+i);
        x2=radius(2)*cosd(angles+i);
        
        for z=1:1:9
            normtemp1 = normpdf(line, x1(z), l);
            norm1 = norm1 + normtemp1;
            normtemp2 = normpdf(line, x2(z), l);
            norm2 = norm2 + normtemp2;
        end
        
        

    end
    norm1 = norm1/sum(norm1);
    norm2 = norm2/sum(norm2);
    
    [pks1,locs1] = findpeaks(norm1',line');
    [pks2,locs2] = findpeaks(norm2',line');
    
    tail1 = [line(line>locs1(2))', norm1(line>locs1(2))'];
    tailtemp1 = flipud(tail1);
    tailtemp1(:,1) = tail1(:,1) - tail1(end,1) + tail1(1,1);
    tail1final = [tail1;tailtemp1];
    tail1final(tail1final(:,2) < max(tail1final(:,2))/2,:) = [];
    halfwidth1 = max(tail1final(:,1)) - min(tail1final(:,1));
    
    tail2 = [line(line>locs2(2))', norm2(line>locs2(2))'];
    tailtemp2 = flipud(tail2);
    tailtemp2(:,1) = tail2(:,1) - tail2(end,1) + tail2(1,1);
    tail2final = [tail2;tailtemp2];
    tail2final(tail2final(:,2) < max(tail2final(:,2))/2,:) = [];
    halfwidth2 = max(tail2final(:,1)) - min(tail2final(:,1));
    
    result(counter,:) = [l, radius(1), locs1(2), halfwidth1, radius(2), locs2(2),...
        halfwidth2, radius(2)-radius(1), locs2(2)-locs1(2)];
    image = figure;
    plot(line', norm2'/max(norm2), 'r', line', norm1'/max(norm1), 'b', 'LineWidth', 2);
    ax = gca;
    ax.YAxis.Visible = 'off';
    ax.FontSize = 16;
    ax.XAxis.Label.String = 'Distance from center (nm)';
    ax.XAxis.Label.FontSize = 20;
    legend('IC138', 'RSP3');
    
    name = [num2str(l),'.tif'];
    print(image, name, '-dtiff', '-r150');
end

close all;
headers = {'precision', 'distance 1', 'peak position 1', 'peak 1 width at 50%', ...
    'distance 2', 'peak position 2', 'peak 2 width at 50%', 'theoretical difference', 'measured difference'};
cd(res_dir);
csvwrite_with_headers('model.csv',result,headers);

image = figure;
plot(result(:,1),result(:,6),'r',result(:,1),result(:,3),'b', 'LineWidth', 2);
ax = gca;
ax.FontSize = 16;
ax.XLim = [0 55];
ax.YLim = [0 200];
ax.XAxis.Label.String = 'Spread (nm)';
ax.XAxis.Label.FontSize = 20;
ax.YAxis.Label.String = 'Peak distance from center (nm)';
ax.YAxis.Label.FontSize = 20;
legend('IC138', 'RSP3');
print(image, 'PeakDistance_vs_spread', '-dtiff', '-r150');

image2 = figure;
plot(result(:,1),result(:,7),'r',result(:,1),result(:,4),'b', 'LineWidth', 2);
ax = gca;
ax.FontSize = 16;
ax.XLim = [0 55];
ax.YLim = [0 200];
ax.XAxis.Label.String = 'Spread (nm)';
ax.XAxis.Label.FontSize = 20;
ax.YAxis.Label.String = 'Peak width (nm)';
ax.YAxis.Label.FontSize = 20;
legend('IC138', 'RSP3');
print(image2, 'PeakWidth_vs_spread', '-dtiff', '-r150');

image3 = figure;
plot(result(:,1),result(:,8),':r',result(:,1),result(:,9),'k', 'LineWidth', 2);
ax = gca;
ax.FontSize = 16;
ax.XLim = [0 55];
ax.YLim = [0 200];
ax.XAxis.Label.String = 'Spread (nm)';
ax.XAxis.Label.FontSize = 20;
ax.YAxis.Label.String = 'Distance between two proteins (nm)';
ax.YAxis.Label.FontSize = 20;
legend('theoretical distance', 'distance between peaks');
print(image3, 'ProteinDistance_vs_spread', '-dtiff', '-r150');


cd(currdir);
close all;
clear variables;
clc;

