clc;
clear variables;
close all;

currdir = pwd;
addpath(pwd);

res_dir = uigetdir();

if exist([res_dir,'/images'],'dir') == 0
    mkdir(res_dir,'/images');
end
im_dir = [res_dir, '/images'];

radius = [50, 91.5];

usedefault = questdlg(strcat('Use default settings: radius1 = ', num2str(radius(1)),...
    '; radius2 =',num2str(radius(2)),'?'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No')
    parameters = inputdlg({'Radius 1:', 'Radius 2:'},...
        'Parameters',1,{num2str(radius(1)), num2str(radius(2))});
    % Redefine extension
    radius(1) = str2double(parameters{1});
    radius(2) = str2double(parameters{2});
else
    parameters{1} = num2str(radius(1));
    parameters{2} = num2str(radius(2));
end

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
    
    if length(pks1)>1
        tail1 = [line(line>locs1(2))', norm1(line>locs1(2))'];
        tailtemp1 = flipud(tail1(2:end,:));
        tailtemp1(:,1) = tail1(1,1) - flipud(abs(tail1(2:end,1) - tail1(1,1)));
        tail1final = [tail1;tailtemp1];
        tail1final = sortrows(tail1final);
        tail1final(tail1final(:,2) < max(tail1final(:,2))/2,:) = [];
        fittail1 = fit(tail1final(:,1),tail1final(:,2),'poly2');
        r1 = roots([fittail1.p1 fittail1.p2 fittail1.p3]);
        m1 = (r1(1) + r1(2))/2;
        h1 = fittail1.p1*m1^2 + fittail1.p2*m1 + fittail1.p3;
        syms x;
        eqn1 = fittail1.p1*x^2 + fittail1.p2*x + fittail1.p3 == h1/2;
        solx1 = double(solve(eqn1,x));
        halfwidth1 = abs(solx1(1)-solx1(2));
    else
        halfwidth1 = 0;
        m1 = 0;
    end
    
    if length(pks2) > 1
    tail2 = [line(line>locs2(2))', norm2(line>locs2(2))'];
    tailtemp2 = flipud(tail2(2:end,:));
    tailtemp2(:,1) = tail2(1,1) - flipud(abs(tail2(2:end,1) - tail2(1,1)));
    tail2final = [tail2;tailtemp2];
    tail2final = sortrows(tail2final);
    tail2final(tail2final(:,2) < max(tail2final(:,2))/2,:) = [];
    fittail2 = fit(tail2final(:,1),tail2final(:,2),'poly2');
    r2 = roots([fittail2.p1 fittail2.p2 fittail2.p3]);
    m2 = (r2(1) + r2(2))/2;
    h2 = fittail2.p1*m2^2 + fittail2.p2*m2 + fittail2.p3;
    syms x;
    eqn2 = fittail2.p1*x^2 + fittail2.p2*x + fittail2.p3 == h2/2;
    solx2 = double(solve(eqn2,x));
    halfwidth2 = abs(solx2(1)-solx2(2));
    else
        m2 = 0;
        halfwidth2 = 0;
    end
    
    result(counter,:) = [l, radius(1), m1, halfwidth1, radius(2), m2,...
        halfwidth2, radius(2)-radius(1), m2-m1];
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

