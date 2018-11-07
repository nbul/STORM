clc;
close all;
clear variables;

currdir = pwd;
addpath(pwd);

excludeangle = 30;
lengthAB = 12;
distance = 1;
usedefault = questdlg(strcat('Use default settings: 2-1 AB distance = ', num2str(distance),'?'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No')
    parameters = inputdlg({'Enter 2-1 AB distance:'},...
        'Parameters',1,{num2str(distance)});
    % Redefine extension
    distance = str2double(parameters{1});
end

Interval = zeros(1,2);
Position = zeros(181,1)+20;
angle = zeros(1,1);
counter = 0;
for k=0:1:180
    P1 = lengthAB*cosd(k)/2;
    all = 0:1:180;
    all = all(all < 180-(k-excludeangle) | all > 180-(k+excludeangle));
    P2 = lengthAB*cosd(k) + lengthAB*cosd(all)/2;
    P12 = P2 - P1;
    Ptrue = P12 < distance + 0.06 & P12 > distance - 0.06;
    
    if sum(Ptrue)>0
        counter = counter + 1;
        Position(k+1) = -lengthAB*cosd(k)/2;
        angle(counter) = k;
    end
end
Interval(1) = min(Position(Position(:,1)<20,1));
Interval(2) = max(Position(Position(:,1)<20,1));
angle2 = [sort(0-angle), angle]';
angle3 = angle2;

if min(angle2) == -180
    angle3(angle2<0) = angle2(angle2<0) + 360;
    angle3 = sort(angle3);
end

D2 = 2*(distance - lengthAB*cosd(angle3)/2)/lengthAB;
D2(D2>1 | D2<-1) = [];
angles2d = acosd(D2);
angles2d = sort(angles2d);

angle2d2 = [sort(0-angles2d); angles2d];
angle2d3 = angle2d2;

if min(acosd(D2)) > (180 - max(acosd(D2)))
    angle2d3(angle2d2<0) = angle2d2(angle2d2<0) + 360;
    angle2d3 = sort(angle2d3);
end


subplot('Position', [0.1,0.5, 0.4,0.4]);
t = linspace(min(angle3),max(angle3));
t2 = linspace(max(angle3)-360, min(angle3));
t3 = linspace(-180,180);
h1 = patch([0,lengthAB*cosd(t),0],[0,lengthAB*sind(t),0],'r', 'FaceAlpha', 0.3, 'EdgeColor', 'k','LineWidth',3,'LineStyle', ':');
hold on;
h2 = patch([0,lengthAB*cosd(t2),0],[0,lengthAB*sind(t2),0],'w','LineWidth',3,'EdgeColor', 'w', 'LineStyle', '-');
hold on;
h3 = patch([0,0.5*cosd(t3),0],[0,0.5*sind(t3),0], 'b', 'EdgeColor', 'b', 'LineStyle', '-');
axis equal


subplot('Position', [0.5,0.5, 0.4,0.4]);
t4 = linspace(min(angle2d3),max(angle2d3));
t5 = linspace(max(angle2d3)-360, min(angle2d3));
h4 = patch([0,lengthAB*cosd(t4),0],[0,lengthAB*sind(t4),0],'b', 'FaceAlpha', 0.3, 'EdgeColor', 'k','LineWidth',3,'LineStyle', ':');
hold on;
h5 = patch([0,lengthAB*cosd(t5),0],[0,lengthAB*sind(t5),0],'w','LineWidth',3,'EdgeColor', 'w', 'LineStyle', '-');
hold on;
axis equal
hold off

hL = subplot('Position', [0.1,0.35, 0.8,0.05]);
poshL = get(hL,'position'); 
lgd = legend(hL,[h1, h4, h3, h2], {'Allowed angles 1st','Allowed angles 2d','antigen',strcat('Distance = ',num2str(distance))},...
    'Location', 'eastoutside','Orientation','horisontal');
lgd.FontSize = 14;
set(lgd,'position',poshL);
axis(hL,'off');      