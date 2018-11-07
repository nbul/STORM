clc;
close all;
clear variables;

currdir = pwd;
addpath(pwd);
filedir = uigetdir();

excludeangle = 30;
length = 12;
usedefault = questdlg(strcat('Use default settings: AB length = ', num2str(length),'?'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No')
    parameters = inputdlg({'Enter AB length:'},...
        'Parameters',1,{num2str(length)});
    % Redefine extension
    length = str2double(parameters{1});
end
cd(filedir);

Distance = csvread('180911.csv');
Distance = Distance(:);
Distance(Distance==0) = [];
Interval = zeros(size(Distance,1),2);
Position = zeros(181,size(Distance,1))+20;
for i=1:size(Distance,1)
    for k=0:1:180
        P1 = length*cosd(k)/2;
        all = 0:1:180;
        all = all(all < 180-(k-excludeangle) | all > 180-(k+excludeangle));
        P2 = length*cosd(k) + length*cosd(all)/2;
        P12 = P2 - P1;
        Ptrue = P12 < Distance(i) + 0.06 & P12 > Distance(i) - 0.06;
        if sum(Ptrue)>0
            Position(k+1,i) = -length*cosd(k)/2;
        end
    end
    Interval(i,1) = min(Position(Position(:,i)<20,i));
    Interval(i,2) = max(Position(Position(:,i)<20,i));
end

csvwrite('intervals.csv',Interval);

cd(currdir);
clear variables;
close all;
clc;