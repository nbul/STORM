clc;
close all;
clear variables;

currdir = pwd;
addpath(pwd);
filedir = uigetdir();

length = 12;
usedefault = questdlg(strcat('Use default settings: AB length = ', num2str(length),'?'),'Settings','Yes','No','Yes');
if strcmp(usedefault, 'No')
    parameters = inputdlg({'Enter AB length:'},...
        'Parameters',1,{num2str(length)});
    % Redefine extension
    length = str2double(parameters{1});
end
cd(filedir);

Data = csvread('AB.csv');
Data(Data(:,1)==0,:) = [];
Distance = Data(:,1) - Data(:,2);
Interval = zeros(size(Data,1),2);
Position = zeros(181,size(Data,1))+20;
for i=1:size(Data,1)
    for k=0:1:180
        P1 = length*cosd(k)/2;
        P2 = length*cosd(k) + length*cosd(0:1:180)/2;
        P12 = P2 - P1;
        Ptrue = P12 < Distance(i) + 0.05 & P12 > Distance(i) - 0.05;
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