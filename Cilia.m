clc
clear variables
close all

Column_intensity = 18;

currdir = pwd;
addpath(pwd);
filedir = uigetdir();
files = dir(strcat(filedir,'/*.csv'));
cd(filedir);
mkdir(filedir,'/results');
resultdir = [filedir, '/results'];

Signal = struct([]);
coordinates = struct([]);
coordinates2 = struct([]);
Center = zeros(numel(files),2);
rotation = zeros(numel(files),1);

for i = 1:numel(files)
    %% Reading data
    Number1 = [num2str(i),'.csv'];
    Signal{i} = csvread(Number1, 1, 1);
    coordinates{i} = [Signal{i}(:,3) Signal{i}(:,4) Signal{i}(:,Column_intensity)];
    
    %% Splitting proportianally to the signal
    t = 1;
    for o=1:length(Signal{i})
        m=ceil(Signal{i}(o,18)/100);
        for f=1:m
            coordinates2{i}(t,:) = coordinates{i}(o,:);
            t=t+1;
        end
    end
    
    %% Determining rotation angle and shift
        mu = mean(coordinates2{i}(:,1:2),1);
        X_minus_mu = coordinates2{i}(:,1:2)-...
            repmat(mu, size(coordinates2{i}(:, 1)));
        Sigma2=(X_minus_mu'*X_minus_mu)/size(coordinates2{i}(:,1:2),1);
        [V, D]=eig(Sigma2);
        rotation(i) = atan2d(V(2,2), V(2,1));
        Center(i,:) = (mean(coordinates2{i}(:,1:2),1) +...
            mean(coordinates2{i}(:,1:2),1))/2;
end

%% Moving, rotating and summing up
coordinates_rot = struct([]);
coordinates_all = zeros(1,3);

for i=1:numel(files)
    
    R = [cosd(-rotation(i)) -sind(-rotation(i)); sind(-rotation(i)) cosd(-rotation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(coordinates2{i}))';
    
    coordinates_rot{i} = coordinates2{i}(:,1:2) - Center_mat_border;
    coordinates_rot{i} = (R* coordinates_rot{i}')';
    coordinates_rot{i}(:,3) = coordinates2{i}(:,3);
    coordinates_all = [coordinates_all;coordinates_rot{i}(:,1:3)]; %#ok<AGROW>
end
coordinates_all(coordinates_all(:,3) ==0,:) =[];

%% Projection and distribution

bin_size = 8;
dist_length = max(abs(coordinates_all(:,2)));
binrange = -dist_length : bin_size : dist_length;
bincenter=binrange(1:(end-1)) + bin_size/2;

Dist = zeros(length(bincenter),numel(files)+1);
for i=1:numel(files)
    Dist(:,i)=histcounts(coordinates_rot{i}(:,2),binrange)';
    Dist(:,i)=  Dist(:,i) / sum( Dist(:,i));
end

Dist(:,numel(files)+1) = histcounts(coordinates_all,binrange)';
Dist(:,numel(files)+1)=  Dist(:,numel(files)+1) / sum( Dist(:,numel(files)+1));

for i = 2:length(bincenter)-1
    Dist(i,numel(files)+1) = (Dist(i,numel(files)+1) + Dist(i-1,numel(files)+1) +...
        Dist(i+1,numel(files)+1))/3;
end
Dist(1,:) = [];
Dist(end,:) = [];
Dist(:,numel(files)+1) = (Dist(:,numel(files)+1) + flipud(Dist(:,numel(files)+1)))/2;
image = figure;
plot(bincenter(2:end-1)',Dist(:,numel(files)+1), [-mean(abs(locs)), mean(abs(locs))], pks, 'o', 'LineWidth', 2);
ax = gca;
ax.YAxis.Visible = 'off';
ax.FontSize = 16;
ax.XAxis.Label.String = 'Distance from center (nm)';
ax.XAxis.Label.FontSize = 20;

Image = figure;
plot(bincenter(2:end-1)',Dist(:,numel(files)+1));

[pks,locs, w, p] = findpeaks(Dist(:,numel(files)+1),bincenter(2:end-1)', 'MinPeakProminence', 0.005);

line = bincenter(2:end-1)';
norm1 = Dist(:,numel(files)+1);
tail = [line(line>=locs(2)), norm1(line>=locs(2))];
tailtemp = flipud(tail(2:end,:));
tailtemp(:,1) = tail(1,1) - flipud(abs(tail(2:end,1) - tail(1,1)));
tailfinal = [tail;tailtemp];
tailfinal = sortrows(tailfinal);
tailfinal(tailfinal(:,2) < max(tailfinal(:,2))/2,:) = [];
fittail = fit(tailfinal(:,1),tailfinal(:,2),'poly2');
r = roots([fittail.p1 fittail.p2 fittail.p3]);
m = (r(1) + r(2))/2;
h = fittail.p1*m^2 + fittail.p2*m + fittail.p3;
syms x;
eqn = fittail.p1*x^2 + fittail.p2*x + fittail.p3 == h/2;
solx = double(solve(eqn,x));
halfwidth = abs(solx(1)-solx(2));



