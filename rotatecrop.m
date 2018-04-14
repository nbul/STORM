coordinates_rot2 = struct([]);
coordinates_rot3 = struct([]);
coordinates_cr_all = zeros(1,3);
Center_cr = zeros(numel(files),2);
cd(filedir);
if exist([filedir,'/cropped'],'dir') == 0
    mkdir(filedir,'/cropped');
end
cr_im_dir = [filedir, '/cropped'];

if exist([filedir,'/croppedresult'],'dir') == 0
    mkdir(filedir,'/croppedresult');
end
cr_dir = [filedir, '/croppedresult'];


for i = 1:numel(files)
    coordinates_rot2{i} = Ring{i};
    R = [cosd(-rotation(i)) -sind(-rotation(i)); sind(-rotation(i)) cosd(-rotation(i))];
    Center_mat_border = repmat([Center(i,1); Center(i,2)], 1, length(coordinates_rot2{i}))';
    
    coordinates_rot2{i}(:,3:4) = coordinates_rot2{i}(:,3:4) - Center_mat_border;
    coordinates_rot2{i}(:,3:4) = (R* coordinates_rot2{i}(:,3:4)')';
    
     %% Cropping
    image13 = figure;
    h3 = scatter(coordinates_rot2{i}(:,3),coordinates_rot2{i}(:,4),10, [0 0 0], 'filled');
    axis equal;
    xlim([min(coordinates_rot2{i}(:,3))-50 max(coordinates_rot2{i}(:,3))+50]);
    ylim([min(coordinates_rot2{i}(:,4))-50 max(coordinates_rot2{i}(:,4))+50]);
    uiwait(msgbox('Select 4 points between which to crop','modal'));
    [x_user, y_user] = ginput(4); % Let the user select an x-value from which to crop.
    coordinates_rot2{i}(coordinates_rot2{i}(:,4) > max(y_user),:) = [];
    coordinates_rot2{i}(coordinates_rot2{i}(:,4) < min(y_user),:) = [];
    coordinates_rot2{i}(coordinates_rot2{i}(:,3) > max(x_user),:) = [];
    coordinates_rot2{i}(coordinates_rot2{i}(:,3) < min(x_user),:) = [];
    
    
    close all;
    cd(cr_im_dir);
    csvwrite([num2str(i),'_cropped.csv'], coordinates_rot2{i});
    %% Splitting proportianally to the signal
    t=1;
    for o=1:length(coordinates_rot2{i})
        m=ceil(coordinates_rot2{i}(o,18)/100);
        for f=1:m
            coordinates_rot3{i}(t,:) = coordinates_rot2{i}(o,:);
            t=t+1;
        end
    end
    %% recenter
    Center_cr(i,:) = (mean(coordinates_rot3{i}(:,3:4),1));
    Center_mat_cr = repmat([Center_cr(i,1); Center_cr(i,2)], 1, length(coordinates_rot3{i}))';    
    coordinates_rot3{i}(:,3:4) = coordinates_rot3{i}(:,3:4) - Center_mat_cr;
    %% Rotated image
    Nbins = 64;
    xmin = min(coordinates_rot3{i}(:,3))-50; xmax = max(coordinates_rot3{i}(:,3))+50;
    ymin = min(coordinates_rot3{i}(:,4))-50; ymax = max(coordinates_rot3{i}(:,4))+50;
    H = hist3([coordinates_rot3{i}(:,3) coordinates_rot3{i}(:,4)],...
        {linspace(xmin, xmax, Nbins), linspace(ymin, ymax, Nbins)});
    
    H = imadjust(double(H/max(H(:))));
    H = imgaussfilt(H,1);
    H = imadjust(H);
    imwrite(H.',[num2str(i),'_cropped.tif']);
    
    coordinates_cr_all = [coordinates_cr_all;coordinates_rot3{i}(:,[3:4,18])];
end

coordinates_cr_all(coordinates_cr_all(:,3) ==0,:) =[];


%% Individual distributions and distances
bin_size = 2;
dist_cr_length = max(max(abs(coordinates_cr_all(:,1:2))));
binrange_cr = 0 : bin_size : (dist_cr_length*2);
bincenter_cr=binrange_cr(1:(end-1)) + bin_size/2;
Dist_cr = zeros(length(bincenter_cr),numel(files)+1);
Dist2_cr = zeros(length(bincenter_cr),numel(files)+1);
curve_cr = struct([]);
gof_cr = struct([]);
curve1_cr = struct([]);
gof1_cr = struct([]);
curve2_cr = struct([]);
gof2_cr = struct([]);
image5 = figure;
Distance_cr = zeros(numel(files),7);
for i=1:numel(files)   
    Dist_cr(:,i)=histcounts(coordinates_rot3{i}(:,4),binrange_cr-dist_cr_length)';
    Dist_cr(:,i)=  Dist_cr(:,i) / sum( Dist_cr(:,i));
    Dist2_cr(:,i)=histcounts(coordinates_rot3{i}(:,3),binrange_cr-dist_cr_length)';
    Dist2_cr(:,i)=  Dist2_cr(:,i) / sum( Dist2_cr(:,i));
    options = fitoptions('gauss2','Lower', [max(Dist_cr(:,i))/5 0 15 max(Dist_cr(:,i))/5 dist_cr_length 15],...
        'Upper', [max(Dist_cr(:,i)) dist_cr_length dist_cr_length max(Dist_cr(:,i)) dist_cr_length*2 dist_cr_length]);
    [curve1_cr{i},gof1_cr{i}] = fit(bincenter_cr',Dist_cr(:,i),'gauss2',options);
    [curve2_cr{i},gof2_cr{i}] = fit(bincenter_cr',Dist_cr(end:-1:1,i),'gauss2',options);
    if gof1_cr{i}.rsquare>gof2_cr{i}.rsquare
        curve_cr{i} = curve1_cr{i};
        gof_cr{i} = gof1_cr{i};
    else
        curve_cr{i} = curve2_cr{i};
        gof_cr{i} = gof2_cr{i};
        Dist_cr(:,i) = Dist_cr(end:-1:1,i);
    end
    subplot(4, ceil(numel(files)/4),i);
    plot(bincenter_cr', Dist_cr(:,i), 'o',...
        bincenter_cr', curve_cr{i}(bincenter_cr'),'r');
    title(num2str(gof_cr{i}.rsquare));
    if gof_cr{i}.rsquare>cutoff
        Distance_cr(i,7) = 1;
    end 
    
    Distance_cr(i,1) = curve_cr{i}.b1;
    Distance_cr(i,2) = curve_cr{i}.c1;
    Distance_cr(i,3) = curve_cr{i}.b2;
    Distance_cr(i,4) = curve_cr{i}.c2;
    Distance_cr(i,5) = gof_cr{i}.rsquare;
    Distance_cr(i,6) = abs(curve_cr{i}.b1-curve_cr{i}.b2);
end

Filenames = 1:1:numel(files);
Distance2_cr = [Filenames' Distance_cr];
cd(cr_dir);
print(image5, 'Distributions.tif', '-dtiff', '-r150');

%% Building final image and summary distribution

final_lines_cr = zeros(length(bincenter_cr)+2,length(bincenter_cr)+2);
for n=1:length(coordinates_cr_all)
   final_lines_cr(ceil((coordinates_cr_all(n,2)+max(max(abs(coordinates_cr_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_cr_all(n,1)+max(max(abs(coordinates_cr_all(:,1:2)))))/bin_size)+1) = ...
       final_lines_cr(ceil((coordinates_cr_all(n,2)+max(max(abs(coordinates_cr_all(:,1:2)))))/bin_size)+1,...
       ceil((coordinates_cr_all(n,1)+max(max(abs(coordinates_cr_all(:,1:2)))))/bin_size)+1) + coordinates_cr_all(n,3);
end
image6 = figure;
imshow(final_lines_cr, [0, max(max(final_lines_cr))]);
print(image6, 'summed_image.tif', '-dtiff', '-r150');


for i=1:numel(files)
    if Distance_cr(i,7) == 1
        Dist_cr(:,numel(files)+1) = Dist_cr(:,numel(files)+1) + Dist_cr(:,i);
        Dist2_cr(:,numel(files)+1) = Dist2_cr(:,numel(files)+1) + Dist2_cr(:,i);
    end
end
Dist_cr(:,numel(files)+1)=  Dist_cr(:,numel(files)+1) / sum( Dist_cr(:,numel(files)+1));
Dist2_cr(:,numel(files)+1)=  Dist2_cr(:,numel(files)+1) / sum( Dist2_cr(:,numel(files)+1));

final_cr = mtimes(Dist_cr(:,numel(files)+1),Dist2_cr(:,numel(files)+1)');

print(image5, 'Distributions.tif', '-dtiff', '-r150');
image6 = figure;
imshow(final_cr, [0, max(max(final_cr))]);
print(image6, 'averaged_image.tif', '-dtiff', '-r150');

image7 = figure;
[curve1_cr{numel(files)+1},gof1_cr{numel(files)+1}] = fit(bincenter_cr',Dist_cr(:,numel(files)+1),'gauss2',options);
[curve2_cr{numel(files)+1},gof2_cr{numel(files)+1}] = fit(bincenter_cr',Dist_cr(end:-1:1,numel(files)+1),'gauss2',options);
if gof1_cr{numel(files)+1}.rsquare>=gof2{numel(files)+1}.rsquare
    curve_cr{numel(files)+1} = curve1_cr{numel(files)+1};
    gof_cr{numel(files)+1} = gof1_cr{numel(files)+1};
else
    curve_cr{numel(files)+1} = curve2_cr{numel(files)+1};
    gof_cr{numel(files)+1} = gof2_cr{numel(files)+1};
    Dist_cr(:,numel(files)+1) = Dist_cr(end:-1:1,numel(files)+1);
end
plot(bincenter_cr', Dist_cr(:,numel(files)+1), 'o',...
        bincenter_cr', curve_cr{numel(files)+1}(bincenter_cr'),'r');
title(num2str(gof_cr{numel(files)+1}.rsquare));
print(image7, 'average_distribution.tif', '-dtiff', '-r150');

Distance2_cr(numel(files)+1,2) = curve_cr{numel(files)+1}.b1;
Distance2_cr(numel(files)+1,3) = curve_cr{numel(files)+1}.c1;
Distance2_cr(numel(files)+1,4) = curve_cr{numel(files)+1}.b2;
Distance2_cr(numel(files)+1,5) = curve_cr{numel(files)+1}.c2;
Distance2_cr(numel(files)+1,6) = gof_cr{numel(files)+1}.rsquare;
Distance2_cr(numel(files)+1,7) = abs(curve_cr{numel(files)+1}.b1-curve_cr{numel(files)+1}.b2);
Otput_summary = 'Distances.csv';
headers2 = {'Ring', 'Center1', 'Width1', 'Center2', 'Width2', 'gof','Distance','Usage'};
csvwrite_with_headers(Otput_summary,Distance2_cr, headers2);

