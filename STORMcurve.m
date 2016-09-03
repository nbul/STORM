%% Curve fitting

%% Assign memory
curve2 = struct([]);
gof2 = struct([]);
image2 = figure;

%% Fitting options and equation
fo = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', [curve{numel(files)+1}.a1 sigma1(numel(files)+1) curve{numel(files)+1}.b1]);
ft = fittype('2*pi*A*exp(-(x*x + R*R)/(2*sigma*sigma))*besseli(0, (R*x/(sigma*sigma)))',...
    'coefficients',{'A','sigma', 'R'},'independent','x','options',fo);

%% Fitting curve and plotting
for i=1:numel(files)
    [curve2{i},gof2{i}] = fit(binsnew'*bin_size,p{i},ft);
    subplot(4, ceil(numel(files)/4),i);
    plot(binsnew'*bin_size, p{i}, 'o', binsnew'*bin_size, curve2{i}(binsnew'*bin_size));
    title(num2str(gof2{i}.rsquare));
end

%% Fitting average curve and plotting
[curve2{numel(files)+1},gof2{numel(files)+1}] = fit(binsnew'*bin_size,p{numel(files)+1},ft);
subplot(4, ceil(numel(files)/4),numel(files)+1);
plot(binsnew'*bin_size, p{numel(files)+1}, 'o', binsnew'*bin_size, curve2{numel(files)+1}(binsnew'*bin_size));
title(num2str(gof2{numel(files)+1}.rsquare));

cd(resultdir);
print(image2, 'summarized_distributions_modified.tif', '-dtiff', '-r150');
cd(currdir);