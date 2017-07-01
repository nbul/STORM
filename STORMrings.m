%% Save means and SD and make rings

peaks2 = zeros(numel(files)+1,6);
Number2 = 1:(numel(files)+1);
ClusterNumber = [ClusterNumber'; 0];
peaks = [Number2' radius' sigma1' pvalue' ClusterNumber usage'];
for i=1:(numel(files)+1)
peaks2(i,:) = [i curve2{i}.R curve2{i}.sigma gof2{i}.rsquare ClusterNumber(i) usage(i)];
end

cd(resultdir);
headers = {'ring', 'radius', 'width', 'p-value', 'Number of clusters','Used in sum?'};
csvwrite_with_headers('radiuses.csv', peaks, headers);
csvwrite_with_headers('radiuses_curve.csv', peaks2, headers);
cd(currdir);


%% Making summarized ring
% creating empty image and specifying its center
final_ring = zeros(length(bincenter)*2,length(bincenter)*2);
final_ring_fitted = zeros(length(bincenter)*2,length(bincenter)*2);
[columnsInImage, rowsInImage] = meshgrid(1:(length(bincenter)*2), 1:(length(bincenter)*2));
centerX2 = length(bincenter);
centerY2 = length(bincenter);

%Drowing the ring
sum_intensity = [sum_intensity; zeros(ceil(sqrt(2*length(bincenter)*length(bincenter))-length(sum_intensity))+1, 1)];
for l=1:(length(bincenter)*2)
    for n = 1:(length(bincenter)*2)
            final_ring(l,n) = abs(sum_intensity(ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)/max(sum_intensity));
            final_ring_fitted(l,n) = abs(curve2{numel(files)+1}((ceil(sqrt((l-centerX2)*(l-centerX2)+(n-centerY2)*(n-centerY2)))+1)*bin_size)...
                /max(curve2{numel(files)+1}(binsnew'*bin_size)));
    end
end

%Saving image
image3 = figure;
imshow(final_ring);
image4 = figure;
imshow(final_ring_fitted);
cd(resultdir);
print(image3, 'summarized_ring.tif', '-dtiff', '-r150');
print(image4, 'summarized_ring_fitted.tif', '-dtiff', '-r150');
cd(currdir);