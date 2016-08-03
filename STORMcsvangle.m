%% This script tests periodicity of clusters

%% Pairwise angles between clusters (clean)
Angles_all = zeros(1,1);
for i=1:numel(files)
    CounterAngles = 0;
    for k=1:(length(Centroid{i})-1)
        V1 = [Centroid{i}(k,1)-centerX(i) Centroid{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid{i})
            V2 = [Centroid{i}(m,1)-centerX(i) Centroid{i}(m,2)-centerY(i) 0];
            CounterAngles = CounterAngles +1;
            Angles{i}(CounterAngles,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
    Angles_all = [Angles_all; Angles{i}];
end

binrangeAngles = [0 : binAnlges : 180];
bincenterAngles=binrangeAngles(1:(end-1)) + bin_size/2;
[NAngles, binsAngles] = histc(Angles_all,binrangeAngles);

%% Pairwise angles between clusters (all)
Angles_all2 = zeros(1,1);
for i=1:numel(files)
    CounterAngles2 = 0;
    for k=1:(length(Centroid2{i})-1)
        V1 = [Centroid2{i}(k,1)-centerX(i) Centroid2{i}(k,2)-centerY(i) 0];
        for m=(k+1):length(Centroid2{i})
            V2 = [Centroid2{i}(m,1)-centerX(i) Centroid2{i}(m,2)-centerY(i) 0];
            CounterAngles2 = CounterAngles2 +1;
            Angles2{i}(CounterAngles2,1) = atan2d(norm(cross(V1,V2)),dot(V1,V2));
        end
    end
    Angles_all2 = [Angles_all2; Angles2{i}];
end

NAngles2 = histc(Angles_all2,binrangeAngles);

Data_angles = [binrangeAngles' NAngles NAngles2];
cd('STORMcsv/');
mkdir('result');
cd('result/');
csvwrite('angles.csv', Data_angles)
cd('../../');
[test, peakMag] = peakfinder(NAngles);

anglerange = [0 : bin_size : 360];
profile_all = zeros(180,1);
all_peaks = [0, 0];
for i=1:numel(files)
    for m=1:length(Ring{i}(:,3))
        if Ring{i}(m,3)>centerX(i)
            angle{i}(m)= atand((Ring{i}(m,3)-centerX(i))./(Ring{i}(m,4)-centerY(i)));
        else angle{i}(m)= 180 +atand((Ring{i}(m,3)-centerX(i))./(Ring{i}(m,4)-centerY(i)));
        end
    end
    for k=1:length(angle{i})
        if angle{i}(k)<0 angle{i}(k) = 270 - angle{i}(k);
        end
    end
    [Allangles{i}, Anglesind{i}] = histc(angle{i},anglerange);
    Allangles{i} = Allangles{i}(1:(end-1));
    %figure; bar(anglerange(1:(end-1)), Allangles{i});
    [PeakPos{i}, PeakMag{i}] = peakfinder(Allangles{i});
    [strongpeak(i) peakindex(i)]= max(PeakMag{i});
    shiftposition(i) = PeakPos{i}(peakindex(i));
    
    for s=0 : length( Allangles{i})-1
        if s<shiftposition(i)
             pr_reshaped{i}(180+s+1-shiftposition(i)) = Allangles{i}(s+1);
         else
             pr_reshaped{i}(s+1-shiftposition(i)) = Allangles{i}(s+1);
         end
     end
     %figure; bar(anglerange(1:(end-1)), pr_reshaped{i});
     [PeakPos2{i}, PeakMag2{i}] = peakfinder(pr_reshaped{i});
     all_peaks = [all_peaks, PeakPos2{i}*bin_size];
     
     profile_all = profile_all + pr_reshaped{i}';
end
all_peaks_range = [0 : binAnlges : 360];
all_peaks = all_peaks(all_peaks~=0);
[all_peaks_hist,all_peaks_ind] = histc(all_peaks,all_peaks_range);
figure; bar(all_peaks_range, all_peaks_hist);
figure; bar(anglerange(1:(end-1)), profile_all);
[PeakPos{numel(files)+1}, PeakMag{numel(files)+1}] = peakfinder(profile_all);