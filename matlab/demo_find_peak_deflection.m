%% load data and metadata from binary file. change the file path if necessary
[meta, data] = read_pillar_bin('I:\mingxi\casminus_1.3um_2dhex_cell03_R3D\casminus_1.3um_2dhex_cell03_R3D.bin');
pixel_size = meta.pixel_size;

%% test one specified pillar, change the value if neccessary
pillar_id = 2547;

%% compute the deflections, use the Euclidean distance
dx = data.DX(:, pillar_id);
dy = data.DY(:, pillar_id);
d2 = sqrt(dx.^2 + dy.^2);

%% convert the pixel unit to nm by multiplying the pixel_size
d = pixel_size*d2;
f = 1 : meta.nframes;
% plot(f, d);
% title(['pillar ID=',num2str(pillar_id)])
% xlabel('frame ID')
% ylabel('deflection (nm)')

%% donoising
dd = sgolayfilt(d,3,11);
subplot(211),plot(f, d);
title(['Original signal, pillar ID=',num2str(pillar_id)])
xlabel('frame ID')
ylabel('deflection (nm)')
axis tight
subplot(212),plot(f, dd);
title(['De-noised signal,pillar ID=',num2str(pillar_id)])
xlabel('frame ID')
ylabel('deflection (nm)')

%% find peaks
peak_threshold = 200;
[max_d, max_frame_id] = max(dd);
min_peak_proms = min(max_d, peak_threshold)/2;
%plot
findpeaks(dd,f,...    
    'MinPeakProminence',min_peak_proms,...
    'MinPeakWidth', 5, ...
    'MinPeakDistance',10,...
    'Annotate','extents');
title(['De-noised signal,pillar ID=',num2str(pillar_id)])
xlabel('frame ID')
ylabel('deflection (nm)')
%find again to return the peaks
[pks,locs,widths,proms] = findpeaks(dd,f,...    
    'MinPeakProminence',min_peak_proms,...
    'SortStr', 'descend',...    
    'MinPeakWidth', 5, ...
    'MinPeakDistance',10,...
    'Annotate','extents');

%% display the maximum peak capped with threshold if existed.
found = pks<peak_threshold;
pks1=pks(found);
locs1=locs(found);
widths1=widths(found);
proms1=proms(found);
if(numel(pks1)>0)
    lc=locs1(1);
    pk=pks1(1);
    fprintf('The maximum peak:%f (<%f)\n location(frame ID):%d\n prominence:%f\n widths(half-height):%f\n', pk, peak_threshold, lc, proms1(1), widths1(1));
    hold on; plot(lc,pk,'o'); hold off
else
    disp('no peaks found!')
end    

%% display the maximum peak without thresholding
% if(numel(pks)>0)
%     lc=locs(1);
%     pk=pks(1);
%     fprintf('The maximum peak:%f\n location(frame ID):%d\n prominence:%f\n widths(half-height):%f\n', pk, lc, proms(1), widths(1));
%     hold on; plot(lc,pk,'x'); hold off
% else
%     disp('no peaks found!')
% end
