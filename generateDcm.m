dcm = load("output/dcm.mat").dcm;

test = load("datafiles\ctd_iso_ALL.mat").ctd;

%%
maxPcm = nan(329,1); minPcm = nan(329,1);
meanPcm = nan(329,1); stdPcm = nan(329,1);

for i = 1:329
    %fcm(i,:) = test.fcm;
    tmp = test(i).pcm;
    if ~isempty(tmp)
        maxPcm(i) = max(tmp,[],'omitnan');
        minPcm(i) = min(tmp,[],'omitnan');
        meanPcm(i) = 2*round(mean(tmp/2,'omitnan'));
        stdPcm(i) = std(tmp,[],"omitnan");
    end
end

%% save
save output/dcm.mat meanPcm minPcm maxPcm stdPcm -append;