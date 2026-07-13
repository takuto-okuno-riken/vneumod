
function plotMarmosetDBSrestSimPerm
    % parameters
    algo = 'var';
    atlasSize = 2; %3; %
    smooth = 's34';
    nuisance = 'gmacomp'; %'aro'; %
    atlas = 'BMA2';

    path = 'results/';
    dlabels = {'BA4,6','cortex','subcortex','cerebellum','ctx,cereb','All'};

    % check surrogate permutation seed and AUC (ROI 4525). fixed power 0.15 (work)
    checkDBSpermseedPw015_06(algo, atlasSize, smooth, nuisance, atlas, path, dlabels);
end

function [gV, aV] = loadAtlasfiles()
    info=niftiinfo('data\warped_spmT_000mirror_WarpedRsfMRI.nii.gz');
    gV = niftiread(info);
    gV = adjustVolumeDir(gV, info);

    atlasinfo = niftiinfo('data/BMA2_labels_100muRsfMRI.nii.gz');
    aV = niftiread(atlasinfo);
    aV = adjustVolumeDir(aV, atlasinfo);
end

function checkDBSpermseedPw015_06(algo, atlasSize, smooth, nuisance, atlas, path, dlabels)
    dbsrois =  [2000 2001 2002 2003]; % STH voxels (sz=2)
    tuM = 8;  % GLM tukey-taper size
    cubename = [atlas 'marmo' num2str(atlasSize)];

    % atlas of cube clusters
    atlas = ['data/' cubename 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);
    atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect

    [gV, aV] = loadAtlasfiles();

    n = single(max(atlasV(:)));
    BPth = 0.05;
%    BPth = 1e-4; % 0.0001
%    BPth = 0.05 / n; % Bonferroni correction
    Tth = abs(tinv(BPth,n-1));
    gTruth = abs(gV)>Tth;

    surrNum = 30;
    permNum = 100;

    rs = nan(length(dbsrois),permNum,6);
    aucs = nan(length(dbsrois),permNum,6);
    for i = 1:length(dbsrois)
        dbsroi = dbsrois(i);
        permstrs = {};
        for p = 1:permNum
            sessionName = ['pd36' cubename smooth nuisance '_' num2str(dbsroi) 'sr' num2str(surrNum) 'pr' num2str(p)];
            fname = [path sessionName '_2nd-Tukey' num2str(tuM) '.nii.gz'];
            if exist(fname,'file')
                info = niftiinfo(fname);
                V = niftiread(info);
                V = adjustVolumeDir(V, info);
            else
                continue;
            end
            [rs(i,p,:), aucs(i,p,:)] = calcPartCorrAndAUC(V, gV, aV, gTruth);
            permstrs{p} = num2str(p);
        end

        % check DBS power effect
        figure; plot(squeeze(rs(i,:,:))); legend(dlabels); ylabel('R'); xlabel('permNum'); xticks(1:permNum); xticklabels(permstrs); % plot correlation
        title(['correlation roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]);

        figure; plot(squeeze(aucs(i,:,:))); legend(dlabels); ylabel('AUC'); xlabel('permNum'); xticks(1:permNum); xticklabels(permstrs); % plot AUC
        title(['detection roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]);
    end

    figure; boxplot(aucs(:,:,6)'); ylabel('AUC'); xlabel('DBS ROI'); xticklabels(dbsrois);
    title(['(all) Tth=' num2str(Tth)]); 
    figure; boxplot(rs(:,:,6)'); ylabel('R'); xlabel('DBS ROI'); xticklabels(dbsrois);
    title(['(all) Tth=' num2str(Tth)]);

    scores = (aucs-0.5)*2 + rs;
    sc6 = squeeze(scores(:,:,6));
    figure; boxplot(sc6'); hold on; plot(sc6,'Color',[0.7 0.7 0.7]); ylabel('score');
    title(['(all) Tth=' num2str(Tth)]); xlabel('DBS ROI'); xticks(1:length(dbsrois)); xticklabels(dbsrois);

    % sort and find top 10 (all, score)
    [sc6des,sidx] = sort(sc6,2,'descend');
    figure; hold on;
    for i=1:10
        idx = sidx(2,i);
        plot(aucs(:,idx,6));
    end
    title(['All, Top10 scores Tth=' num2str(Tth)]); xlabel('DBS ROI'); xticks(1:length(dbsrois)); xticklabels(dbsrois);

    % sort and find top 10 (BA4,6, auc)
    au1 = squeeze(aucs(:,:,1));
    [sc6des,aidx] = sort(au1,2,'descend');
    figure; hold on;
    for i=1:10
        idx = aidx(2,i);
        plot(aucs(:,idx,1));
    end
    title(['BA4,6 Top10 AUCs Tth=' num2str(Tth)]); xlabel('DBS ROI'); xticks(1:length(dbsrois)); xticklabels(dbsrois);
    figure; hold on;
    for i=1:10
        idx = aidx(2,i);
        plot(rs(:,idx,1));
    end
    title(['BA4,6 Top10 Rs Tth=' num2str(Tth)]); xlabel('DBS ROI'); xticks(1:length(dbsrois)); xticklabels(dbsrois);
end

function [rs, aucs] = calcPartCorrAndAUC(V, gV, aV, gTruth)
    rs = nan(6,1);
    aucs = nan(6,1);
    idx = {};
    idx{1} = find(((aV>=31 & aV<=37) | (V>=10031 & V<=10037)) & gTruth); % brodmann 4,6
    idx{2} = find(((aV>=26 & aV<=150) | (aV>=10026 & aV<=10150)) & gTruth); % cortex
    idx{3} = find(((aV>=151 & aV<=632) | (aV>=10151 & aV<=10632)) & gTruth); % subcortex
    idx{4} = find(((aV>=787 & aV<=831) | (aV>=10787 & aV<=10831)) & gTruth); % cerebellum
    idx{5} = unique([idx{2}(:); idx{4}(:)]); % cortex & cerebellum
    idx{6} = find(aV>0 & gTruth); % cortex & subcortex & cerebellum

    for j = 1:length(idx)
        X = V(idx{j});
        G = gV(idx{j});
        nanidx = find(isnan(X) | isnan(G));
        X(nanidx) = []; G(nanidx) = [];

        rs(j) = corr(X,G,'type','pearson','rows','pairwise');
        [~, ~, aucs(j)] = calcROCcurvePlMi(X, G);
    end
end

function [X, Y, AUC] = calcROCcurvePlMi(T, G)
    n = length(T(:));
    X = nan(n,1);
    Y = nan(n,1);
    [~,I] = sort(abs(T(:)),1,'descend');
    TG = T(:) .* G(:);
    tc = sum(TG>0);
    fc = n - tc;
    q = T(I);
    r = G(I);
    x = 0; y = 0;
    for j=1:n
        if (q(j) > 0 && r(j) > 0) || (q(j) < 0 && r(j) < 0)
            y = y + 1/tc;
        else
            x = x + 1/fc;
        end
        X(j) = x;
        Y(j) = y;
    end
    if tc==0 || fc == 0
        X(end+1)=1; Y(end+1)=1;
    end
    AUC = trapz(X(:),Y(:));
end

function [X, Y, auc] = plotROCcurvePlMi(T, G)
    [X, Y, auc] = calcROCcurvePlMi(T, G);
    hold on;
    plot(X, Y);
    plot([0 1], [0 1],':','Color',[0.5 0.5 0.5]);
    hold off;    
    ylim([0 1]);
    xlim([0 1]);
    daspect([1 1 1]);
    title('ROC curve');
    xlabel('False Positive Rate')
    ylabel('True Positive Rate')
end

