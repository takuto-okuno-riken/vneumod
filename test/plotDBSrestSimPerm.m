
function plotDBSrestSimPerm
    % parameters
    algo = 'var';
    atlasSize = 2; %3; %
    smooth = 's34';
    nuisance = 'gmacomp'; %'aro'; %
    atlas = 'Allen';

    path = 'results/';
    dlabels = {'BA4,6','cortex','subcortex','cerebellum','ctx,cereb','All'};

    % check surrogate permutation seed and AUC (ROI 4525). fixed power 0.15 (work)
    checkDBSpermseedPw015_06(algo, atlasSize, smooth, nuisance, atlas, path, dlabels);
end

function [gV, aV, bmV, bcV] = loadAtlasfiles()
    info=niftiinfo('D:\work\DBS\DataForFigures\Figure3\Figure 3A\Optimal\spmT_000mirror.nii.gz');
%    info=niftiinfo('D:\work\DBS\DataForFigures\Figure3\Figure 3B\Optimal_5.5V\spmT_000mirror.nii.gz');  % supress DMN
%    info=niftiinfo('D:\work\DBS\DataForFigures\FigureS2\S2A_Optimal_contact+voltage\spmT_000mirror.nii.gz');
    gV = niftiread(info);
    gV = adjustVolumeDir(gV, info);
    
    atlasinfo = niftiinfo('D:\work\defaultmode\human\itksnap_annotation_full2mm.nii.gz');
    aV = niftiread(atlasinfo);
    aV = adjustVolumeDir(aV, atlasinfo);
    brodminfo = niftiinfo('D:\work\defaultmode\human\brodmann_MMrep2mm.nii.gz');
    bmV = niftiread(brodminfo);
    bmV = adjustVolumeDir(bmV, brodminfo);
    brodcinfo = niftiinfo('D:\work\defaultmode\human\brodmannConnfMRI.nii.gz');
    bcV = niftiread(brodcinfo);
    bcV = adjustVolumeDir(bcV, brodcinfo);
%{
    % find target voxel in cerebellum clus1 & 2
    V = gV; V(:) = 0;
    V(bcV>=107) = 1;
    V(bcV>=107&bcV<=110) = 2;
    niftiwrite(V,'D:\work\defaultmode\human\brodmannCrus12_2mm.nii',info,'Compressed',true);
%}
end

function checkDBSpermseedPw015_06(algo, atlasSize, smooth, nuisance, atlas, path, dlabels)
    dbsrois =  [4525];%[4501, 4502, 4503]; % STH anteri, post, other (sz=2)
    tuM = 8;  % GLM tukey-taper size
    cubename = [atlas 'Cube' num2str(atlasSize)];

    % atlas of cube clusters
    % need to run testDBS5.m first
    atlas = ['data/' cubename 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);
    atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect

    [gV, aV, bmV, bcV] = loadAtlasfiles();

    n = single(max(atlasV(:)));
%    BPth = 0.05 / n; % Bonferroni correction
    Tth = abs(tinv(1e-4,n-1));
    gTruth = abs(gV)>Tth;

    surrNum = 40;
    permNum = 100;

    aucth = 0.70;
%    tidx = find(dbsroi==4502);
    tidx = find(dbsrois==4525);
    targets = [];
    rs = nan(length(dbsrois),length(permNum),6);
    aucs = nan(length(dbsrois),length(permNum),6);
    for i = 1:length(dbsrois)
        dbsroi = dbsrois(i);
        permstrs = {};
        for p = 1:permNum
            sessionName = ['ppmi81CX' cubename smooth nuisance '_' num2str(dbsroi) 'sr' num2str(surrNum) 'pr' num2str(p)];
            fname = [path sessionName '_2nd-Tukey' num2str(tuM) '.nii.gz'];
            if exist(fname,'file')
                info = niftiinfo(fname);
                V = niftiread(info);
                V = adjustVolumeDir(V, info);
            else
                continue;
            end
            [rs(i,p,:), aucs(i,p,:)] = calcPartCorrAndAUC(V, gV, bmV, aV, bcV, gTruth);
            permstrs{p} = num2str(p);
        end

        % check DBS power effect
        figure; plot(squeeze(rs(i,:,:))); legend(dlabels); ylabel('R'); xlabel('permNum'); xticks(1:permNum); xticklabels(permstrs); % plot correlation
        title(['correlation roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]);

        figure; plot(squeeze(aucs(i,:,:))); legend(dlabels); ylabel('AUC'); xlabel('permNum'); xticks(1:permNum); xticklabels(permstrs); % plot AUC
        title(['detection roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]);
    end

    figure; boxplot(aucs(:,:,6)'); ylabel('AUC');
    title(['(all) roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]); 
    figure; boxplot(rs(:,:,6)'); ylabel('R')
    title(['(all) roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]);

    % sort and find top 10
    scores = abs(aucs-0.5)*2 + abs(rs);
    sc6 = squeeze(scores(:,:,6));
    [sc6des,idx] = sort(sc6,'descend');
    figure; plot(sc6des); ylabel('score');
    title(['(all) roi=' num2str(dbsroi) ' Tth=' num2str(Tth)]); 
end

function [rs, aucs] = calcPartCorrAndAUC(V, gV, bmV, aV, bcV, gTruth)
    rs = nan(6,1);
    aucs = nan(6,1);
    idx = {};
    idx{1} = find((bmV==4 | bmV==6) & aV>0 & gTruth); % brodmann 4,6
    idx{2} = find(bmV>0 & aV>0 & gTruth); % cortex
    idx{3} = find(bmV==0 & ~(bcV>=107) & aV>0 & gTruth); % subcortex
    idx{4} = find(bcV>=107 & aV>0 & gTruth); % cerebellum
    idx{5} = find((bcV>=107 | bmV>0) & aV>0 & gTruth); % cortex & cerebellum
    idx{6} = find(aV>0 & gTruth); % cortex & subcortex & cerebellum
    for j = 1:6
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

