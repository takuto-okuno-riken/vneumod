% calc DBS simulated VAR surrogate of PD rs-fMRI cube ROI time-series
% testDBSrestC.m, testDBSrestSurr.m run first.

function testDBSrestSimAddMul
    algos = {'var','pc80','pcvar','pc999','pc9999','rdg01','rdg05','rdg800','rdg8000','las02','las05','las08','mkvar','mk50','mk80'}; % 'sigvar','sigmvar', % sigvar not work. sigmvar, soso.
    smooth = 's34';
    nuisance = 'gmacomp'; % 'aro'; %
    atlasSizes = [3, 2];
    usegpus = [false, false];
    lags = [1,2,3,4,5];
%    dbsroi = [45];%, 12, 777]; % STN, GPi, non
%    dbsroi = [4501, 4502, 4503]; % STN anteri, center lateral, center medial (sz=2) (4502==4525)
%    dbsroi = [4504]; % STN limbic asso (sz=2)
%    dbsroi = [4511:4515]; % STN anteri, center, post (sz=2)
%    dbsroi = [4521:4527]; % STN full voxel (sz=2)
    dbsroi = 4525; % STN center lateral
%    dbsroi = [20011:20019]; % PPN full voxel (sz=2)
%    dbsroi = 20017; % PPN ventral posterior 
%    dbsroi = [1221:1240]; % GPi full voxel (sz=2)
%    dbsroi = 1230; %[1230, 1236]; % GPi posterior medial
%    dbsroi = [3511:3546]; % Vim/Vop full voxel (sz=2)
    dbsroi = 3518; % Vim-Vo posterior 3527; % Vim
%    dbsroi = [21001:21008]; % PSA (Vim bottom) full voxel (sz=2)
%    dbsroi = [1201, 1202, 1203]; % GPi post, anteri (sz=3)
%    dbsroi = [1211, 3502, 20001]; % GPi post, VLp (VIM), PPN rostral (sz=2)
%    dbsroi = [3503, 3504, 3402, 3403]; % VIM, VIM, Vop, Vop (sz=2)
%    dbsroi = [1212, 3403, 20002]; % GPi asso, Vop, PPN asso (sz=2)
    side = [0, 1]; % both, left
    dtype = ''; % for pd1s 'hcp1s'; % 'pd'; % 'hc'; % 
    sbjmax = ''; % for others '61'; %'30'; % for pd30 
    mtype = ''; % for whole brain 'Ecp'; % except cerebellum & pons 
    kfold = 1;%  10; %if 1 no fold.

    % 12 -- left, right GPi;
    % 114 -- left, right SN (substantia nigra)
    % 45 -- left, right STH (subthalamic nucleus)

    for a=1 %1:3
        algo = algos{a};
        for sz=2 %1:1 % 2:2 %
            atlasSize = atlasSizes(sz);
            for p=1:1
                net = []; CXall = {};
                for pr=[15,38]%,44,47,59,60,69,74,81,97]
%                    permstr = num2str(590*pr);
%                    permstr = sprintf('pm%02d',pr);
%                    permstr = sprintf('ng%02d',pr);
%                    permstr = sprintf('sb%02d',pr);
                    permstr = sprintf('sp%02d',pr);
                    for tg=1:length(dbsroi)
                        for d=1:1
                            for kf=1:kfold
                                [net, CXall] = checkDbsVarSurrogateAlgos(algo, atlasSize, lags(p), usegpus(sz), dbsroi(tg), side(d), smooth, nuisance, dtype, sbjmax, mtype, permstr, kfold, kf, net, CXall);
                            end
                        end
                    end
                end
            end
        end
    end
end

function [net, CXall] = checkDbsVarSurrogateAlgos(algo, atlasSize, lag, usegpu, dbsroi, side, smooth, nuisance, dtype, sbjmax, mtype, permstr, kfold, kf, net, CXall)
    res = 16; % HRF sampling resolution
    sp = 8;   % HRF sampling starting point
    TR = 1.0;

    tuM = 8;  % GLM tukey-taper size
    Pth = 1e-4; % GLM height threshold
    cth = 30;   % GLM cluster threshold
    isRtoL = false;
    contnames = {'STH DBS'}; % GLM contrust name
    contrasts = {[1 0]'}; % GLM contrust

    surrNum = 40; % surrogate number

    % atlas of cube clusters
    % need to run testDBS5.m first
    atlas = ['data/allen' mtype 'Cube' num2str(atlasSize) 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);
    atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect
    % back ground image template
    backf = 'data/MNI152_T1_2mm_brain.nii.gz';
    backinfo = niftiinfo(backf);
    backV = niftiread(backinfo);

    cubename = ['Allen' mtype 'Cube' num2str(atlasSize)];
    % calc var net values (regression)
    if lag==1
        lagStr = '';
    else
        lagStr = ['L' num2str(lag)];
    end

    % prepare small brain mask volume for estimateSmoothFWHM
    if atlasSize > 1
        estiV = resamplingNifti3DVolume(atlasV, atlasSize, atlasSize, 'max');
    else
        estiV = atlasV;
    end

    % load ROI index 
    if any(ismember(dbsroi, [21001:21008]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'PSA.mat']; % sub parts version
    elseif any(ismember(dbsroi, [3511:3546]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'VL.mat']; % sub parts version
    elseif any(ismember(dbsroi, [1221:1240]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'GPi.mat']; % sub parts version
    elseif any(ismember(dbsroi, [20011:20019]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'Ppn.mat']; % sub parts version
    elseif any(ismember(dbsroi,[4521:4527]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'Stn3.mat']; % sub parts version
    elseif dbsroi == 4504
        ROIidxf = ['data/allen' mtype 'StnlimassoROIindex' num2str(atlasSize) '.mat']; % sub parts version
    elseif dbsroi == 3503 || dbsroi == 3504 || dbsroi == 3402 || dbsroi == 3403
        ROIidxf = ['data/allen' mtype 'VopvimROIindex' num2str(atlasSize) '.mat']; % sub parts version
    elseif dbsroi == 1212 || dbsroi == 3401 || dbsroi == 20002
        ROIidxf = ['data/allen' mtype 'PpnassoROIindex' num2str(atlasSize) '.mat']; % sub parts version
    elseif dbsroi == 20001
        ROIidxf = ['data/allen' mtype 'PpnmotROIindex' num2str(atlasSize) '.mat']; % sub parts version
    elseif dbsroi == 1211 || dbsroi == 3502
        ROIidxf = ['data/allen' mtype 'StnvimppnROIindex' num2str(atlasSize) '.mat']; % sub parts version
    elseif any(ismember(dbsroi,[1201:1203]))
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) 'GPi.mat']; % sub parts version
    else
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) '.mat']; % original
    end
    load(ROIidxf);
    if side > 0
        dbsidx = find(ROIidx==dbsroi & SDidx==side); % find left target index
        sidestr = ['-' num2str(side)];
    else
        dbsidx = find(ROIidx==dbsroi); % find target index
        sidestr = '';
    end

    % block number of surrogate frame length
    if strcmp(dtype,'pd') || strcmp(dtype,'hc')
        TR = 2.5;
    else
        TR = 1;
    end

    % load cube ROI time-series (CX from testDBSrestSurr.m)
    if kfold==1
        kfolddir = '';
        kfoldstr = '';
    else
        kfolddir = [num2str(kfold) 'fold'];
        kfoldstr = [num2str(kfold) 'fold'  '-' num2str(kf)];
    end
    if isempty(net)
        fname = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurr' algo lagStr cubename smooth nuisance dtype sbjmax kfoldstr '.mat'];
        load(fname);
    end
    if isempty(CXall)
        cxfname = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurrCX' cubename smooth nuisance dtype sbjmax '.mat'];
        cf=load(cxfname);
        CXall = cf.CX;
    end
    step = length(CXall) / kfold;
    CX = CXall;
    if kfold > 1, CX(round(1+step*(kf-1)):round(step*kf)) = []; end

    % find time series range
    cxlen = length(CX);
    X3 = nan(size(CX{1},1),size(CX{1},2),cxlen,'single');
    for nn=1:cxlen
        X = CX{nn};
        X3(1:size(X,1),1:size(X,2),nn) = X;
    end
    frames = size(CX{1},2);
    m = nanmean(X3(:));
    s = nanstd(X3(:),1);
    clear X3;

    % load surrogate permutation (seed)
    perm = []; C = []; Err = [];
    permf = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurrPerm' permstr algo lagStr cubename smooth nuisance dtype sbjmax kfoldstr '.mat'];
    if exist(permf,'file')
        load(permf);
    else
        uxtime = uint32(posixtime(datetime('now')));
        if length(permstr)>2 && strcmp(permstr(1:2),'pe')
            % noise permutation in each subject
            rng(uxtime);
            perm = randperm(frames);
        elseif length(permstr)>2 && strcmp(permstr(1:2),'sb')
            % ordered residual with same subject
            perm = (1:frames) + (str2double(permstr(3:4))-1)*frames;
        elseif length(permstr)>2 && strcmp(permstr(1:2),'sp')
            % ordered residual with subject permutation
            rng(uxtime);
            rp = randperm(cxlen);
            for i=1:cxlen
                perm = [perm, (1:frames) + (rp(i)-1)*frames];
            end
        end
        if ~isempty(perm)
            save(permf,'uxtime','perm','-v7.3');
        end
    end
    if length(permstr)>2 && strcmp(permstr(1:2),'ng')
       dist = 'gaussian';
    else
       dist = 'residuals';
    end

    dbspws = 0.15; %0.2; %0.6; %[0.25 0.5 0.6 0.75 1]; % for atlasSize=2
%    dbspws = [0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5]; % for atlasSize=2
    srframes = 160; dbsoffsec = 28; dbsonsec = 22;
    for p = 1:length(dbspws)
        dbspw = dbspws(p);
        sessionName = ['testdbsSim' num2str(dbsroi) sidestr 'AddMul' num2str(srframes) '-' num2str(dbsoffsec) '-' num2str(dbsonsec) '-' num2str(dbspw) 'pw' num2str(surrNum) 'sr'  permstr algo lagStr cubename smooth nuisance dtype sbjmax kfoldstr];

        if exist(['results/dbs' dtype num2str(atlasSize) 'nii' kfolddir '/' sessionName '2nd-mix-Tukey' num2str(tuM) '.nii.gz'],'file'), continue; end

        % DBS block design 30s on, 60s off (off, on, off, on, ...)
        CA = cell(1,surrNum); Chrf = cell(1,surrNum); CM = cell(1,surrNum); 
        dt = TR / res;
        [t, hrf] = getGlmHRF(dt); % human's HRF;

        for i=1:surrNum
            n = size(CX{i},1);
            bmax = floor((srframes * TR) / (dbsoffsec+dbsonsec));
            ons = []; dur = [];
            if strcmp(dtype,'pd') || strcmp(dtype,'hc')
                for j=0:bmax-1
                    ons(j+1) = dbsoffsec + j*(dbsoffsec+dbsonsec);
                    dur(j+1) = dbsonsec;
                end
            else
                for j=0:bmax-1
                    ons(j+1) = dbsoffsec + j*(dbsoffsec+dbsonsec);
                    dur(j+1) = dbsonsec;
                end
            end
            onsets{1} = ons;
            durations{1} = dur;
            [Chrf{i}, U] = getGlmHRFDesignMatrix(onsets, durations, srframes, TR, res, sp, hrf);
            block = nan(n,srframes,'single');
            block(dbsidx,:) = repmat(Chrf{i}' * dbspw * s,length(dbsidx),1);
            CA{i} = block;
            block(dbsidx,:) = repmat(1 - 0.5 * Chrf{i}',length(dbsidx),1);
            CM{i} = block;
        end
    
        % simulated VAR surrogate
        S = [];
        fname = ['results/dbs' dtype num2str(atlasSize) '/' sessionName 'Surr.mat'];
        if exist(fname,'file')
            load(fname);
        end
        if length(S) < surrNum
            if isempty(S), S = cell(surrNum,1); end

            % surrogate time-series
            for i=1:surrNum
                X = [];
                for k=0:cxlen-1
                    X=[X, CX{mod((i-1)*2+k,cxlen)+1}];
                    if size(X,2) >= srframes, break; end
                end
                if i > length(S) || isempty(S{i})
                    tc = tic;
                    if ~isnan(str2double(permstr)) % ordered residual from nBaset
                        nBaset = frames*(i-1) + str2double(permstr);
                    elseif length(permstr)>2 && strcmp(permstr(1:2),'pe') % residual permutation in each subject
                        nBaset = {perm, frames*(i-1)};
                    elseif length(permstr)>2 && strcmp(permstr(1:2),'sb') % ordered residual with same subject
                        nBaset = {perm, 0};
                    elseif length(permstr)>2 && strcmp(permstr(1:2),'sp') % ordered residual with subject permutation
                        nBaset = {perm(frames*(i-1)+1:frames*i), 0};
                    else
                        nBaset = perm; % first time, perm is empty (full residual permutation or gaussian error)
                    end
                    if ~isempty(Err), dist ='residuals'; end % use Err as is

                    [S{i}, C, Err, gperm] = surrogateDbsMVAR(X(:,1:srframes), [], [], [], net, CA{i}, CM{i}, dist, 1, NaN, nBaset, C, Err, usegpu);
                    disp(['done t=' num2str(toc(tc)) 'sec']);

                    if isempty(nBaset), perm = gperm; end
                end
            end

            % save data
%            save(fname,'S','-v7.3');
        end
        if ~exist(permf,'file') && ~isempty(perm)
            if length(permstr)<=2 || ~strcmp(permstr(1:2),'ng'), Err = []; end %  Err is big. do not save if not necessary.
            save(permf,'perm','Err','-v7.3');
        end

        % calc 1st-level GLM
%%{
        bmatC = cell(1,surrNum);
        for i=1:surrNum
            f = struct();
            % calc 1st level GLM
            Xorg = Chrf{i};
            Xt = [Xorg, ones(size(Xorg,1),1)];
            [f.B2, RSS, df, X2is, tRs, R] = calcGlmTukey(S{i}', Xt, tuM);

            [recel, f.FWHM] = estimateSmoothFWHM(R, RSS, df, estiV);
            bmatC{i} = f;
        end
%}
        % calc 2nd-level estimation
        B1 = [];
        X2 = [];
        FWHMs = [];
        for i=1:surrNum
            f = bmatC{i};
            % 2nd-level Y vector
            B2 = f.B2(:,[1,2]); % include design and intercept (we need more than 8 length for tukey taper)
            B1 = [B1; B2'];
            FWHMs = [FWHMs; f.FWHM];
    
            % 2nd-level design matrix
            X2 = [X2; eye(size(B2,2))];
        end
        clear f;
        B1(isnan(B1)) = 0; % there might be nan
        FWHMs = mean(FWHMs,1); % let's take the mean of FWHM.
    
%        betaBmat = ['results/dbs' dtype num2str(atlasSize) '/' sessionName '2nd-mix-Tukey' num2str(tuM) '.mat'];
%        if exist(betaBmat,'file')
            % load beta volumes
%            load(betaBmat);
%        else
            % calc 2nd-level estimation
            [B, RSS, df, X2is, tRs, R] = calcGlmTukey(B1, X2, tuM);
    
            [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, estiV);
    
            % output beta matrix
%            save(betaBmat,'B','RSS','X2is','tRs','recel','FWHM','df','-v7.3');
%        end

        % GLM contrast images
        Ts = calcGlmContrastImage(contrasts, B, RSS, X2is, tRs);

        % GLM contrast image
        thParam = {df, Pth};
        clParam = {cth, FWHM}; % clustering parameter for GLM contrast
        [Tth, Vts, Vfs, Tmaxs, Tcnts] = plotGlmContrastImage(contnames, Ts, thParam, clParam, atlasV, (atlasSize==1), isRtoL, backV, ...
            [sessionName '2nd-mix-Tukey' num2str(tuM)]); %, rangePlus, rangeMinus, [], [], []);

        % save T-value NIfTI volume
%        saveNifti(backf,Vts,['results/dbs' dtype num2str(atlasSize) 'nii' kfolddir '/'],[sessionName '2nd-mix-Tukey' num2str(tuM) 'th']);
        saveNifti(backf,Vfs,['results/dbs' dtype num2str(atlasSize) 'nii' kfolddir '/'],[sessionName '2nd-mix-Tukey' num2str(tuM)]);
    end
end

%% should save one file
function saveNifti(tfmri, V2s, path, outname)
    info = niftiinfo(tfmri);
    info.ImageSize = info.ImageSize(1:3);
    info.PixelDimensions = info.PixelDimensions(1:3);
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    for j=1:length(V2s)
        fname = [path outname '.nii'];
        V = adjustVolumeDir(V2s{j}, info); 
        niftiwrite(V,fname,info,'Compressed',true);
    end
end

