% generate BOLD signal time-series based on cube atlas
% by PD subjects (male/female) of resting-state fMRI
% (testDBSrestC.m)

function generateSubjectTimeSeries
    % parameters
    smooth = 34;
    filter = ''; %'hf';%
    nuisance = 'gmacomp'; % 'aro'; %
    dtype = ''; % for pd1s % 'hcp1s'; % 'pd'; % 'hc'; %  'prod'; % 

    % atlas of cube clusters
    % need to run makeCubeAtlas.m first
    atlasSizes = [10, 9, 8, 7, 6, 5, 4, 3, 2, 1];
    for sz=9
        atlasSize = atlasSizes(sz);
        generateSubjectTimeSeriesByCubeAtlas(atlasSize,smooth,filter,nuisance,dtype);
    end
end

function generateSubjectTimeSeriesByCubeAtlas(atlasSize,smooth,filter,nuisance,dtype)
    path = ['results/dbs' dtype num2str(atlasSize) '/'];
    cmap = turbo;
%    cmap = hot;
    rmFrame = 10;  % number of removing frames
    TR = 1.0;

    % gaussian filter
    FWHM = [smooth/10 smooth/10 smooth/10]; % voxel size;
    sigma = FWHM / sqrt(8*log(2));
    filterSize = 2*ceil(2*sigma)+1;
    hpfTh = 1 / 128; % highpass filter threthold (Hz)

    if strcmp(dtype,'prod')
        rsfmribase = {['H:\PPMI\prod1s']};
        maxfiles = inf;
    elseif strcmp(dtype,'pd') % TR=2.5
        rsfmribase = {['H:\PPMI\pd']};
        maxfiles = 200;
    elseif strcmp(dtype,'hc') % TR=2.5
        rsfmribase = {['H:\PPMI\hc']};
        maxfiles = inf;
    elseif strcmp(dtype,'hcp1s')
        rsfmribase = {['H:\PPMI\hcp1s']};
        maxfiles = inf;
    else
        rsfmribase = {['H:\PPMI\pd1s']};
        maxfiles = inf;
    end

    % load cube atlas 
    atlas = ['data/allenCube' num2str(atlasSize) 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);
    atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect
    % get downsampled ROI & Network index (testSeed5.m)
    ROIidxf = ['data/allenROIindex' num2str(atlasSize) '.mat'];
    load(ROIidxf);
    
    % for Nuisance Signal Regression
    csfF = 'data\csf.nii.gz';
    csfinfo = niftiinfo(csfF);
    csfV = niftiread(csfinfo);
    csfV = adjustVolumeDir(csfV, csfinfo);
    csfV = single(csfV) / 255; % to [0 1] range
    wmF = 'data\white.nii.gz';
    wminfo = niftiinfo(wmF);
    wmV = niftiread(wminfo);
    wmV = adjustVolumeDir(wmV, wminfo);
    wmV = single(wmV) / 255; % to [0 1] range
    gsF = ['data/itksnap_annotation_full2mmMask.nii.gz']; % MNI152_T1_2mm_brain_mask.nii.gz';
    gsinfo = niftiinfo(gsF);
    gsV = niftiread(gsinfo);
    gsV = adjustVolumeDir(gsV, gsinfo);
    gsV(gsV>=1) = 1;
    gsV(gsV<1) = 0;
    % ICA-AROMA mask files
%{
    arocsfF = 'D:\work\ICA-AROMA-master\mask_csf.nii.gz';
    info = niftiinfo(arocsfF);
    arocsfV = niftiread(info);
    arocsfV = adjustVolumeDir(arocsfV, info);
    aroegF = 'D:\work\ICA-AROMA-master\mask_edge.nii.gz';
    info = niftiinfo(aroegF);
    aroegV = single(niftiread(info));
    aroegV = adjustVolumeDir(aroegV, info);
    arooutF = 'D:\work\ICA-AROMA-master\mask_out.nii.gz';
    info = niftiinfo(arooutF);
    arooutV = niftiread(info);
    arooutV = adjustVolumeDir(arooutV, info);
    icabase = {'D:\work\melodic\pd1s34ASingles'};
%}
    % calc beta values (regression)
    cubename = ['AllenCube' num2str(atlasSize)];
    sessionName = ['testdbs' cubename];
    outcnt = 0;

    % load demographics for pd1s
    democount = [];

    outpath = ['results/dbs' num2str(atlasSize)];
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    % nifti output path
%{
    niftioutpath = ['results/dbs' num2str(atlasSize) 'nii']; %'';
    if ~exist(niftioutpath,'dir')
        mkdir(niftioutpath);
    end
%}
    for j=1:length(rsfmribase)
        if strcmp(dtype,'hcp1s')
            listing = dir([rsfmribase{j} '\wu*.nii.gz']);
        else
            listing = dir([rsfmribase{j} '\wau*.nii.gz']);
        end
        for i=1:length(listing)
            if listing(i).isdir == 1, continue; end
            if outcnt >= maxfiles, continue; end
            rsfmri = listing(i).name;
            names = split(rsfmri,'.');
            ids = split(names{1},'_');
            if strcmp(dtype,'hcp1s')
                id = extractAfter(ids{1},2);
            else
                id = extractAfter(ids{1},3);
            end
            date = ids{end-1};

            % check output file first
            xmatf = [path sessionName 's' num2str(smooth) filter nuisance '_' id '_' date '.mat'];
            if exist(xmatf,'file')
                disp(['file found (skipped) : ' xmatf]);
%{
                if ~isempty(niftioutpath)
                    load(xmatf);
                    V = getNifti4DFromRoiTS(X', atlasV);
                    V(isnan(V)) = 0;
                    info = niftiinfo(gsF);
                    info.ImageSize(4) = size(V,4);
                    info.PixelDimensions(4) = TR;
                    info.Datatype = 'single';
                    info.BitsPerPixel = 32;
                    V = adjustVolumeDir(V, info);
                    fname = [niftioutpath '/' sessionName 's' num2str(smooth) filter nuisance '_' id '_' date '.nii'];
                    niftiwrite(V,fname,info,'Compressed',true);
                end
%}
                continue;
            end
            
            % read fMRI volume
            fmri = [rsfmribase{j} '\' rsfmri];
            if ~exist(fmri,'file')
                disp(['file not found (skipped) : ' fmri]);
                continue;
            end
            disp(['loading : ' fmri]);
            mriinfo = niftiinfo(fmri);
            V = single(niftiread(mriinfo));
            V = adjustVolumeDir(V, mriinfo);
            V(isnan(V)) = 0;
    
            if (strcmp(dtype,'') || strcmp(dtype,'hcp1s') || strcmp(dtype,'prod')) && size(V,4) < 300
                disp(['too small frames = ' num2str(size(V,4))]);
                continue;
            elseif (strcmp(dtype,'pd') || strcmp(dtype,'hc')) && size(V,4) < 200
                disp(['too small frames = ' num2str(size(V,4))]);
                continue;
            end

            % read rp_*.txt (6 head motion parameters)
            if strcmp(nuisance, 'aro')
                rpf = [rsfmribase{j} '/rp_' names{1}(4:end) '.txt'];
                if ~exist(rpf,'file')
                    disp(['file not found (skipped) : ' rpf]);
                    continue;
                end
                disp(['loading : ' rpf]);
                M = readmatrix(rpf);
                Md = [zeros(1,6); diff(M,1,1)];
                icapath = [icabase{1} '/s34r' names{1}];
            end

            % cut first frames
            disp([num2str(size(V,4)) ' frames. cut first ' num2str(rmFrame) ' frames.']);
            V = V(:,:,:,rmFrame+1:end);
            if strcmp(nuisance, 'aro')
                M = M(rmFrame+1:end,:);
                Md = Md(rmFrame+1:end,:);
            end

            % gaussian filter
            disp(['sz=' num2str(smooth) ', sigma=' num2str(sigma(1)) ', flSz=' num2str(filterSize(1))]);
            for k=1:size(V,4)
                V(:,:,:,k) = imgaussfilt3(V(:,:,:,k), sigma, 'FilterSize', filterSize);
            end

            % nuisance regression
            if length(nuisance) > 0
                if strcmp(nuisance, 'gmacomp')
                    % get Nuisance time-series (Global Mean, CSF comps, WM comps)
                    Sd = getNuisanceMeanTimeSeries(V, [], [], []);
                    aComp = getNuisanceaCompCor(V, csfV, wmV, Sd);
                    Xn = [Sd, aComp];
                elseif strcmp(nuisance, 'aro')
                    % get Nuisance time-series by ICA-AROMA, coeffs are calculated by analyzeICAaromaHcp.m
                    Xn = getNuisanceICAaroma(TR, arocsfV, aroegV, arooutV, M, Md, icapath, 0.4, 0.02, [-0.37612, 0.98889]);
                end
%                figure; imagesc([Xn], [-4, 12]); colorbar;
                V = getNuisanceRegressionOut(V, Xn, gsV);
            end

            % ROI time-series from rs-fMRI
            X = getRoiTSFromNifti4D(V, atlasV, 'mean');
            X = X - nanmean(X,2);
%            X = convert2SigmoidSignal(X); % both work, convert or not convert
            X = X';

            xstd = std(X(:),1);
            if ~strcmp(dtype,'hcp1s') && xstd > 10
                disp(['std is too big =' num2str(xstd) ', ' id '_' date]);
                continue;
            end

            % high pass filter as preprocessing step (M.W.Woolrich, 2001) type.
            if strcmp(filter,'hf')
                disp(['apply highpass filter (' num2str(hpfTh) ' Hz) : tfMRI and design matrix...']);
                TR = mriinfo.PixelDimensions(4);
                X = highpass(X,hpfTh,1/TR);
            end

            % output X matrix
            save(xmatf,'X','-v7.3');
            outcnt = outcnt + 1;

            % sample plot of voxels
%{
            figure;
            for a=1:10
                hold on; plot(squeeze(X(:,15+a))); hold off;
            end
            title(['some voxels of ' id ' ' date]);
%}
            % show raster plot
%            figure; imagesc(X.'); title(['laster prot of ' id ' ' date]);
            % show histogram 
%            figure; histogram(X); title(['histogram of ' id ' ' date]);
        end
    end
%{
    female = length(find(democount(:,2)==0));
    male = length(find(democount(:,2)==1));
    mold = mean(democount(:,3)); sold = std(democount(:,3),1);
    disp(['male=' num2str(male) ', female=' num2str(female) ', ' num2str(mold) '±' num2str(sold) ' years old']);
%}
%{
    % for unity
    fid = fopen(['results/dbs/' sessionName '.dat'],'w');
    for i=1:size(Zall,1)
        fwrite(fid,single(Zall(i,1:1190*20)),'single');
    end
    fclose(fid);
    return;
%}
end
