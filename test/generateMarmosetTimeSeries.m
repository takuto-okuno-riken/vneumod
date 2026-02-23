% generate BOLD signal time-series based on marmoset ROI atlas
% by PD marmoset of resting-state fMRI

function generateMarmosetTimeSeries
    % parameters
    smooth = 34;
    filter = ''; %'hf';%
    nuisance = 'gmacomp'; % 'aro'; %
    dtype = '';

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
    rmFrame = 5;  % number of removing frames
    TR = 1.0;

    % gaussian filter
    FWHM = [smooth/10 smooth/10 smooth/5]; % voxel size;
    sigma = FWHM / sqrt(8*log(2));
    filterSize = 2*ceil(2*sigma)+1;
    hpfTh = 1 / 128; % highpass filter threthold (Hz)

    rsfmribase = {['G:\marmoset\pd_org']};
    maxfiles = inf;

    % load cube atlas 
    atlas = ['data/BMA2marmo' num2str(atlasSize) 'atlas.nii' ];
    atlasinfo = niftiinfo([atlas '.gz']);
    atlasV = niftiread(atlasinfo);
    atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect
    
    % for Nuisance Signal Regression
    csfF = 'data\csf_BM2RsfMRI.nii.gz';
    csfinfo = niftiinfo(csfF);
    csfV = niftiread(csfinfo);
    csfV = adjustVolumeDir(csfV, csfinfo);
    csfV = single(csfV) / 255; % to [0 1] range
    wmF = 'data\white_BM2RsfMRI.nii.gz';
    wminfo = niftiinfo(wmF);
    wmV = niftiread(wminfo);
    wmV = adjustVolumeDir(wmV, wminfo);
    wmV = single(wmV) / 255; % to [0 1] range
    gsF = ['data/BMA2_T2wi_100muRsfMRI.nii.gz'];
    gsinfo = niftiinfo(gsF);
    gsV = niftiread(gsinfo);
    gsV = adjustVolumeDir(gsV, gsinfo);
    gsV(gsV>=1) = 1;
    gsV(gsV<1) = 0;

    % calc beta values (regression)
    cubename = ['BMA2marmo' num2str(atlasSize)];
    sessionName = ['x' cubename];
    outcnt = 0;

    outpath = ['results/ts' num2str(atlasSize)];
    if ~exist(outpath,'dir')
        mkdir(outpath);
    end
    for j=1:length(rsfmribase)
        listing = dir([rsfmribase{j} '\w*.nii.gz']);
        for i=1:length(listing)
            if listing(i).isdir == 1, continue; end
            if outcnt >= maxfiles, continue; end
            rsfmri = listing(i).name;
            names = split(rsfmri,'.');
            ids = split(names{1},'_');
            id = ids{1};

            % check output file first
            xmatf = [outpath '/' sessionName 's' num2str(smooth) filter nuisance '_' id '.mat'];
            if exist(xmatf,'file')
                disp(['file found (skipped) : ' xmatf]);
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
    
            % cut first frames
            disp([num2str(size(V,4)) ' frames. cut first ' num2str(rmFrame) ' frames.']);
            V = V(:,:,:,rmFrame+1:end);

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
            X = X';

            xstd = std(X(:),1);
            if xstd > 10
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
    % for unity
    fid = fopen(['results/dbs/' sessionName '.dat'],'w');
    for i=1:size(Zall,1)
        fwrite(fid,single(Zall(i,1:1190*20)),'single');
    end
    fclose(fid);
    return;
%}
end
