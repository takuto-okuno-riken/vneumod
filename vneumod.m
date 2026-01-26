%%
% virtual neuromodulation surrogate command line tool

function vneumod(varargin)

    % set version number
    versionNumber = '0.1';

    % add script path
    if ~isdeployed % checking MATLAB mode or stand-alone mode.
        [st,ind] = dbstack('-completenames');
        relpath = st(ind).file;
        [exedir,exename,ext] = fileparts(relpath);
        if exist([exedir '/util'],'dir')
            addpath([exedir '/util']);
            addpath([exedir '/lib']);
        end
    end

    % get exe file full path
    global exePath;
    global exeName;
    [exePath, exeName, ext] = exeFilename();

    % init command line input
    handles.commandError = 0;
    handles.matFiles = {};
    handles.outpath = 'results';
    handles.cx = [];
    handles.model = [];
    handles.pymodel = [];
    handles.roisidx = [];
    handles.atlas = [];
    handles.targatl = []; % target atlas
    handles.rois = [];
    handles.out = 1;
    handles.glm = 0;
    handles.nocache = 0;

    handles.surrNum = 40;
    handles.srframes = 160;
    handles.vnparam = [28,22,0.15];
    handles.tr = 1;
    handles.hrfparam = [16,8];

    % load command line input
    i = 1;
    while true
        if i > size(varargin, 2)
            break;
        end
        switch varargin{i}
            case {'-c','--cx'}
                handles.cx = varargin{i+1};
                i = i + 1;
            case {'-m','--model'}
                handles.model = varargin{i+1};
                i = i + 1;
            case {'--pymodel'}
                handles.pymodel = varargin{i+1};
                i = i + 1;
            case {'-i','--roiidx'}
                handles.roisidx = varargin{i+1};
                i = i + 1;
            case {'-a','--atlas'}
                handles.atlas = varargin{i+1};
                i = i + 1;
            case {'-t','--targatl'}
                handles.targatl = varargin{i+1};
                i = i + 1;
            case {'-r','--roi'}
                rois = varargin{i+1};
                if contains(rois,':')
                    s = split(rois,':');
                    handles.rois = str2num(s{1}):str2num(s{2});
                else
                    handles.rois = str2num(rois); % single or comma case 
                end
                i = i + 1;
            case {'-o','--out'}
                handles.out = str2num(varargin{i+1});
                i = i + 1;
            case {'--surrnum'}
                handles.surrNum = str2num(varargin{i+1});
                i = i + 1;
            case {'--srframes'}
                handles.srframes = str2num(varargin{i+1});
                i = i + 1;
            case {'--vnparam'}
                handles.vnparam(1) = str2num(varargin{i+1});
                handles.vnparam(2) = str2num(varargin{i+2});
                handles.vnparam(3) = str2num(varargin{i+3});
                i = i + 3;
            case {'--tr'}
                handles.tr = str2num(varargin{i+1});
                i = i + 1;
            case {'--hrfparam'}
                handles.hrfparam(1) = str2num(varargin{i+1});
                handles.hrfparam(2) = str2num(varargin{i+2});
                i = i + 2;
            case {'--outpath'}
                handles.outpath = varargin{i+1};
                i = i + 1;
            case {'--glm'}
                handles.glm = 1;
            case {'--nocache'}
                handles.nocache = 1;
            case {'-h','--help'}
                showUsage();
                return;
            case {'--version'}
                disp([exeName ' version : ' num2str(versionNumber)]);
                return;
            otherwise
                if strcmp(varargin{i}(1), '-')
                    disp(['bad option : ' varargin{i}]);
                    i = size(varargin, 2);
                    handles.commandError = 1;
                else
                    handles.matFiles = [handles.matFiles varargin{i}];
                end
        end
        i = i + 1;
    end

    % check command input
    if handles.commandError
        showUsage();
        return;
    elseif isempty(handles.cx)
        disp('no cells of subject time-series (CX) file. please specify .mat file.');
        showUsage();
        return;
    elseif isempty(handles.model) && isempty(handles.pymodel) 
        disp('no group surrogate model file. please specify .mat file.');
        showUsage();
        return;
    elseif isempty(handles.roisidx) && isempty(handles.targatl) 
        disp('no modulation target index or target atlas file. please specify .mat or nifti file.');
        showUsage();
        return;
    elseif ~isempty(handles.targatl) && isempty(handles.atlas)
        disp('when target atlas is specified, an atlas file must be specified.');
        showUsage();
        return;
    elseif handles.glm>0 && isempty(handles.atlas)
        disp('when GLM output is specified, an atlas file must be specified.');
        showUsage();
        return;
    end

    % process input files
    processInputFiles(handles);
end

%%
% show usage function
function showUsage()
    global exePath;
    global exeName;
    
    disp(['usage: ' exeName ' [options] [permfile.mat ...]']);
    disp('  -c, --cx name       set cells of subject time-series (<filename>.mat)');
    disp('  -m, --model name    set (VAR) group surrogate model (<filename>_gsm_var.mat)');
    disp('  -a, --atlas name    set cube atlas nifti file (<filename>.nii.gz)');
    disp('  -t, --targatl name  set modulation target atlas nifti file (<filename>.nii.gz)');
    disp('  -i, --roiidx name   set modulation target ROIidx file (<filename>.mat)');
    disp('  -r, --roi num       set modulation target ROI <num> or <range text>');
    disp('  -o, --out num       set output perm & surrogate files number <num> (default:1)');
    disp('  --surrnum num       output surrogate sessions per one file <num> (default:40)');
    disp('  --srframes num      output surrogate frames <num> (default:160)');
    disp('  --vnparam num num num  set virtual neuromodulation params <num num num> (default:28 22 0.15)');
    disp('  --tr num               set TR (second) of fMRI time-series <num> (default:1)');
    disp('  --hrfparam num num     set HRF (for convolution) params <num num> (default:16 8)');
    disp('  --glm               output GLM result nifti file.');
    disp('  --outpath path      output files <path> (default:"results")');
    disp('  --nocache           do not output surrogate file');
    disp('  --pymodel path      set (VAR) group surrogate model <path> by vneumodpy (Python)');
    disp('  --version           show version number');
    disp('  -h, --help          show command line help');
end

%%
% process input files (mail rutine)
%
function processInputFiles(handles)
    global exePath;
    global exeName;

    [path,savename,ext] = fileparts(handles.cx); % get savename
    disp(['set savename=' savename])

    % load modulation target ROI index
    dbsidxs = {};
    if ~isempty(handles.roisidx)
        disp(['load ROIidx file: ' handles.roisidx])
        load(handles.roisidx);
        if length(ROIidx) ~= size(CX{1},1)
            disp(['error: node number does not match. ROIidx(' num2str(length(ROIidx)) ') vs. CX(' num2str(size(CX{1},1)) ')']);
            return;
        end
        if isempty(handles.rois)
           rois = unique(ROIidx);
           rois(isnan(rois)) = [];
           rois(rois==0) = [];
           handles.rois = rois;
        end
        for i=1:length(handles.rois)
            dbsidxs{i} = find(ROIidx==handles.rois(i)); % find target index
        end
    end
    if ~isempty(handles.targatl)
        disp(['load target atlas file: ' handles.targatl])
        % load cube atlas file
        atlasinfo = niftiinfo(handles.atlas);
        atlasV = niftiread(atlasinfo);
        atlasV = adjustVolumeDir(atlasV, atlasinfo);
        % load modulation target atlas file
        atargatlinfo = niftiinfo(handles.targatl);
        targV = niftiread(atargatlinfo);
        targV = adjustVolumeDir(targV, atargatlinfo);
        if isempty(handles.rois)
           rois = unique(targV);
           rois(isnan(rois)) = [];
           rois(rois==0) = [];
           handles.rois = rois;
        end
        for i=1:length(handles.rois)
            idx = find(targV==handles.rois(i)); % find target index
            dbsidx = unique(atlasV(idx));
            dbsidx(isnan(dbsidx)) = [];
            dbsidxs{i} = dbsidx;
        end
    end
    if isempty(dbsidxs)
        disp(['error: empty modulation target. bad ROI num=' num2str(handles.rois)]);
        return;
    end
    
    % load cells of subject time-series file
    disp(['load subject time-series file: ' handles.cx])
    load(handles.cx);

    % load model file
    if ~isempty(handles.model)
        disp(['load model file: ' handles.model])
        load(handles.model);
    elseif ~isempty(handles.pymodel) 
        net = loadMvarNetworkPy(handles.pymodel);
    else
        disp('error: empty model file.');
        return;
    end

    % init
    vnpm = handles.vnparam;
    hrfpm = handles.hrfparam;
    isMatf = ~isempty(handles.matFiles);
    if isMatf
        N = length(handles.matFiles);
    else
        N = handles.out;
    end
    
    % process each file
    for i = 1:N
        perm = [];
        if isMatf
            % load permutation mat file
            fname = handles.matFiles{i};
            if ~exist(fname,'file')
                disp(['file is not found. ignoring : ' fname]);
                continue;
            end
            [path,name,ext] = fileparts(fname);
            if strcmp(ext,'.mat')
                f = load(fname);
                if isfield(f,'perm')
                    disp(['load perm file : ' fname]);
                    perm = f.perm;
                end
            end
            if isempty(perm)
                disp(['perm file is not found. ignoring : ' fname]);
                continue;
            end
        end
        permfname = [handles.outpath '/perm' num2str(i) '_' savename '.mat'];
        if exist(permfname,'file')
            disp(['load perm file : ' permfname]);
            load(permfname); % load prev result
        end
        if isempty(perm)
            % generating subject permutation
            [perm, uxtime] = getVnmSubjectPermutation(CX);
            save(permfname, 'perm', 'uxtime', '-v7.3');
            disp(['save perm file : ' permfname]);
        end

        % loop for neuromodulation target rois
        for j=1:length(handles.rois)
            S = [];
            roi = handles.rois(j);

            % get modulation (add & mul) time-series for vertual neuromodulation
            disp(['generate modulation (add & mul) time-series, target roi=' num2str(roi) ', srframes=' num2str(handles.srframes) ', dbsoffsec=' num2str(vnpm(1)) ', dbsonsec=' num2str(vnpm(2)) ', dbspw=' num2str(vnpm(3))]);
            disp(['convolution params tr=' num2str(handles.tr) ', res=' num2str(hrfpm(1)) ', sp=' num2str(hrfpm(2))]);
            [CA, Chrf, CM] = getVnmAddMulSignals(CX, dbsidxs{j}, handles.surrNum, handles.srframes, vnpm(1), vnpm(2), vnpm(3), handles.tr, hrfpm(1), hrfpm(2));

            % load or calc virtual neuromodulation surrogate
            sessionName = [savename '_' num2str(roi) 'sr' num2str(handles.surrNum) 'pr' num2str(i)];
            outfname = [handles.outpath '/' sessionName '.mat'];
            if exist(outfname,'file')
                disp(['load surrogate file : ' outfname]);
                load(outfname); % load prev result
            end
            if isempty(S)
                % calc virtual neuromodulation VAR surrogate
                disp(['calc virtual neuromodulation surrogate. roi=' num2str(roi)  ', surrNum=' num2str(handles.surrNum)]);
                S = calcVnmVarSurrogate(net, CX, CA, CM, perm, handles.surrNum, handles.srframes);
        
                % save VAR surrogate file
                if handles.nocache == 0
                    save(outfname,'S','-v7.3');
                    disp(['save virtual neuromodulation surrogate file : ' outfname]);
                end
            end

            % calcurate 2nd level GLM and save nifti file
            if handles.glm > 0
                % load atlas file
                atlasinfo = niftiinfo(handles.atlas);
                atlasV = niftiread(atlasinfo);
                atlasV = adjustVolumeDir(atlasV, atlasinfo); % this does not affect
        
                tuM = 8;  % GLM tukey-taper size
                betaBmat = [handles.outpath '/' sessionName '_2nd-Tukey' num2str(tuM) '.mat'];
                if exist(betaBmat,'file')
                    % load beta volumes
                    disp(['load 2nd level GLM result file : ' betaBmat]);
                    load(betaBmat);
                else
                    % calc 1st-level GLM
                    disp(['calc 1st-level GLM...']);
                    surrNum = length(S);
                    bmatC = cell(1,surrNum);
                    for k=1:surrNum
                        f = struct();
                        % calc 1st level GLM
                        Xorg = Chrf{k};
                        Xt = [Xorg, ones(size(Xorg,1),1)];
                        [f.B2, RSS, df] = calcGlmTukey(S{k}', Xt, tuM);
%                        [recel, f.FWHM] = estimateSmoothFWHM(R, RSS, df, estiV);
                        bmatC{k} = f;
                    end

                    % calc 2nd-level estimation
                    disp(['calc 2nd-level GLM...']);
                    B1 = [];
                    X2 = [];
                    FWHMs = [];
                    for k=1:surrNum
                        f = bmatC{k};
                        % 2nd-level Y vector
                        B2 = f.B2(:,[1,2]); % include design and intercept (we need more than 8 length for tukey taper)
                        B1 = [B1; B2'];
%                        FWHMs = [FWHMs; f.FWHM];

                        % 2nd-level design matrix
                        X2 = [X2; eye(size(B2,2))];
                    end
                    clear f;
                    B1(isnan(B1)) = 0; % there might be nan
%                    FWHMs = mean(FWHMs,1); % let's take the mean of FWHM.

                    % calc 2nd-level estimation
                    [B, RSS, df, X2is, tRs] = calcGlmTukey(B1, X2, tuM);
%                    [recel, FWHM] = estimateSmoothFWHM(R, RSS, df, estiV);
    
                    % output beta matrix
                    if handles.nocache == 0
                        save(betaBmat,'B','RSS','X2is','tRs','df','-v7.3');
                        disp(['save 2nd level GLM result file : ' betaBmat]);
                    end
                end
        
                % GLM contrast images
                contrasts = {[1 0]'}; % GLM contrust
                Ts = calcGlmContrastImage(contrasts, B, RSS, X2is, tRs);
                V2 = getNifti4DFromRoiTS(Ts{1}, atlasV);
        
                % save T-value NIfTI volume
                saveNifti(handles.atlas,V2,[handles.outpath '/'],[sessionName '_2nd-Tukey' num2str(tuM)]);
            end
        end
    end
end

%% should save one file
function saveNifti(tfmri, V2, path, outname)
    info = niftiinfo(tfmri);
    info.ImageSize = info.ImageSize(1:3);
    info.PixelDimensions = info.PixelDimensions(1:3);
    info.Datatype = 'single';
    info.BitsPerPixel = 32;
    fname = [path outname '.nii'];
    V = adjustVolumeDir(V2, info); 
    niftiwrite(V,fname,info,'Compressed',true);
    disp(['save nifti file : ' fname]);
end
