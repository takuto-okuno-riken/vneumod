% generate CX file of PD rs-fMRI cube ROI time-series
% generateSubjectTimeSeries.m run first.

function generateCXfile
    algos = {'var','pc80','pcvar','pc999','pc9999','rdg01','rdg05','rdg800','rdg8000','las02','las05','las08','mkvar','mk50','mk80'}; % 'sigvar','sigmvar', % sigvar not work. sigmvar, soso.
    smooth = 's34';
    nuisance = 'gmacomp'; % 'aro'; % 
    atlasSizes = [3, 2];
    lags = [1,2,3,4,5];
    dtype = ''; % for pd1s 'hcp1s'; % 'hcp'; % 'pd'; % 'hc'; %  
    sbjmax = ''; % for others '61'; %'30'; % for pd30 
    mtype = ''; % whole brain 'Ecp'; % except cerebellum & pons 

    for a=1
        algo = algos{a};
        for sz=2
            atlasSize = atlasSizes(sz);
            for p=1:1
                generateCXfileAlgos(algo, atlasSize, lags(p), smooth, nuisance, dtype, sbjmax, mtype);
            end
        end
    end
end

function generateCXfileAlgos(algo, atlasSize, lag, smooth, nuisance, dtype, sbjmax, mtype)
    mIdx = [];
    if strcmp(mtype,'Ecp')
        % need to run testDBS5.m first
        ROIidxf = ['data/allen' mtype 'ROIindex' num2str(atlasSize) '.mat'];
        load(ROIidxf);

        if isempty(mIdx)
            ainfo = niftiinfo(['data/allen' 'Cube' num2str(atlasSize) 'atlas.nii.gz']);
            aV = niftiread(ainfo);
            einfo = niftiinfo(['data/allen' mtype 'Cube' num2str(atlasSize) 'atlas.nii.gz']);
            eV = niftiread(einfo);
            evmax = max(eV(:));
            mIdx = nan(1,evmax);
            for i=1:evmax
                A = aV(eV==i);
                mIdx(i) = mode(A);
            end
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx','mIdx');
        end
    end

    cubename = ['Allen' mtype 'Cube' num2str(atlasSize)];
    % calc var net values (regression)
    if lag==1
        lagStr = '';
    else
        lagStr = ['L' num2str(lag)];
    end

    cxname = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurrCX' cubename smooth nuisance dtype sbjmax '.mat'];
    if ~exist(cxname,'file')
        CX = {};
        if strcmp(dtype, 'hcp')
            listing = dir(['results/glm2/testglm6' 'AllenCube' num2str(atlasSize) 'rest' smooth nuisance '*.mat']);
            perm = randperm(length(listing));
            for i=1:length(listing)
                if ~isempty(sbjmax) && str2double(sbjmax) < i, continue; end
                j = perm(i);
                load(['results/glm2/' listing(j).name]);
                CX{i} = X.';
            end
        else
            listing = dir(['results/dbs' dtype num2str(atlasSize) '/testdbs' 'AllenCube' num2str(atlasSize) smooth nuisance '*.mat']);
            perm = randperm(length(listing));
            for i=1:length(listing)
                if ~isempty(sbjmax) && str2double(sbjmax) < i, continue; end
                j = perm(i);
                matf = ['results/dbs' dtype num2str(atlasSize) '/' listing(j).name];
                load(matf);
                if ~isempty(mIdx)
                    X = X(:,mIdx);
                end
                xs = std(X(:),1);
                if ~strcmp(dtype, 'hcp1s') && xs > 10, disp('bad std X'); end
                if strcmp(extractBefore(algo,4),'sig')
                    X = convert2SigmoidSignal(X'); % normalized [0 1] range by each subject.
                    if strcmp(extractBefore(algo,5),'sigm')
                        X = X - 0.5; % normalized [-0.5 0.5] range
                    end
                    X = X';
                else
%                    X = X / xs; % normalize in each subject (not necessary)
                end
                % show histogram 
%                figure; histogram(X); title(['histogram of ' matf]);
                CX{i} = X.';
            end
        end
        save(cxname, 'CX', 'perm', '-v7.3');
    end
end
