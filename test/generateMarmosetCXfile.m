% generate CX file of PD rs-fMRI cube ROI time-series
% generateSubjectTimeSeries.m run first.

function generateMarmosetCXfile
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

    cubename = ['BMA2' mtype 'marmo' num2str(atlasSize)];
    % calc var net values (regression)
    if lag==1
        lagStr = '';
    else
        lagStr = ['L' num2str(lag)];
    end

    cxname = ['results/pd36' cubename smooth nuisance dtype sbjmax '.mat'];
    if ~exist(cxname,'file')
        CX = {};
        listing = dir(['results/ts' dtype num2str(atlasSize) '/x' 'BMA2' mtype 'marmo' num2str(atlasSize) smooth nuisance '*.mat']);
        perm = randperm(length(listing));
        for i=1:length(listing)
            if ~isempty(sbjmax) && str2double(sbjmax) < i, continue; end
            j = perm(i);
            matf = ['results/ts' dtype num2str(atlasSize) '/' listing(j).name];
            load(matf);
            if ~isempty(mIdx)
                X = X(:,mIdx);
            end
            xs = std(X(:),1);
            if ~strcmp(dtype, 'hcp1s') && xs > 10, disp('bad std X'); end
%            X = X / xs; % normalize in each subject (not necessary)

            % show histogram 
%            figure; histogram(X); title(['histogram of ' matf]);
            CX{i} = X.';
        end
        save(cxname, 'CX', 'perm', '-v7.3');
    end
end
