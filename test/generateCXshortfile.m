% generate CX short file of PD rs-fMRI cube ROI time-series
% generateSubjectTimeSeries.m, generateCXfile.m run first.

function generateCXshortfile
    smoothes = {'s24','s26','s28'}; %'s34'; % 34 for group. check 14 to 36 for individual.
    nuisance = 'gmacomp'; % 'aro'; % 
    atlasSizes = [3, 2];
    lags = [1,2,3,4,5];
    dtype = '1005'; %'0722'; %'0622'; % ''; % for pd1s 'hcp1s'; % 'hcp'; % 'pd'; % 'hc'; %  
    sbjmax = '';
    sessidx = [1]; %[1 2 3];
    mtype = ''; % whole brain 'Ecp'; % except cerebellum & pons 

    for sz=2
        atlasSize = atlasSizes(sz);
        for p=1:1
            for s=1:length(smoothes)
                generateCXfileAlgos(atlasSize, smoothes{s}, nuisance, dtype, sbjmax, sessidx, mtype);
            end
        end
    end
end

function generateCXfileAlgos(atlasSize, smooth, nuisance, dtype, sbjmax, sessidx, mtype)
    cubename = ['Allen' mtype 'Cube' num2str(atlasSize)];
    cxname = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurrCX' cubename smooth nuisance dtype sbjmax '.mat'];
    if ~exist(cxname,'file')
        disp(['file not found : ' cxname]);
        return;
    end
    load(cxname);

    str = [];
    for i=1:length(sessidx)
        str = [str num2str(sessidx(i))];
    end
    cxname = ['results/dbs' dtype num2str(atlasSize) '/testdbsSurrCX' cubename smooth nuisance dtype sbjmax str '.mat'];
    CX = CX(sessidx);
    [a,idx] = ismember(perm,sessidx);
    perm(idx==0) = [];
    save(cxname, 'CX', 'perm', '-v7.3');
end
