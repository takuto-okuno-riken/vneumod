%%
% Generate PCA components time-series by PCA with Cell inputs of X and Ex.
% Time lags is fixed to 1.
% returns PCA components time-series and max component number for explainedTh.
% coeff and pcmu are PCA coefficient and mean values.
% input:
%  CX              cells of multivariate time-series matrix {node x time-series}
%  CexSignal       cells of multivariate time-series matrix {exogenous input x time-series} (default:{})
%  explainedTh     explained threshold for PCA components (default:1.0)
%  verbose         show verbose log (default:false)
%  maxSigLen       maxmum time-series length for input CX (default:[])
%  cachename       unique cache name (default:[])
%  cachepath       cache path (default:'results/cache')

function [CX, maxComp, coeff, pcmu] = generatePCAcompTimeSeries(CX, CexSignal, explainedTh, verbose, maxSigLen, cachename, cachepath)
    if nargin < 7, cachepath = 'results/cache'; end
    if nargin < 6, cachename = []; end
    if nargin < 5, maxSigLen = []; end
    if nargin < 4, verbose = false; end
    if nargin < 3, explainedTh = 1.0; end
    if nargin < 2, CexSignal = {}; end
    cxNum = length(CX);
    nodeNum = size(CX{1},1);
    if ~isempty(CexSignal)
        exNum = size(CexSignal{1},1);
    else
        exNum = 0;
    end
    if ~isempty(maxSigLen)
        for i=1:cxNum
            CX{i} = CX{i}(:,1:maxSigLen);
        end
        for i=1:exNum
            CexSignal{i} = CexSignal{i}(:,1:maxSigLen);
        end
    end
    lags = 1;
    expTh = explainedTh * 100;

    % set vector auto-regression (VAR) inputs
    if verbose, disp('set VAR inputs'); end
    allInLen = 0;
    for i=1:cxNum
        allInLen = allInLen + size(CX{i},2) - lags;
    end
%    Y = single(nan(allInLen,nodeNum));
    Xti = single(nan(allInLen,nodeNum*lags));
    Xe = single(nan(allInLen,exNum));
    cxRange = cell(1,cxNum);
    xts = 1;
    for i=1:cxNum
        % set node input
        X = flipud(CX{i}.'); % need to flip signal

        sLen = size(X,1);
        sl = sLen-lags;
%        Y(xts:xts+sl-1,:) = X(1:sl,:);
        Xti(xts:xts+sl-1,:) = X(2:sl+1,:);
        if exNum > 0
            E = flipud(CexSignal{i}.');
            Xe(xts:xts+sl-1,:) = E(1:sl,:);
        end
        cxRange{i} = [xts, xts+sl-1];
        xts = xts + sl;
    end
    clear CX;
    clear X; % memory clear
    clear Xt; % memory clear

    % apply the Principal Component Regress function
    if verbose, disp('apply PCA'); end
    pcacache = [cachepath '/pca-cache-' cachename '-' num2str(size(Xti,1)) 'x' num2str(size(Xti,2)) '-full.mat'];
    if exist(pcacache,'file') 
        load(pcacache);
    else
        Xti = double(Xti);
        [coeff,score,~,~,explained,pcmu] = pca(Xti); % relation : Xti == score * coeff.' + repmat(mu,size(score,1),1); time x node
%{      
        Xti = tall(double(Xti)); % no memory
        [~,score,~,~,explained] = pca(Xti); % relation : Xti == score * coeff.' + repmat(mu,size(score,1),1);
%}
%{
        numObs = 10000;
        [coeff,~,latent,~,~,mu]=pca(Xti(1:numObs,:),Centered=true,NumComponents=10000);
        IncrementalMdl = incrementalPCA(Coefficients=coeff, Latent=latent, Means=mu, NumObservations=numObs);
        % relation : Xti == score{i} * coeff{i}.' + repmat(mu{i},size(score{i},1),1);
        n = numel(Xti(numObs+1:end,1));
        nchunk = floor(n/numObs);
        topEV = zeros(nchunk,1);
        
        % Incremental fitting
        for j = 1:nchunk
            ibegin = min(n,1000 + numObs*(j-1) + 1);
            iend = min(n,1000 + numObs*j);
            IncrementalMdl = fit(IncrementalMdl,Xti(ibegin:iend,:));
            topEV(j) = IncrementalMdl.ExplainedVariance(1);
        end
%}
%        Xti = single(gather(Xti)); % no memory
        clear Xti;
        if ~isempty(cachename)
            save(pcacache,'coeff','score','explained','pcmu','-v7.3');
        end
%{
        if ~exist(icoeffcache,'file') 
            icoeff = invQR(coeff');
            save(icoeffcache,'icoeff','-v7.3');
        end
%}
    end

    % find 99% component range
    expTotal = 0;
    maxComp = size(score,2);
    for j=1:length(explained)
        expTotal = expTotal + explained(j);
        if expTotal >= expTh
            maxComp = j;
            break;
        end
    end

    CX = cell(1,cxNum);
    for i=1:cxNum
        r = cxRange{i};
        CX{i} = score(r(2):-1:r(1),1:maxComp)'; % reduce memory. !caution! score (Xti) has flipped signal. need to be opposit.
    end
end
