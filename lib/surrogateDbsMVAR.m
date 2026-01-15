%%
% Simulated DBS surrogate multivariate signal generation by multivariate VAR
% based on autoregressive (AR) surrogates (R. Liegeois et al., 2017)
% returns surrogated signals (Y)
% input:
%  X            multivariate time series matrix (node x time series)
%  exSignal     multivariate time series matrix (exogenous input x time series) (optional)
%  nodeControl  node control matrix (node x node) (optional)
%  exControl    exogenous input control matrix for each node (node x exogenous input) (optional)
%  net          mVAR network
%  A            DBS stimulus multivariate time series matrix (add) (node x time series)
%  M            DBS stimulus multivariate time series matrix (multi) (node x time series)
%  dist         distribution of noise to yield surrogates ('gaussian'(default), 'residuals')
%  surrNum      output number of surrogate samples (default:1)
%  yRange       range of Y value (default:[Xmin-Xrange/5, Xmax+Xrange/5])
%  nBaset       noise base time (default:empty (randperm))
%  Cin          coefficient matrix of VAR (default:[])
%  Errin        Err matrix of VAR (default:[])
%  usegpu       use gpu for regress function (default:false)

function [Y, C, Err, perm] = surrogateDbsMVAR(X, exSignal, nodeControl, exControl, net, A, M, dist, surrNum, yRange, nBaset, Cin, Errin, usegpu)
    if nargin < 12, usegpu = false; end
    if nargin < 11, Errin = []; end
    if nargin < 10, Cin = []; end
    if nargin < 9, nBaset = []; end
    if nargin < 8, yRange = NaN; end
    if nargin < 7, surrNum = 1; end
    if nargin < 6, dist = 'gaussian'; end

    nodeNum = size(X,1);
    sigLen = size(X,2);
    exNum = size(exSignal,1);
    inputNum = nodeNum + exNum;
    lags = net.lags;

    % set node input
    Xorg = [X; exSignal];

    % set control 3D matrix (node x node x lags)
    [~,~,control] = getControl3DMatrix(nodeControl, exControl, nodeNum, exNum, lags);

    % calc Y range
    if isnan(yRange)
        t = max(X(:)); d = min(X(:));
        r = t - d;
        yRange = [d-r/5, t+r/5];
    end

    rvlen = length(net.rvec{1});
    if isempty(Errin)
        for i=2:nodeNum, if rvlen > length(net.rvec{i}), rvlen = length(net.rvec{i}); end; end
        Err = single(nan(nodeNum,rvlen));
        for i=1:nodeNum
            Err(i,:) = net.rvec{i}(1:rvlen);
        end
    else
        Err = Errin;
    end
    % get coefficient matrix
    if isempty(Cin)
        C = single(zeros(nodeNum,inputNum*lags+1));
        for i=1:nodeNum
            idx = find(control(i,:,:)==1);
            C(i,[idx(:).' end]) = net.bvec{i};
        end
    else
        C = Cin;
    end

    if strcmp(dist,'gaussian')
        P  = mean(Err.');
        EC = cov(Err.',1);
        Err = (mvnrnd(P,EC,size(Err,2)))';
    end
    noise = Err;

    S2 = ones(inputNum*lags+1,1);
    % use gpu array
    if usegpu
        % check max grid size
        maxGridSize = gpuDevice().MaxGridSize(1);
        if size(C,1)*size(C,2) < maxGridSize, C = gpuArray(C); end
        if size(noise,1)*size(noise,2) < maxGridSize, noise = gpuArray(noise); end
        S2 = gpuArray(S2);
    end

    % set noise permutation
    if ~isempty(nBaset)
        if iscell(nBaset)
            perm = nBaset{1};
            perm = perm + nBaset{2};
        elseif length(nBaset)==1
            perm = nBaset + 1:(size(noise,2)-lags);
        else
            perm = nBaset;
        end
    else
        rng('shuffle');
        perm = randperm(size(noise,2)-lags); % full residual/gaussian error permutation. this one is too unstable.
    end

    Y = nan(nodeNum,sigLen,surrNum);
    for k=1:surrNum
        disp(['surrogate sample : ' num2str(k)]);
        S = single(Xorg);
        % random is not so robust. use fixed end frame signals
        S(1:nodeNum,1:lags) = S(1:nodeNum,sigLen-lags+1:sigLen); % Initialization of the AR surrogate
        perm2 = circshift(perm,-sigLen*(k-1),2); % circle shift perm for each surrogate

        for t=lags+1:sigLen
            T = S(:,t); % next output

            for p=1:lags
                S2(1+inputNum*(p-1):inputNum*p) = S(:,t-p);
            end
            T(1:nodeNum) = C * S2 + noise(:,perm2(t-lags));

            % DBS stimulation
            if ~isempty(M)
                Mt = M(:,t);
                idx = find(~isnan(Mt));
                if ~isempty(idx)
                    T(idx) = T(idx) .* Mt(idx);
                end
            end
            if ~isempty(A)
                At = A(:,t);
                idx = find(~isnan(At));
                if ~isempty(idx)
                    T(idx) = T(idx) + At(idx);
                end
            end

            % fixed over shoot values
            if ~isempty(yRange)
                T(T < yRange(1)) = yRange(1);
                T(T > yRange(2)) = yRange(2);
            end
            S(:,t) = T;
        end
        Y(:,:,k) = S(1:nodeNum,:);
    end
end
