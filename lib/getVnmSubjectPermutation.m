%%
% Generate subject permutation time-series for virtual neuromodulation surrogate
% returns permutated time-series (perm), uxtime and residual length (reslen).
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  lags            number of lags for autoregression

function [perm, uxtime, reslen] = getVnmSubjectPermutation(CX, lags)
    perm = [];
    cxlen = length(CX);
    reslen = size(CX{1},2) - lags; % residual length

    % ordered residual with subject permutation
    uxtime = uint32(posixtime(datetime('now')));
    rng(uxtime);
    rp = randperm(cxlen);
    for i=1:cxlen
        perm = [perm, (1:reslen) + (rp(i)-1)*reslen];
    end
end
