%%
% Generate subject permutation time-series for virtual neuromodulation surrogate
% returns permutated time-series (perm).
% input:
%  CX              cells of multivariate time series matrix {node x time series}

function [perm, uxtime] = getVnmSubjectPermutation(CX)
    perm = [];
    cxlen = length(CX);
    frames = size(CX{1},2);

    % ordered residual with subject permutation
    uxtime = uint32(posixtime(datetime('now')));
    rng(uxtime);
    rp = randperm(cxlen);
    for i=1:cxlen
        perm = [perm, (1:frames) + (rp(i)-1)*frames];
    end
end
