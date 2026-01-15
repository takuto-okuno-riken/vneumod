%%
% Calculate virtual neuromodulation surrogate by Vector Auto-Regression 
% returns modulated surrogate time-series (S).
% input:
%  net             mVAR network (struct)
%  CX              cells of multivariate time series matrix {node x time series}
%  CA              cells of Addition time-series {node x time series}
%  CM              cells of Multiplication time-series {node x time series}
%  perm            subject permutation for ordered residual
%  surrNum         output number of surrogate samples
%  srframes        frame length of surrogate time-series

function S = calcVnmVarSurrogate(net, CX, CA, CM, perm, surrNum, srframes)
    dist = 'residuals';
    usegpu = false;

    cxlen = length(CX);
    frames = size(CX{1},2);

    % Virtual Neuromodulation VAR surrogate
    S = cell(surrNum,1);
    C = []; Err = [];
    for i=1:surrNum
        X = [];
        for k=0:cxlen-1
            X=[X, CX{mod((i-1)*2+k,cxlen)+1}];
            if size(X,2) >= srframes, break; end
        end
        if i > length(S) || isempty(S{i})
            tc = tic;
            % ordered residual with subject permutation
            nBaset = {perm(frames*(i-1)+1:frames*i), 0};

            [S{i}, C, Err] = surrogateDbsMVAR(X(:,1:srframes), [], [], [], net, CA{i}, CM{i}, dist, 1, NaN, nBaset, C, Err, usegpu);
            disp(['done t=' num2str(toc(tc)) 'sec']);
        end
    end
end
