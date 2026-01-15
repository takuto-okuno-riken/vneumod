%%
% Get Addition and Multiplication time-series for virtual neuromodulation
% returns cells of Addition time-series (CA), cells of HRF time-series (Chrf), and cells of Multiplication time-series (CM).
% input:
%  CX              cells of multivariate time series matrix {node x time series}
%  dbsidx          target modulation ROI index (for CX)
%  surrNum         output number of surrogate samples
%  srframes        frame length of surrogate time-series
%  dbsoffsec       neuromodulation off duration (seconds)
%  dbsonsec        neuromodulation on duration (seconds)
%  dbspw           neuromodulation power
%  TR              TR of fMRI data (CX)
%  res             HRF sampling resolution (optional)
%  sp              HRF sampling starting point (optional)

function [CA, Chrf, CM] = getVnmAddMulSignals(CX, dbsidx, surrNum, srframes, dbsoffsec, dbsonsec, dbspw, TR, res, sp)
    if nargin < 10, sp = 8; end  % HRF sampling starting point
    if nargin < 9, res = 16; end  % HRF sampling resolution
    if nargin < 8, TR = 1.0; end
    if nargin < 7, dbspw = 0.15; end 
    if nargin < 6, dbsonsec = 22; end 
    if nargin < 5, dbsoffsec = 28; end 
    if nargin < 4, srframes = 160; end 
    if nargin < 3, surrNum = 40; end 

    % get HRF
    dt = TR / res;
    [t, hrf] = getGlmHRF(dt); % human's HRF;

    % find time series range
    cxlen = length(CX);
    X3 = nan(size(CX{1},1),size(CX{1},2),cxlen,'single');
    for nn=1:cxlen
        X = CX{nn};
        X3(1:size(X,1),1:size(X,2),nn) = X;
    end
    frames = size(CX{1},2);
    m = nanmean(X3(:));
    s = nanstd(X3(:),1);
    clear X3;

    % DBS block design (off, on, off, on, ...)
    CA = cell(1,surrNum); Chrf = cell(1,surrNum); CM = cell(1,surrNum); 
    for i=1:surrNum
        n = size(CX{i},1);
        bmax = floor((srframes * TR) / (dbsoffsec+dbsonsec));
        ons = zeros(bmax,1); dur = zeros(bmax,1);
        for j=0:bmax-1
            ons(j+1) = dbsoffsec + j*(dbsoffsec+dbsonsec);
            dur(j+1) = dbsonsec;
        end
        onsets{1} = ons;
        durations{1} = dur;
        [Chrf{i}, U] = getGlmHRFDesignMatrix(onsets, durations, srframes, TR, res, sp, hrf);
        block = nan(n,srframes,'single');
        block(dbsidx,:) = repmat(Chrf{i}' * dbspw * s, length(dbsidx), 1);
        CA{i} = block;
        block(dbsidx,:) = repmat(1 - 0.5 * Chrf{i}', length(dbsidx), 1);
        CM{i} = block;
    end
end
