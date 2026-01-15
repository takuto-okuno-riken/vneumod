%%
% get NIfTI 4D volume from ROI time-series (matrix).
% returns nifti 4D volume (V)(X x Y x Z x frames)
% input:
%  X            ROI time-series (ROIs x frames)
%  atlasV       nifti 3D atlas (X x Y x Z)

function V = getNifti4DFromRoiTS(X, atlasV)
    roiMax = size(X,1);

    xs = size(atlasV,1);
    ys = size(atlasV,2);
    zs = size(atlasV,3);
    % reshape should work. but here not use it.
    % V = reshape(X, xs, ys, zs, size(X,2));
    V = single(nan(xs, ys, zs, size(X,2)));
    for i=1:roiMax
        idxs{i} = find(atlasV(:)==i);
    end
    for t=1:size(X,2)
        A = nan(xs, ys, zs);
        for i=1:roiMax
            A(idxs{i}) = X(i,t);
        end
        V(:,:,:,t) = A;
    end
end

