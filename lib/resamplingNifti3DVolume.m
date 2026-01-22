%%
% resampling NIfTI volume.
% returns re-sampled volume (outV)
% input:
%  V            nifti 3D volume (X x Y x Z)
%  stepXY       XY axes resampling rate
%  stepZ        Z axes resampling rate
%  operation    operation for each plane ('mode'(default),'max','min','mean','median')

function outV = resamplingNifti3DVolume(V, stepXY, stepZ, operation)
    if nargin < 4, operation = 'mode'; end

    xs = floor(size(V,1)/stepXY);
    ys = floor(size(V,2)/stepXY);
    zs = floor(size(V,3)/stepZ);
    outV = zeros(xs,ys,zs,'single');
    if stepXY >= 1 && stepZ >= 1
        % scale down
        for z=1:zs
            for y=1:ys
                for x=1:xs
                    A = V(round(x*stepXY-(stepXY-1)):round(x*stepXY), ...
                          round(y*stepXY-(stepXY-1)):round(y*stepXY), ...
                          round(z*stepZ-(stepZ-1)):round(z*stepZ));
                    switch(operation)
                    case 'mode'
                        m = mode(A(~isnan(A)),'All');
                    case 'min'
                        m = nanmin(A(:));
                    case 'max'
                        m = nanmax(A(:));
                    case 'mean'
                        m = nanmean(A(:));
                    case 'median'
                        m = nanmedian(A(:));
                    end
                    outV(x,y,z) = m;
                end
            end
        end
    else
        % scale up
        for z=1:size(V,3)
            for y=1:size(V,2)
                for x=1:size(V,1)
                    m = V(x,y,z);
                    xx = 1+round((x-1)/stepXY):round(x/stepXY);
                    yy = 1+round((y-1)/stepXY):round(y/stepXY);
                    zz = 1+round((z-1)/stepZ):round(z/stepZ);
                    outV(xx,yy,zz) = m;
                end
            end
        end
    end
end

