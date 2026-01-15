% cube Atlas (generated from resampling) ROIs for DBS simulation
% based on itksnap_annotation_full2mmMask.nii.gz
% (testDBS5.m)

function makeCubeAtlas
%%{
    % whole brain
    name = 'allen';
    idxname = 'allen';%'brodmann'; %
    originalMask = 'data/itksnap_annotation_full2mmMask.nii.gz';
%}
%{
    % whole brain except cerebellum & pons
    name = 'allenEcp';
    idxname = 'allenEcp';%'brodmannEcp'; %
    originalMask = 'data/itksnap_annotation_full2mmMaskEcp.nii.gz';
%}
    % generate atlas of cube ROIs
    for atlasSize = 2 %10:-1:2
        cubename = [name 'Cube' num2str(atlasSize)];
        atlas = ['data/' cubename 'atlas.nii' ];
        if exist([atlas '.gz'],'file')
            atlasinfo = niftiinfo([atlas '.gz']);
            atlasV = niftiread(atlasinfo);
            atlasV = adjustVolumeDir(atlasV, atlasinfo);
        else
            % load mask
            info = niftiinfo(originalMask);
            maskV = niftiread(info);
            maskV = adjustVolumeDir(maskV, info); % need to adjust direction for cubeROI.
            maskV(maskV>0) = 1;
            % resampling mask size to 1/atlasSize
            maskV2 = resamplingNifti3DVolume(maskV, atlasSize, atlasSize, 'max');
            idx = find(maskV2(:)==1);
            maskV2(idx) = 1:length(idx);
            V = resamplingNifti3DVolume(maskV2, 1/atlasSize, 1/atlasSize);
            V(91,109,91) = 0;
            atlasV2 = V .* single(maskV);

            % use info as template for cubeAtlas, but direction is opposite.
            if max(atlasV2(:)) > 65535
                atlasV = int32(atlasV2);
                info.Datatype = 'int32';
                info.BitsPerPixel = 32;
            else
                atlasV = uint16(atlasV2);
                info.Datatype = 'uint16';
                info.BitsPerPixel = 16;
            end
            % set info. info.raw is not necessary to set (niftiwrite() does it)
            info.Description = 'Cube ROI';
            info.TransformName = 'Sform';
            info.PixelDimensions(1:3) = 2;
            info.Transform.T(1:3,1:3) = eye(3) * 2; % Sform. left to right (same as volume)
            info.Transform.T(4,1:3) = [-90 -126 -72]; % flip x offset
            % output nii file
            aV = adjustVolumeDir(atlasV, info); % this does not affect volume
            niftiwrite(aV,atlas,info,'Compressed',true);
        end

        % get downsampled ROI index
        ROIidxf = ['data/' idxname 'ROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            if strcmp(idxname,'brodmann') || strcmp(idxname,'brodmannEcp')
                originalAtlas = ['data/brodmannConnfMRI.nii.gz'];
            else
                originalAtlas = ['data/itksnap_annotation_full2mm.nii.gz'];
            end
            roiinfo = niftiinfo(originalAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV==0) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            % find left or right
            info = niftiinfo(originalMask);
            maskV = niftiread(info);
            maskV = adjustVolumeDir(maskV, info); % need to adjust direction for cubeROI.
            rV = maskV; rV(1:45,:,:) = 1; rV(46:end,:,:) = 2; % left to right (same as volume)
            maskV = maskV .* rV;
            maskV(maskV==0) = nan;
            SDidx = getRoiTSFromNifti4D(maskV, atlasV, 'mode');
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%{
        ROIidxf = ['data/' name 'StnvimppnROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            vimppnAtlas = ['data/stnvimppn2mm.nii.gz'];
            roiinfo = niftiinfo(vimppnAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV==0) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%}
%{
        ROIidxf = ['data/' name 'PpnmotROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            vimppnAtlas = ['D:\work\defaultmode\human\VIMPPN/ppnmot2mm.nii.gz'];
            roiinfo = niftiinfo(vimppnAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV<10000) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%}
%{
        ROIidxf = ['data/' name 'PpnassoROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            vimppnAtlas = ['D:\work\defaultmode\human\VIMPPN/ppnasso2mm.nii.gz'];
            roiinfo = niftiinfo(vimppnAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV<1000) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%}
%{
        ROIidxf = ['data/' name 'VopvimROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            vimppnAtlas = ['D:\work\defaultmode\human\VIMPPN/vopvim2mm.nii.gz'];
            roiinfo = niftiinfo(vimppnAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV<1000) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%}
%{
        ROIidxf = ['data/' name 'StnlimassoROIindex' num2str(atlasSize) '.mat'];
        if exist(ROIidxf,'file')
            load(ROIidxf);
        else
            vimppnAtlas = ['D:\work\defaultmode\human\VIMPPN/stnlimasso2mm.nii.gz'];
            roiinfo = niftiinfo(vimppnAtlas);
            roiV = single(niftiread(roiinfo));
            roiV = adjustVolumeDir(roiV, roiinfo); % need to adjust direction for cubeROI.
            roiV(roiV<1000) = nan;
            ROIidx = getRoiTSFromNifti4D(roiV, atlasV, 'mode');
            [RBs, RIs] = sort(ROIidx);
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        end
%}
        % allen DBS 4 x 4 x 8 mm atlas
%{
        % load mask
        info = niftiinfo(originalMask);
        maskV = niftiread(info);
        maskV = adjustVolumeDir(maskV, info); % need to adjust direction for cubeROI.
        maskV(maskV>0) = 1;
        maskV = circshift(maskV,-1,3);
        % resampling mask size to 1/atlasSize
        maskV2 = resamplingNifti3DVolume(maskV, atlasSize, atlasSize*2, 'max');
        idx = find(maskV2(:)==1);
        maskV2(idx) = 1:length(idx);
        V = resamplingNifti3DVolume(maskV2, 1/atlasSize, 1/(atlasSize*2));

        V(91,109,91) = 0;
        atlasV2 = V .* single(maskV);
        atlasV2 = circshift(atlasV2,1,3);
        % set info. info.raw is not necessary to set (niftiwrite() does it)
        info.Datatype = 'uint16';
        info.BitsPerPixel = 16;
        info.Description = 'Cube ROI';
        info.TransformName = 'Sform';
        info.PixelDimensions(1:3) = 2;
        info.Transform.T(1:3,1:3) = eye(3) * 2; % Sform. left to right (same as volume)
        info.Transform.T(4,1:3) = [-90 -126 -72]; % flip x offset
        % output nii file
        aV = adjustVolumeDir(atlasV2, info); % this does not affect volume
        niftiwrite(uint16(aV),['data/allenDbs448mmAtlas.nii' ],info,'Compressed',true);
%}
        if strcmp(idxname,'allen') && atlasSize == 2
%{
            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) 'Stn.mat'];
            ROIidx([9603,8515,9598,8510]) = 4501; % anterior
            ROIidx([9574,9567]) = 4502; % center lateral
            ROIidx([9573,8480,9568,8475]) = 4503; % center
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
            aV = getNifti4DFromRoiTS(ROIidx,atlasV);
            aV = uint16(adjustVolumeDir(aV, atlasinfo)); % this does not affect volume
            niftiwrite(aV,['data/' cubename 'atlasStn.nii'],atlasinfo,'Compressed',true);

            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) 'Stn2.mat'];
            ROIidx([8515,8510]) = 4511; % anterior STN
            ROIidx([9603,9598]) = 4512; % anterior STN
            ROIidx([8480,8475]) = 4513; % center STN
            ROIidx([9573,9568]) = 4514; % center STN
            ROIidx([9574,9567]) = 4515; % posterior STN
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
            aV = getNifti4DFromRoiTS(ROIidx,atlasV);
            aV = uint16(adjustVolumeDir(aV, atlasinfo)); % this does not affect volume
            niftiwrite(aV,['data/' cubename 'atlasStn2.nii'],atlasinfo,'Compressed',true);
%}
%{
            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) 'Stn3.mat'];
            ROIidx([8514, 8511]) = 4521; % most anterior (ventral)
            ROIidx([8515, 8510, 9603, 9598]) = 4522; % anterior lateral
            ROIidx([8479, 8476]) = 4523; % center medial
            ROIidx([8480, 8475, 9573, 9568]) = 4524; % center 
            ROIidx([9574, 9567]) = 4525; % center lateral
            ROIidx([8448, 8443, 9544]) = 4526; % posterior medial
            ROIidx([9538, 9545]) = 4527; % most posterior (dorsal)
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
            aV = getNifti4DFromRoiTS(ROIidx,atlasV);
            aV = uint16(adjustVolumeDir(aV, atlasinfo)); % this does not affect volume
            niftiwrite(aV,['data/' cubename 'atlasStn3.nii'],atlasinfo,'Compressed',true);
%}
%{
            pinfo = niftiinfo('data/stnvimppn2mm.nii.gz'); % for PPN, Vim, Vim bottom
%            pinfo = niftiinfo('data/itksnap_annotation_full2mm.nii.gz'); % for GPi
            pV = niftiread(pinfo);
            pV = adjustVolumeDir(pV, pinfo);
%            ppidx = unique(atlasV(pV==200)); % PPN voxels
%            baseroi = 20011; partname = 'Ppn';
%            ppidx = unique(atlasV(pV==12)); % GPi voxels
%            baseroi = 1221; partname = 'GPi';
%            ppidx = unique(atlasV(pV==3501|pV==3502)); % Vim/Vop voxels
%            baseroi = 3511; partname = 'VL';
            ppidx = unique(atlasV(pV==210)); % Vim bottom voxels
            baseroi = 21001; partname = 'PSA';

            ppidx(ppidx==0) = []; % remove 0
            % find one x line
            loc = nan(length(ppidx),3);
            for i=1:length(ppidx)
                aidx = find(atlasV==ppidx(i));
                [x,y,z] = ind2sub(size(atlasV),aidx);
                D(1)=max(x); D(2)=max(y); D(3)=max(z);
                D = D + mod(D,2); % to odd position.
                loc(i,:) = [D(1), (110-D(2))*1000+D(3), ppidx(i)]; % flip Y for anteri to posteri
            end
            G = unique(loc(:,2));
            for i=1:length(G) % anterior to posterior, ventral to dorsal
                K = find(loc(:,2)==G(i));
                X = loc(K,1);
                X = abs(X - 47);
                % find left & right voxel pair, but sometimes not
                for j=1:2:45
                    M = (j==X);
                    if ~any(M), continue; end
                    ROIidx(loc(K(M),3)) = baseroi; % medial to lateral.
                    baseroi = baseroi + 1;
                end
            end
            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) partname '.mat'];
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
            aV = getNifti4DFromRoiTS(ROIidx,atlasV);
            aV = uint16(adjustVolumeDir(aV, atlasinfo)); % this does not affect volume
            niftiwrite(aV,['data/' cubename 'atlas' partname '.nii'],atlasinfo,'Compressed',true);
%}
        elseif strcmp(idxname,'allenEcp') && atlasSize == 2 % STN (11 voxels)
            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) 'Stn.mat'];
            ROIidx([6055,4996,6050,4991]) = 4501; % anterior STN
            ROIidx([6026,6019]) = 4502; % posterior STN
            ROIidx([6025,4963,6020,4958]) = 4503; % other STN
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
        elseif strcmp(idxname,'allen') && atlasSize == 3 % STN (one voxel), GPi (3 voxels)
            ROIidxf = ['data/' name 'ROIindex' num2str(atlasSize) 'GPi.mat'];
            ROIidx([3590,3597]) = 1201; % posterior GPi
            ROIidx([3614,3620]) = 1202; % anterior GPi
            ROIidx([3615,3619]) = 1203; % other GPi
            save(ROIidxf,'ROIidx','RBs','RIs','roiinfo','SDidx');
            aV = getNifti4DFromRoiTS(ROIidx,atlasV);
            aV = uint16(adjustVolumeDir(aV, atlasinfo)); % this does not affect volume
            niftiwrite(aV,['data/' cubename 'atlasGPi.nii'],atlasinfo,'Compressed',true);
        end
%}
        % get downsampled Network index
        NetIdxf = ['data/' idxname 'NetIndex' num2str(atlasSize) '.mat'];
        if exist(NetIdxf,'file')
            load(NetIdxf);
        else
            % default mode network (generated by testGLM6restC.m and some ROIs are manually fixed)
            originalAtlas = 'D:\work\spm12\conn\rois\networks.nii'; % Sform. right to left (flip x axis).
            netinfo = niftiinfo(originalAtlas);
            netV = single(niftiread(netinfo));
            netV = adjustVolumeDir(netV, netinfo); % need to adjust direction for cubeROI.
            netV(netV==0) = nan;
            NetIdx = getRoiTSFromNifti4D(netV, atlasV, 'mode');
            NetIdx2 = zeros(size(NetIdx,1),1);

            dmnf = 'data/humanAtlasDMN2a2.nii'; % Sform. right to left (flip x axis).
            dmninfo = niftiinfo([dmnf '.gz']);
            V3 = single(niftiread(dmninfo));
            V3(V3==0) = nan;
            V3 = adjustVolumeDir(V3, dmninfo); % need to adjust direction for cubeROI.
            V3i = getRoiTSFromNifti4D(V3, atlasV, 'mode');
            for i = 1:8
                idx = find(V3i==i);
                NetIdx2(idx) = i;
            end
            [NBs, NIs] = sort(NetIdx2);
            save(NetIdxf,'NetIdx','NetIdx2','NBs','NIs','netinfo');
        end

        % output ROI centroid for unity
%{
        roiNum = max(atlasV(:));
        x = []; y = []; z = [];
        for i=1:roiNum
            BW = logical(zeros(size(atlasV,1),size(atlasV,2),size(atlasV,3)));
            idx = find(atlasV==i);
            BW(idx) = 1;
            s = regionprops3(BW, 'Centroid');
            x = [x; s.Centroid(1,1)];
            y = [y; s.Centroid(1,2)];
            z = [z; s.Centroid(1,3)];
        end
        centroid = ['data/' cubename 'atlasCentroid.csv' ];
        writematrix([x y z], centroid);
%}
    end
end
