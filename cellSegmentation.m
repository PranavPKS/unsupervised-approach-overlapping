function cellSegmentationResult = cellSegmentation(im)
%% initializing parameters
param.timestep=3;  % time step
param.mu=.25/param.timestep;  % coefficient of the distance regularization term R(phi)
param.iter_inner=5;
param.iter_outer=20;
param.lambda=1; % coefficient of the weighted length term L(phi)
param.alfa=8;  % coefficient of the weighted area term A(phi)
param.epsilon=0.5; % papramater that specifies the width of the DiracDelta function
param.sigma=1.5;     % scale parameter in Gaussian kernel
se = strel('disk',2);        

minSizeNucleus = 100;
maxSizeNucleus = 1000;

%% pre-processing
[BW2] = preProcess(im,0);

% % remove spurious regions
th = 700;
CC = bwconncomp(BW2);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = find(numPixels<th);
BW = BW2;
for v1 = 1:size(idx,2)
    BW(CC.PixelIdxList{idx(v1)}) = 0;
end
% figure,imshow(BW);title('BW');

% % update CC
CC = bwconncomp(BW);

% % using median
nucleusSize = [5 5];
medianImg = medfilt2(im,nucleusSize);

%% nucleus segmentation: otsu thesholding in each connected component region
for v1 = 1:CC.NumObjects
    connectCompBinary = (zeros(size(BW)));
    connectCompBinary(CC.PixelIdxList{v1}) = 1;
    connectComp = uint8(zeros(size(BW)));
    connectComp(CC.PixelIdxList{v1}) = medianImg(CC.PixelIdxList{v1});
    cellClumpOrig(:,:,v1) = connectComp;
    temp3 = medianImg(CC.PixelIdxList{v1});
    threshold = myotsuSTw1(double(temp3(:)),0.05);  % only innovation  
    cellClumpNucleusBinary(:,:,v1) = im2bw(connectComp,max(threshold));
    cellClumpNucleusBinary(:,:,v1) = addBorder(cellClumpNucleusBinary(:,:,v1));
    %% if cell clump is crowded, then do these steps
    tempvar0 = cellClumpNucleusBinary(:,:,v1);
    tempvar1 =  imfill(tempvar0,'holes');
    rrr = xor(tempvar0,tempvar1);
    tempvar2 = bwareaopen(rrr, minSizeNucleus);
    tempvar2CC = bwconncomp(tempvar2);
    numPixelsCC1 = cellfun(@numel,tempvar2CC.PixelIdxList);
    [~,idx1] = find(numPixelsCC1>maxSizeNucleus);
    for var1 = 1:length(idx1)
        co = medianImg(tempvar2CC.PixelIdxList{idx1(var1)});
        d = myotsuSTw1(double(co),0.05);
        newTemp = uint8(zeros(size(BW)));
        newTemp(tempvar2CC.PixelIdxList{idx1(var1)}) = medianImg(tempvar2CC.PixelIdxList{idx1(var1)});
        tt = im2bw(newTemp,d);
        
        cellClumpNucleusBinary(:,:,v1) = tt + cellClumpNucleusBinary(:,:,v1);
       
        %{
        threshold = myotsuSTw1(double(temp3(:)),0.015);
        cellClumpNucleusBinary(:,:,v1) = im2bw(connectComp,max(threshold));
        %}
    end
    %% remove small nucleus of less size
    tempvar0 = cellClumpNucleusBinary(:,:,v1);
    tempvar1 = 1 - tempvar0;
    cellClumpNucleusBinary(:,:,v1) = 1-bwareaopen(tempvar1, minSizeNucleus);    
% % %     figure, imshowpair(cellClumpOrig(:,:,v1),cellClumpNucleusBinary(:,:,v1),'montage');
end

% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % % % % % % display what has obtained
% % % % for v0=1:40
% % % % figure,subplot(121),imshow(cellClumpOrig(:,:,v0))
% % % % subplot(122),imshow(cellClumpNucleusBinary(:,:,v0))
% % % % pause
% % % % end
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % testing hole filling and line join removal
% % % % for vp = 1:CC.NumObjects
% % % %     tem = cellClumpNucleusBinary(:,:,vp);
% % % %     holeFilled = imfill(tem,'holes');
% % % % %     figure,imshow(holeFilled)
% % % %     
% % % %     seCutCytoplasm = strel('disk',4); 
% % % %     BWn12 = imerode(holeFilled,seCutCytoplasm);
% % % %     BWn12 = imdilate(BWn12,seCutCytoplasm);
% % % %     figure,imshow(BWn12);title('BWn12');
% % % % end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cytoplasm segmentation
count = 0;      % cytoplasm counter initialized to '0' for each image
for v1 = 1:CC.NumObjects
    temp1 = cellClumpOrig(:,:,v1);
    [croppedRegionOrig, ~] = cropRegion(im, CC.PixelIdxList{v1}, 10);
    temp2 = cellClumpNucleusBinary(:,:,v1);
    [croppedRegionBinary, beginPos] = cropRegion(temp2, CC.PixelIdxList{v1}, 10);
    t1 = 1-croppedRegionBinary;
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %     figure,imshow(croppedRegionOrig)
% % % %     figure,imshow(croppedRegionBinary)
% % % %     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    CC1 = bwconncomp(t1);

    % cut the cytoplasms of two cells joined by a  thin line
    holeFilled = imfill(1-t1,'holes');
% % % %     seCutCytoplasm = strel('disk',4); 
% % % %     BWn12 = imerode(holeFilled,seCutCytoplasm);
% % % %     BWn12 = imdilate(BWn12,seCutCytoplasm);
% % % %     figure,imshow(BWn12);title('BWn12');
% % % %     BWn2 = 1- (BWn12 & (1-BWn1));
% % % %     figure,imshow(BWn2);title('BWn12');
    fillNucleus = imfill(temp2,'holes');
    if CC1.NumObjects <= 2 % only one nucleus is present
        count = count + 1;        
        cellSegmentationResult{count,1} = fillNucleus;
    else        
        numNucleus = CC1.NumObjects - 1;  % white region contains "nucleus + non-cytoplasm region"
        numPixels = cellfun(@numel,CC1.PixelIdxList);
        [~,idx] = sort(numPixels,'descend');        
        
        param.type = 'cytoplasm';
        count1 = 0;
        for v2 = 2:length(idx) % the largest connected component being the background is not considered
            [I2,J2] = ind2sub(size(croppedRegionBinary),CC1.PixelIdxList{v2});
            param.pos = [mean(I2), mean(J2)];
            param.rad = sqrt(length(I2))*1.5;
            res(:,:,v2-1) = cytoplasmSegmentLSM(croppedRegionOrig,holeFilled, param);            
            
            segmentLSM = res(:,:,v2-1)>0;
            erImg = imerode(segmentLSM,se);
            diaImg = imdilate(erImg,se);
% % % % %             figure(10), subplot(121);imshow(croppedRegionOrig);subplot(122);imshow(diaImg);
% % % %             C = imfuse(croppedRegionOrig,diaImg*155,'blend','Scaling','joint');
% % % %             figure, imshow(C); %pause
        
            z = zeros(size(im));
            z(beginPos(1):beginPos(1)+size(diaImg,1)-1,beginPos(2):beginPos(2)+size(diaImg,2)-1) = diaImg;
            
            minSizeCytoplasm = 700;
            if sum(z(:))>minSizeCytoplasm % check if the size of region comply with the size of cytoplasm
                count1 = count1 + 1;
                currentClump{count1,1} = z;
            end            
% % % %             figure,imshow(SegmentationResult{imgNo,1}{count,1});
        end
        
        %% check if all region of cell clump is segmented properly...
        if exist('currentClump')
            for i=1:size(currentClump,1)
                count = count + 1;
                cellSegmentationResult{count,1} = currentClump{i,1};                
            end
         end
        clear res;
        clear currentClump;
        
    end
    
end
end