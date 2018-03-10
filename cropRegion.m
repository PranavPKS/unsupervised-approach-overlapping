function [res, beginPos] = cropRegion(im, pixelIndexes, th)
% Usage:
%     [croppedImage, beginPos] = cropRegion(image, CC.PixelIdxList{index}, border_width);
%     e.g. : [croppedRegion, beginPos] = cropRegion(im, CC.PixelIdxList{3}, 10);

    [I1,J1] = ind2sub(size(im),pixelIndexes);
    r1 = min(I1) - th;
    r2 = max(I1) + th;
    c1 = min(J1) - th;
    c2 = max(J1) + th;
    if r1<1, r1 = 1; end
    if c1<1, c1 = 1; end
    if r2>size(im,1), r2 = size(im,1); end
    if c2>size(im,2), c2 = size(im,2); end
    res = im(r1:r2,c1:c2);    
    beginPos = [r1,c1];
end