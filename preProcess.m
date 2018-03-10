function [BW2] = preProcess(im,showImgs)
if(showImgs),figure,imshow(im);title('original');end
% using  mean, median
nucleusSize = [5 5];
im = medfilt2(im,nucleusSize);if(showImgs),figure,imshow(im);title('median filter');end

%% pre-processing
% histogram equalization
I_eq = adapthisteq(im);if(showImgs),figure,imshow(I_eq);title('hist_eq');end

% % % % image bottom hat
% % se = strel('disk',3);
% % I_eq = imsubtract(imadd(I_eq,imtophat(I_eq,se)), imbothat(I_eq,se));

% image h-maxima transform
% I1 = imhmin(im,60);if(showImgs),figure,imshow(I1);title('hmin');end

% image h-maxima transform
I2 = imhmax(I_eq,40);if(showImgs),figure,imshow(I2);title('hmax');end


%% Scene segmentation

% % % % Otsu thresholding
% % bw = im2bw(I_eq, graythresh(I_eq));
% % figure,imshow(bw)

regmax = imcomplement(imregionalmax(I2));
regmax = imfill(regmax, 'holes');
BW2 = bwareaopen(regmax, 40);
if(showImgs),figure,imshow(BW2);end

end