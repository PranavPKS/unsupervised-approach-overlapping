%{
path = 'E:\face and emotion databases\Selected Cohn-canade images\only\disgust\';
data1 = ReadImgs(path,'*.png');
%}

function [X,images] = ReadImgs(Folder,ImgType)
% Author: S L Happy
    Folder = strcat(Folder,'/');
    Imgs = dir([Folder '/' ImgType]);
    x = pwd;
    cd (Folder)
    imInfo = imfinfo(Imgs(1).name);
    rows = imInfo.Width;
    cols = imInfo.Height;
    NumImgs = size(Imgs,1);
    X = zeros([rows cols NumImgs]);
    for i=1:NumImgs,
      imInfo = imfinfo(Imgs(i).name);
      image = imread(imInfo.Filename);
      images{i} = imInfo.Filename;
      if (size(image,3) == 1)
        X(:,:,i) = imresize(image, [rows cols]);
      else
        p = rgb2gray(image);        
        X(:,:,i) = imresize(p, [rows cols]);
      end
    end
    cd (x)
end