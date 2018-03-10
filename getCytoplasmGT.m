function [cytogt] = getCytoplasmGT(path)
% get cytoplasm ground truth
%
% Input:
%   path - path to the training image folder containing the stack of images.
%
% Output:
%     cytogt - cytoplasm ground truth ( cell datatype)
% 
% Usage:
%   To display the annotation of image frame014:
%     path_1 = '..\ISBI 2015 challenge\Training_R1_01Dec2014\Training\seg_frame004_png';
%     [cytogt] = getCytoplasmGT(path_1);


    path = strcat(path,'/');
    d = dir(path);
    idCell = 1;
    for i = 1:size(d,1)
        if d(i,1).isdir == 0
            cytogt{idCell,1} = imread(strcat(path, d(i,1).name));
            idCell = idCell + 1;
        end
    end
end