close all;
clear;
clc;

addpath('./DRLSE_v0/DRLSE_v0')
addpath('../ISBI 2015 challenge/EvaluationCode/EvaluationCode')

%% load cytoplasm ground truth

f1 = '../ISBI 2015 challenge/Training_R1_01Dec2014/Training';
f2 = '../ISBI 2015 challenge/Training_R2_Jan2015';

folders = {'seg_frame004_png','seg_frame011_png','seg_frame013_png',...
    'seg_frame014_png','seg_frame007_png','seg_frame010_png',...
    'seg_frame016_png','seg_frame017_png'};
GT_Stack = {};

for i=1:8
    if i<=4
        GT_Stack{i,1}.path = strcat(f1,'/',folders{1,i});
    else
        GT_Stack{i,1}.path = strcat(f2,'/',folders{1,i});
    end
end

for i=1:8
CytoGroundTruth{i,1} = getCytoplasmGT(GT_Stack{i,1}.path);
end

save CytoGroundTruth.mat CytoGroundTruth -v7.3


%% load image and segment cytoplasm

% get the paths to the image stacks

folders = {'frame004_stack','frame011_stack','frame013_stack','frame014_stack','frame007_stack','frame010_stack','frame016_stack','frame017_stack'};
imgStack = {};

for i=1:8 % for each folder(8 folders)
    if i<=4
        imgStack{i,1}.path = strcat(f1,'/',folders{1,i});
    else
        imgStack{i,1}.path = strcat(f2,'/',folders{1,i});
    end
end

tic
for imgNo=1:8
    % load image
    [~, images] = ReadImgs(imgStack{imgNo,1}.path,'*.png'); %images wil contain the stack of all images from differennt focal plane
    opt.WSize = 69; % options/ parameters passed to fstack_mod
    opt.Alpha=0.2;
    opt.Sth=13;
    im = fstack_mod(images,opt);        % gives you the "merged version" of the EDF image
    
    cellSegmentationResult = cellSegmentation(im);  %perform segmentaion
    
    SegmentationResult{imgNo,1} = cellSegmentationResult;
end
toc
save SegmentationResult.mat SegmentationResult -v7.3

%
%% evaluate the results
load('SegmentationResult.mat');
load('CytoGroundTruth.mat');

[meanDice70, ...
meanFNR70_object, ...
meanTPR70_pixel, ...
meanFPR70_pixel, ...
stdDice70, ...
stdTPR70_pixel, ...
stdFPR70_pixel, ...
stdFNo_70] ...
    = evaluateCytoSegmentation(CytoGroundTruth, SegmentationResult);

%}