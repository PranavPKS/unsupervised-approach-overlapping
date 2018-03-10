close all;
clear;
clc;

addpath('./DRLSE_v0/DRLSE_v0')
addpath('../ISBI 2015 challenge/EvaluationCode/EvaluationCode')

%% load image and segment cytoplasm

% get the paths to the image stacks
f1 = '../ISBI 2015 challenge/Testing_R1_01Dec2014/Testing';
f2 = '../ISBI 2015 challenge/Testing_R2_Jan2015';

folders = {'frame000_stack','frame001_stack','frame002_stack','frame015_stack',...
    'frame018_stack','frame019_stack','frame020_stack','frame021_stack','frame022_stack'};
imgStack = {};

for i=1:9
    if i<=4
        imgStack{i,1}.path = strcat(f1,'/',folders{1,i});
    else
        imgStack{i,1}.path = strcat(f2,'/',folders{1,i});
    end
end

tic
for imgNo=1:8
    % load image
    [~, images] = ReadImgs(imgStack{imgNo,1}.path,'*.png');
    opt.WSize = 69;
    opt.Alpha=0.2;
    opt.Sth=13;
    im = fstack_mod(images,opt);
    
    cellSegmentationResult = cellSegmentation(im);
    
    SegmentationResult{imgNo,1} = cellSegmentationResult;
end
toc
save SegmentationResultTestset.mat SegmentationResult -v7.3
