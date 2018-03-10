function T = myotsuSTw1(I,w1)
% Author: S L Happy
% Input:
%     I - array/matrix in double datatype
%     w1 - propability of class 1 (dark region in segmentation)
% Use:
%     y = imread('cameraman.tif');
%     th = myotsuSTw1(double(y),0.07)*256;
% create histogram
nbins = 256; 
[x,~] = hist(I(:),1:nbins);
% calculate probabilities
p = x./sum(x);
% initialisation
prob1 = zeros(1,nbins); 
mu1 = zeros(1,nbins); 
mu2 = zeros(1,nbins);

for t = 1:nbins,
    prob1(t) = sum(p(1:t));
end
if w1==1
    p1 = p;
else
    th_temp = find(prob1>w1);
    th = th_temp(1);
%     th = round(w1*nbins);
    x1 = [p(1:th)*(1-prob1(th)) p(th+1:end)*prob1(th)];
    p1 = x1./sum(x1);
end
for t = 1:nbins,
    prob1(t) = sum(p1(1:t));
    mu1(t) = sum(p1(1:t).*[1:t])/prob1(t);
    mu2(t) = sum(p1(t+1:nbins).*[t+1:nbins])/(1-prob1(t));
end
sigma = ((mu1-mu2).^2).*(prob1.*(1-prob1));
idx = find(sigma == max(sigma));
if ~isempty(idx)
    T = idx(1)/nbins;
else
    T = 0;
end
% figure, bar(sigma)