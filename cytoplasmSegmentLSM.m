function res = cytoplasmSegmentLSM(Img, holeFilled, param)

% Img=double(Img);
Img=double(Img.*uint8(holeFilled));

%% parameter setting
seCutCytoplasm = strel('disk',15);
holeFilled = imerode(holeFilled,seCutCytoplasm); % truncate the boundary regions

if strcmp('cytoplasm',param.type)
%     disp('cytoplasm segmentation....')
    timestep=param.timestep;  % time step
    mu=param.mu;  % coefficient of the distance regularization term R(phi)
    iter_inner=param.iter_inner;
    iter_outer=param.iter_outer;
    lambda=param.lambda; % coefficient of the weighted length term L(phi)
    alfa=param.alfa;  % coefficient of the weighted area term A(phi)
    epsilon=param.epsilon; % papramater that specifies the width of the DiracDelta function
    sigma=param.sigma;     % scale parameter in Gaussian kernel
    initialLSF = double(drawDisc(size(Img), param.pos, param.rad).*holeFilled)*8 - 2;
%     figure(12), imagesc(initialLSF)
elseif strcmp('nucleus',param.type)
    disp('nucleus segmentation....')
    timestep=param.timestep;  % time step
    mu=param.mu;  % coefficient of the distance regularization term R(phi)
    iter_inner=param.iter_inner;
    iter_outer=param.iter_outer;
    lambda=param.lambda; % coefficient of the weighted length term L(phi)
    alfa=param.alfa;  % coefficient of the weighted area term A(phi)
    epsilon=param.epsilon; % papramater that specifies the width of the DiracDelta function
    sigma=param.sigma;     % scale parameter in Gaussian kernel
    initialLSF = double(param.mask);
else
    timestep=5;  % time step
    mu=.25/timestep;  % coefficient of the distance regularization term R(phi)
    iter_inner=10;
    iter_outer=50;
    lambda=5; % coefficient of the weighted length term L(phi)
    alfa=1.5;  % coefficient of the weighted area term A(phi)
    epsilon=1.5; % papramater that specifies the width of the DiracDelta function
    sigma=1.5;     % scale parameter in Gaussian kernel
    initialLSF = double(drawDisc(size(Img), param.pos, param.rad))*8 - 2;
end

G=fspecial('gaussian',5,sigma);
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

phi=initialLSF;
% figure(1);
% mesh(-phi);   % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% title('Initial level set function');
% view([-80 35]);


potential=2;  
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end


% start level set evolution
for n=1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
end

% refine the zero level contour by further level set evolution with alfa=0
alfa=0;
phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);
%{
figure, 
subplot(131),imshow(uint8(Img));
subplot(132), imshow(phi);
subplot(133), imshow(initialLSF);
%}
res = phi>0;
res = res.*Img;