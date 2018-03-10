function Im = fstack_mod(Scene, Options)
%fstack
%Focus stacking.
%
%SINTAX:
%     [~, images] = ReadImgs('E:\MIP\Assignments\Project\ISBI 2015 challenge\Testing_R2_Jan2015\frame019_stack\','*.png');
%     opt.WSize = 69;
%     opt.Alpha=0.2;
%     opt.Sth=13;
%     Im = fstack_mod(images,opt);
%
%DESCRIPTION:
%Generate all-in-focus (AIF) image from focus sequence.
%
%OUTPUTS:
% Im,       is the AIF image
%INPUTS:
% images,   is a cell array where each cell is an image 
%           (or a string with the path of an image).
% Options and their values (default in perenthesis)
% may be any of the following:
%   'FSize',      Size of focus measure window (9).
%   'Focus',      Vector with the focus of each frame.
%   'Alpha',      Parameter of the AIF algorithm (0.2).
%   'Sth',        Parameter of the AIF algorithm (13).%   
%
%For further details, see:
%Pertuz et. al. "Generation of all-in-focus images by
%noise-robust selective fusion of limited depth-of-field
%images" IEEE Trans. Image Process (in print).
%
%Said Pertuz
%Sep14/2011



%Parse inputs:


if ~isfield(Options,'Size')
    Im1 = (imread(Scene{1}));    
    Options.Size = [size(Im1) size(Scene,2)];
end
    M = Options.Size(1);
    N = Options.Size(2);
    P = Options.Size(3);
if ~isfield(Options,'Focus')
    Options.Focus = 1:Options.Size(3);
end


%********* Read images and compute fmeasure **********
%Initialize:
FM = zeros(Options.Size);

ImagesG = zeros(Options.Size);

%Read:
% fprintf('Fmeasure     ')
for p = 1:P
    Im = (imread(Scene{p}));
    ImagesG(:,:,p) = Im;
    FM(:,:,p) = gfocus(im2double(Im), Options.WSize);
%     fprintf('\b\b\b\b\b[%2.0i%%]',round(100*p/P))
end

%********** Compute Smeasure ******************
% fprintf('\nSMeasure     ')
[u s A Fmax] = gauss3P(Options.Focus, FM);
%Aprox. RMS of error signal as sum|Signal-Noise|
%instead of sqrt(sum(Signal-noise)^2):
Err = zeros(M,N);
for p = 1:P
    Err = Err + abs( FM(:,:,p) - ...
        A.*exp(-(Options.Focus(p)-u).^2./(2*s.^2)));
    FM(:,:,p) = FM(:,:,p)./Fmax;
%     fprintf('\b\b\b\b\b[%2.0i%%]',round(100*p/P))
end
H = fspecial('average', Options.WSize);
inv_psnr = imfilter(Err./(P*Fmax), H, 'replicate');

S = 20*log10(1./inv_psnr);
% S(isnan(S))=min(S(:));
% fprintf('\nWeights      ')
Phi = 0.5*(1+tanh(Options.Alpha*(S-Options.Sth)))/...
   Options.Alpha;
% Phi = 0.5*(1+tanh(Options.Alpha*(S-Options.Sth)));
% Phi(isnan(Phi)) = min(Phi(:));
% Phi = medfilt2(Phi, [3 3]);
%********** Compute weights: ********************
fun = @(phi,fm) 0.5+0.5*tanh(phi.*(fm-1));
for p = 1:P    
    FM(:,:,p) = feval(fun, Phi, FM(:,:,p));
%     fprintf('\b\b\b\b\b[%2.0i%%]',round(100*p/P))
end

%********* Fuse images: *****************
% fprintf('\nFusion ')
FMn = sum(FM,3);
    Im = uint8(sum((ImagesG.*FM), 3)./FMn);
% fprintf('[100%%]\n')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u s A Ymax] = gauss3P(x, Y)
%Internal parameter:
STEP = 2;
%%%%%%%%%%%%%%%%%%%
[M,N,P] = size(Y);
[Ymax, I] = max(Y,[ ], 3);
[IN,IM] = meshgrid(1:N,1:M);
Ic = I(:);
Ic(Ic<=STEP)=STEP+1;
Ic(Ic>=P-STEP)=P-STEP;
Index1 = sub2ind([M,N,P], IM(:), IN(:), Ic-STEP);
Index2 = sub2ind([M,N,P], IM(:), IN(:), Ic);
Index3 = sub2ind([M,N,P], IM(:), IN(:), Ic+STEP);
Index1(I(:)<=STEP) = Index3(I(:)<=STEP);
Index3(I(:)>=STEP) = Index1(I(:)>=STEP);
x1 = reshape(x(Ic(:)-STEP),M,N);
x2 = reshape(x(Ic(:)),M,N);
x3 = reshape(x(Ic(:)+STEP),M,N);
y1 = reshape(log(Y(Index1)),M,N);
y2 = reshape(log(Y(Index2)),M,N);
y3 = reshape(log(Y(Index3)),M,N);
c = ( (y1-y2).*(x2-x3)-(y2-y3).*(x1-x2) )./...
    ( (x1.^2-x2.^2).*(x2-x3)-(x2.^2-x3.^2).*(x1-x2) );
b = ( (y2-y3)-c.*(x2-x3).*(x2+x3) )./(x2-x3);
s = sqrt(-1./(2*c));
u = b.*s.^2;
a = y1 - b.*x1 - c.*x1.^2;
A = exp(a + u.^2./(2*s.^2));
% Fmax = x(I);
% IDX = imag(u)~=0|isnan(u);
% u(IDX) = Fmax(IDX);
% A(IDX) = 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FM = gfocus(Image, WSize)
MEANF = fspecial('average',[WSize WSize]);
U = imfilter(Image, MEANF, 'replicate');
FM = (Image-U).^2;
FM = imfilter(FM, MEANF, 'replicate');
end