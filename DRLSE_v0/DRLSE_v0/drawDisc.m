function res = drawDisc(sz, position, radius)

% Usage:
% centerX = 320;
% centerY = 240;
% position = [centerY, centerX];
% sz = [480,640];
% radius = 100;
% drawDisc(sz, position, radius)

imageSizeX = sz(2);
imageSizeY = sz(1);
[columnsInImage rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
res = (rowsInImage - position(1)).^2 ...
    + (columnsInImage - position(2)).^2 <= radius.^2;
% imshow(res) ;
end