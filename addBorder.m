function img = addBorder(img)
% Author: S L Happy

if sum(img(1,2:size(img,2)-1)')>0
    img(1,:)=1;
end
if sum(img(size(img,1),2:size(img,2)-1)')>0
    img(size(img,1),:)=1;
end
if sum(img(2:size(img,1)-1,1))>0
    img(:,1)=1;
end
if sum(img(2:size(img,1)-1,size(img,2)))>0
    img(:,size(img,2))=1;
end

end