function res = pad_array(img,pad)
% matlab inbuilt function padarray can do this and more
res = zeros(size(img)+2*pad);
res(pad(1)+1:end-pad(1),pad(2)+1:end-pad(1)) = img;
end