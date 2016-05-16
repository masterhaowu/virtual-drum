Obj = VideoReader('sex.MOV');
%Itest = imread('lol.jpg');
I4= read(Obj, [27 27]);
%I4 = rgb2ycbcr(Itest);
COL = 1080;
ROW = 1920;
M1 = [0;0];
M2 = [0 0];
Tb = 0;
Tr = 0;
p = 0;
I2 = I4;


%I put a filter here so the original image will be smoother

for i = 2:COL-1
    for j = 2:ROW-1
        I2(i,j,1) = (double(I4(i,j,1)) + double(I4(i-1,j,1)) + double(I4(i+1,j,1)) + double(I4(i,j+1,1)) + double(I4(i,j-1,1)) + double(I4(i-1,j+1,1)) + double(I4(i+1,j+1,1)) + double(I4(i+1,j-1,1)) + double(I4(i-1,j-1,1)))/9;
        I2(i,j,2) = (double(I4(i,j,2)) + double(I4(i-1,j,2)) + double(I4(i+1,j,2)) + double(I4(i,j+1,2)) + double(I4(i,j-1,2)) + double(I4(i-1,j+1,2)) + double(I4(i+1,j+1,2)) + double(I4(i+1,j-1,2)) + double(I4(i-1,j-1,2)))/9;
        I2(i,j,3) = (double(I4(i,j,3)) + double(I4(i-1,j,3)) + double(I4(i+1,j,3)) + double(I4(i,j+1,3)) + double(I4(i,j-1,3)) + double(I4(i-1,j+1,3)) + double(I4(i+1,j+1,3)) + double(I4(i+1,j-1,3)) + double(I4(i-1,j-1,3)))/9;
        
    end
end
imshow(I4)
figure
imshow(I2)
