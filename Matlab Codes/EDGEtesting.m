Image = imread('lab2.jpg');
I = rgb2gray(Image);
BW1 = edge(I,'sobel');
BW2 = edge(I,'canny');
BW3 = edge(I,'log');
figure;
imshow(BW2)
