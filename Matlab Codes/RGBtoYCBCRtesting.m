Image = imread('me.jpg');
I = rgb2ycbcr(Image);
%imshow(I)
for i = 1:800
    for j = 1:600
        if(I(i,j,2)>80 && I(i,j,2)<135)
            if(I(i,j,3)>131 && I(i,j,3)<185)
                I(i,j,1) = 0;
                I(i,j,2) = 0;
                I(i,j,3) = 0;
            end
        end
    end
end
imshow(I)