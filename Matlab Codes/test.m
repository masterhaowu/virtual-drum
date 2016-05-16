Image = imread('me.jpg');
I = rgb2ycbcr(Image);
for i = 1:960
    for j = 1:720
        if(I(i,j,2)>100 && I(i,j,2)<135)
            if(I(i,j,3)>141 && I(i,j,3)<165)
                I(i,j,1) = 0;
                I(i,j,2) = 0;
                I(i,j,3) = 0;
            end
        end
    end
end
imshow(I)