Image = imread('skin1.jpg');
I = rgb2ycbcr(Image);
%I = im2double(Image);
%imshow(I)
b = 0;
r = 0;
sum = 0;
COL = 800;
ROW = 600;
for i = 1:COL
    for j = 1:ROW
        b = b + double(I(i,j,2));
        r = r + double(I(i,j,3));
        sum = sum + 1;
        
    end
end
b = b/sum;
r = r/sum;

db = 0;
dr = 0;
c11 = 0;
c12 = 0;
c21 = 0;
c22 = 0;
c = [0 0;0 0];
for i = 1:COL
    for j = 1:ROW
        db = double(I(i,j,2)) - b;
        dr = double(I(i,j,3)) - r;
        c11 = c11 + db*db;
        c12 = c12 + db*dr;
        c21 = c21 + dr*db;
        c22 = c22 + dr*dr;
    end
end
c11 = c11/sum;
c12 = c12/sum;
c21 = c21/sum;
c22 = c22/sum;
c = [c11 c12;c21 c22];
ic = inv(c);

Itest = imread('lab.jpg');
I2 = rgb2ycbcr(Itest);

M1 = [0;0];
M2 = [0 0];
Tb = 0;
Tr = 0;
p = 0;
for i = 1:COL
    for j = 1:ROW
        Tb = double(I2(i,j,2));
        Tr = double(I2(i,j,3));
        Tb = Tb - b;
        Tr = Tr - r;
        M1 = [Tb;Tr];
        M2 = [Tb Tr];
        p = exp((-0.5)*M2*ic*M1);
        if(p>0.4)
            I2(i,j,1) = 0;
            I2(i,j,2) = 0;
            I2(i,j,3) = 0;
        end
        
    end
end





imshow(I2)