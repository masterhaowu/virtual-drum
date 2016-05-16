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

Itest = imread('lol3.jpg');
I4 = rgb2ycbcr(Itest);
I2 = I4;
COL = 800;
ROW = 1060;
M1 = [0;0];
M2 = [0 0];
Tb = 0;
Tr = 0;
p = 0;
for i = 2:COL-1
    for j = 2:ROW-1
        I2(i,j,1) = (double(I4(i,j,1)) + double(I4(i-1,j,1)) + double(I4(i+1,j,1)) + double(I4(i,j+1,1)) + double(I4(i,j-1,1)))/5;
        I2(i,j,2) = (double(I4(i,j,2)) + double(I4(i-1,j,2)) + double(I4(i+1,j,2)) + double(I4(i,j+1,2)) + double(I4(i,j-1,2)))/5;
        I2(i,j,3) = (double(I4(i,j,3)) + double(I4(i-1,j,3)) + double(I4(i+1,j,3)) + double(I4(i,j+1,3)) + double(I4(i,j-1,3)))/5;
        
    end
end
for i = 1:COL
    for j = 1:ROW
        Tb = double(I2(i,j,2));
        Tr = double(I2(i,j,3));
        Tb = Tb - b;
        Tr = Tr - r;
        M1 = [Tb;Tr];
        M2 = [Tb Tr];
        p = exp((-0.5)*M2*ic*M1);
        if(p>0.3)
            Itest(i,j,1) = 0;
            Itest(i,j,2) = 0;
            Itest(i,j,3) = 0;
        else
            Itest(i,j,1) = 255;
            Itest(i,j,2) = 255;
            Itest(i,j,3) = 255;
        end
        
    end
end

I3 = rgb2gray(Itest);

DIV = 10;
COLS = COL/DIV;
ROWS = ROW/DIV;
SDIV = zeros(COLS,ROWS);
temp = 0;
MAX = DIV*DIV*255;
valc = 0;
valr = 0;
for i = 1:COLS
    for j = 1:ROWS
        
        
        for ii = 1:DIV
            for jj = 1:DIV
                valc  = (i-1)*DIV + ii;
                valr = (j-1)*DIV + jj;
                temp = temp + double(I3(valc,valr));
            end
        end
        SDIV(i,j) = temp/MAX;
        temp = 0;
        
        
    end
end

for i = 2:COLS-1
    for j = 2:ROWS-1
        blackcheck = 0;
        if(SDIV(i,j)<0.7)
            blackcheck = 1;
        else
            if(SDIV(i-1,j)<0.7)
                blackcheck = 1;
            end
            if(SDIV(i+1,j)<0.7)
                blackcheck = 1;
            end
            if(SDIV(i,j-1)<0.7)
                blackcheck = 1;
            end
            if(SDIV(i,j+1)<0.7)
                blackcheck = 1;
            end
        end
    if(blackcheck == 0)
        for ii = 1:DIV
            for jj = 1:DIV
                valc  = (i-1)*DIV + ii;
                valr = (j-1)*DIV + jj;
                I3(valc,valr) = 255;
            end
        end
    end
        
    end
end



imshow(I3)