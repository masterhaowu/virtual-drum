clc;
close;
clear all;
Obj = VideoReader('sss.MOV');

Analysis = read(Obj, [27 41]);

Ans = zeros(5,4);
Frame_prev2 = zeros(1080,1920,3);
Frame_prev1 = zeros(1080,1920,3);
Frame_prev = zeros(1080,1920,3);




for I = 1:5
    Ioriginal = Analysis(:,:,:,I*3);
    Itest = Ioriginal;
    I4 = rgb2ycbcr(Itest);
    
    


ic = [0.06 0.0708; 0.0708 0.0995];
b = 112.3835;
r = 147.3064;

%End of constant definition

% Itest = imread('lab.jpg');
% Itest = Analysis(:,:,:,1);
% I2 = rgb2ycbcr(Itest);

IMAGE_SIZE = size(I4);
COL = IMAGE_SIZE(1);
ROW = IMAGE_SIZE(2);
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



%here i divide the image to smaller pieces and to completely remove the piece that has less
%than certain amount of black pixels while its neighbor also does not have
%many black pixels. (This is to make sure i dont accidently remove the part that is next to the hand but ends up in a region that are mostly white)

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


HALFROW = ROW/2;
%TMiddle = 5*COL; 
Tall = 0; %Tall is all black pixels in the middle rows
Ti = 0; %Ti is black pixels in the middle cols from col 0 to col i
Stop = 0;
CutCOL = 0;
for i=1:COL
    if(I3(i,HALFROW-2) == 0)
        Tall = Tall + 1;
    end
    if(I3(i,HALFROW-1) == 0)
        Tall = Tall + 1;
    end
    if(I3(i,HALFROW) == 0)
        Tall = Tall + 1;
    end
    if(I3(i,HALFROW+1) == 0)
        Tall = Tall + 1;
    end
    if(I3(i,HALFROW+2) == 0)
        Tall = Tall + 1;
    end
end
for i=1:COL
    if(I3(i,HALFROW-2) == 0)
        Ti = Ti + 1;
    end
    if(I3(i,HALFROW-1) == 0)
        Ti = Ti + 1;
    end
    if(I3(i,HALFROW) == 0)
        Ti = Ti + 1;
    end
    if(I3(i,HALFROW+1) == 0)
        Ti = Ti + 1;
    end
    if(I3(i,HALFROW+2) == 0)
        Ti = Ti + 1;
    end
    Pmiddle = Ti/Tall;
    if((Stop == 0)&&(Pmiddle>0.95))
        CutCOL = i;
        Stop = 1;
    end
end

%if the percentage of the total black pixels in the middle rows that we selected are low, we
%assume the face is not in the shot and thus do not attempt to remove the
%face
Pmiddle = Ti/(5*CutCOL);
if(Pmiddle<0.1)
    CutCOL = 1;
end

imshow(flipud(I3))
figure
%Transfer of data
Frame_prev2 = Frame_prev1;
Frame_prev1 = Frame_prev;
Frame_prev = double(Ioriginal);

%Find location of black pixels
Black_pix = zeros(2, 1080*1920);
c = 1;
for i = 1:COL
    for j = 1:ROW
        if I3(i,j) == 0
            Black_pix(:,c) = [i, j];
            c = c+1;
        end
    end
end

%motion detection
c = 1;
if I > 2
    while Black_pix(1,c) ~= 0
        i = Black_pix(1,c);
        j = Black_pix(2,c);
        if Frame_prev(i,j,1) > Frame_prev1(i,j,1) - 10 && Frame_prev(i,j,1) <...
                Frame_prev1(i,j,1) + 10 ...
            && Frame_prev(i,j,2) > Frame_prev1(i,j,2) - 10 ...
                && Frame_prev(i,j,2) < Frame_prev1(i,j,2) + 10 && Frame_prev(i,j,3) ...
                > Frame_prev1(i,j,3) - 10 && Frame_prev(i,j,3) < Frame_prev1(i,j,3)...
                + 10 ...
            && Frame_prev(i,j,1) > Frame_prev2(i,j,1) - 10 ...
                && Frame_prev(i,j,1) < Frame_prev2(i,j,1) + 10 && Frame_prev(i,j,2) ...
                > Frame_prev2(i,j,2) - 10 ...
                && Frame_prev(i,j,2) < Frame_prev2(i,j,2) + 10 && Frame_prev(i,j,3) ...
                > Frame_prev2(i,j,3) - 10 && Frame_prev(i,j,3) < Frame_prev2(i,j,3) + 10
            I3(i,j) = 255;
            fprintf('was this gone through');
        end
        c = c +1;
    end
end

%See effect of motion detection
imshow(flipud(I3))
figure

%now I am going to remove the head.

LLL = HALFROW; %the left point of the head
RRR = HALFROW; %the right point of the head
whitecounter = 0;
whitebar = floor(ROW/10); %this is defined as how many white pixels do we considered as white bar 
%whitebarT = whitebar*255;
for i = 1:CutCOL
    Rc = 0; %Right white counter
    Lc = 0; %Left white counter
    StopR = 0; %Stop the loop once we find the RRR
    StopL = 0; %Stop the loop once we find the LLL
    for j = 1:HALFROW - whitebar - 1
        Rbar = HALFROW + j;
        Lbar = HALFROW - j;
        if(I3(i,Rbar)==255 && StopR==0)
            Rc = Rc+ 1;
        end
        if(I3(i,Rbar)==0)
            Rc = 0;
        end
        if(Rc>whitebar && StopR==0)
            RRR = Rbar;
            StopR = 1;
        end
        
        if(I3(i,Rbar)==255 && StopL==0)
            Lc = Lc+ 1;
        end
        if(I3(i,Lbar)==0)
            Lc = 0;
        end
        if(Lc>whitebar && StopL==0)
            LLL = Lbar;
            StopL = 1;
        end
        
        
    end
    %RRR
    %now we have our RRR and LLL for col i. We paint everything between RRR
    %and LLL to white
    for k = LLL:RRR
        I3(i,k) = 255;
    end
end

end