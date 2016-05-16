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

Obj = VideoReader('sex.MOV');
%Itest = imread('lol.jpg');
Itest = read(Obj, [27 27]);
Iori = Itest;
I4 = rgb2ycbcr(Itest);
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

IPoB = zeros(COL,ROW);

for i = 1:COL
    for j = 1:ROW
        Tb = double(I2(i,j,2));
        Tr = double(I2(i,j,3));
        Tb = Tb - b;
        Tr = Tr - r;
        M1 = [Tb;Tr];
        M2 = [Tb Tr];
        p = exp((-0.5)*M2*ic*M1);
        IPoB(i,j) = p;
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
imshow(I3)
figure


%here i divide the image to smaller pieces and to completely remove the piece that has less
%than certain amount of black pixels while its neighbor also does not have
%many black pixels. (This is to make sure i dont accidently remove the part that is next to the hand but ends up in a region that are mostly white)

DIV = 20;
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
                IPoB(valc,valr) = IPoB(valc,valr)/3;
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
for i=1:COL/2
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
for i=1:COL/2
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
        %I3(i,k) = 255;
        IPoB(i,k) = 0;
    end
end


 

%Hand Mean Program Starts
% Ans = zeros(1,4);

%Segmentation for right hand
Y = 0;
X = 0;
YS = 0;
XS = 0;
Number = 0;
%First Distance filter
for i = 1:COL
    for j = 1:HALFROW
        if I3(i,j,1) == 0
            X = X+ j;
            Y = Y+ i;
            XS = XS + j^2;
            YS = YS + i^2;
            Number = Number + 1;
            
        end
    end
end

MeanX = X/Number;
MeanY = Y/Number;
VarianceX = XS/Number - MeanX^2;
VarianceY = YS/Number - MeanY^2;

% for i = 1:COL             %rectangular mean/variance cutout
%     for j = 467:(467+ROW)
%         if j < (MeanX - 1.5*sqrt(VarianceX)) || j > (MeanX + 1.5*sqrt(VarianceX)) || i < (MeanY - sqrt(VarianceY)) || i > (MeanY + sqrt(VarianceY))
%         I4(i,j) = 255;
%         end
%     end
% end
 
for i = 1:COL                   %round distance cutout
    for j = 1:HALFROW
        if ((i-MeanY)^2/(2.5*VarianceX) + (j-MeanX)^2/(2.5*VarianceY)) > 1
            I3(i,j) = 255;
        end
    end
end

% imshow(I3)
% figure

%Obtaining the new mean
Y = 0;
X = 0;
YS = 0;
XS = 0;
Number = 0;

for i = 1:COL
    for j = 1:HALFROW
        if I3(i,j,1) == 0
            X = X+ j;
            Y = Y+ i;
            %XS = XS + j^2;
            %YS = YS + i^2;
            Number = Number + 1;
            
        end
    end
end
RIGHTMeanX = X/Number;
RIGHTMeanY = Y/Number;

Ans(I,1) = RIGHTMeanX;
Ans(I,2) = RIGHTMeanY;

% for i = 1:COL
%     for j = 1:ROW
%         if (sqrt((i-MeanY)^2 + (j-MeanX)^2) > 170)
%             I4(i,j) = 255;
%         end
%     end
% end


%Left hand segmentation
Y = 0;
X = 0;
YS = 0;
XS = 0;
Number = 0;


%First Distance filter
for i = 1:COL
    for j = (ROW-HALFROW):ROW
        if I3(i,j,1) == 0
            X = X+ j;
            Y = Y+ i;
            XS = XS + j^2;
            YS = YS + i^2;
            Number = Number + 1;
            
        end
    end
end

MeanX = X/Number;
MeanY = Y/Number;
VarianceX = XS/Number - MeanX^2;
VarianceY = YS/Number - MeanY^2;

% for i = 1:COL             %rectangular mean/variance cutout
%     for j = 467:(467+ROW)
%         if j < (MeanX - 1.5*sqrt(VarianceX)) || j > (MeanX + 1.5*sqrt(VarianceX)) || i < (MeanY - sqrt(VarianceY)) || i > (MeanY + sqrt(VarianceY))
%         I4(i,j) = 255;
%         end
%     end
% end

%I33 = I3;

for i = 1:COL                   %round distance cutout
    for j = (ROW-HALFROW):ROW
        if ((i-MeanY)^2/(2.5*VarianceX) + (j-MeanX)^2/(2.5*VarianceY)) > 1
            %I3(i,j) = 255;
            %Itest(i,j,1) = 255;
            %Itest(i,j,2) = 255;
            %Itest(i,j,3) = 255;
            %IPoB(i,j) = IPoB(i,j)/3;
        end
    end
end
Igray = rgb2gray(Iori);
BW1 = edge(Igray,'canny');
% for i = 1:COL
%     for j = 1:ROW
%         if(I3(i,j)==0 && BW1(i,j)==0)
%             I3(i,j) = 255;
%         end
%     end
% end


Iedge = I3;
IPoE = zeros(COL,ROW);
IPoA = zeros(COL,ROW);

for i = 2:COL-1
    for j = 2:ROW-1
        IPoE(i,j) = (double(BW1(i,j)) + double(BW1(i-1,j)) + double(BW1(i+1,j)) + double(BW1(i,j+1)) + double(BW1(i,j-1)) + double(BW1(i-1,j+1)) + double(BW1(i+1,j+1)) + double(BW1(i+1,j-1)) + double(BW1(i-1,j-1)))/9;
        IPoA(i,j) = IPoE(i,j)*IPoB(i,j);
        
    end
end

for i=1:COL
    for j=1:ROW
        if(IPoA(i,j)>0.015)
            I3(i,j) = 0;
        else
            I3(i,j) = 255;
        end
    end
end
 imshow(I3)
 figure
