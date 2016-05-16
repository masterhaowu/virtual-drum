clc;
close;
clear all;
Obj = VideoReader('sss.MOV');

Sframe = 27;
Eframe = Sframe + 15;

Analysis = read(Obj, [Sframe Eframe]);

%Motion = read(Obj, [Sframe-5 Eframe]);

Ans = zeros(5,4);
Frame_prev2 = zeros(1080,1920,3);
Frame_prev1 = zeros(1080,1920,3);
Frame_prev = zeros(1080,1920,3);
for I = 1:1
    Cframe = Sframe + I - 1;
    Iori = Analysis(:,:,:,I*3);
    %p1frame = Sframe - 1;
    Ip1 = read(Obj, [Sframe-1 Sframe-1]);
    
    Itest = Iori;
    I4 = rgb2ycbcr(Itest);
    
    ic = [0.06 0.0708; 0.0708 0.0995];
    b = 112.3835;
    r = 147.3064;
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
    
    IPoM = zeros(COL,ROW);
    for i = 1:COL
        for j = 1:ROW
            IPoM(i,j) = (abs(double(Ip1(i,j,1)) - double(Iori(i,j,1))) + abs(double(Ip1(i,j,2)) - double(Iori(i,j,2))) + abs(double(Ip1(i,j,3)) - double(Iori(i,j,3))))/300;
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
            IPoA(i,j) = IPoB(i,j)*IPoM(i,j);

        end
    end
IPoA2 = IPoA;
    for i=1:COL
        for j=1:ROW
            if(IPoA(i,j)<0.01)
                IPoA2(i,j) = 0;
            end
        end
    end
    
 
    for i=1:COL
        for j=1:ROW
            I3(i,j) = IPoA2(i,j)*IPoE(i,j);
            if(I3(i,j)>0.005)
                I6(i,j) = 0;
            else
                I6(i,j) = 255;
            end
        end
    end
     imshow(I6)
     figure
     imshow(IPoA2)
%      figure
%      imshow(IPoB)
%      figure
%      imshow(IPoE)
%      figure
%      imshow(IPoM)


    
    
    
end