function [Imstar0] = starsimulator(alfa,delta,roll,magmin)
    % Writer: Nugraha Setya Ardi
    %
    % This function generates Imstar0 matrix taking alfa, delta, roll, and
    % threshold for star's magnitude.
    % Input:
    %   1. alfa: in deg
    %   2. delta: in deg
    %   3. roll: in deg
    %   4. magmin: threshold for star's magnitude. Stars dimmer than this
    %   value will be excluded
    % Output:
    %   Imstar0: w x l matrix of star image
    % Other parameters must be tuned manually.
    
    load('catalog.mat')

    alfa0 = alfa;
    delta0 = delta;
    phi0 = roll;

    MAfull = catalog(:,3);
    rax = catalog(:,2);
    dey = catalog(:,1);

    [m1,n1]=size(rax);

    %star sensor focal length
    f = 0.003044898

    %star sensor pixel
    l = 3280;
    w = 2464;

    %lenght per pixel
    myu = 1.12*10^(-6);

    %Star sensor FOV
    FOVy = 2*atan((myu*w/2)/f)*180/pi;
    FOVx = 2*atan((myu*l/2)/f)*180/pi;
    R = 0.5*sqrt(FOVx^2+FOVy^2)*pi/180;

    %starsensor ROI
    ROI = 50%round((5/1024)*w);

    %star sensor treshold
    treshold = 0;

    %CCD background intensity
    background = 0;

    %Image star treshold
    treshold = 1;

    %Rotation matrix from star sensor coor to celestial coor
    a0 = alfa0*pi/180;
    d0 = delta0*pi/180;
    ph0 = phi0*pi/180;

    M1 = [cos(a0-pi/2) -sin(a0-pi/2) 0;sin(a0-pi/2) cos(a0-pi/2) 0;0 0 1];
    M2 = [1 0 0;0 cos(d0+pi/2) -sin(d0+pi/2);0 sin(d0+pi/2) cos(d0+pi/2)];
    M3 = [cos(ph0) -sin(ph0) 0;sin(ph0) cos(ph0) 0;0 0 1];
    M = M1*M2*M3;

    %creating star image

    for i = 1:m1
        alfastart = a0-R/cos(d0);
        alfaend = a0+R/cos(d0);
        deltastart = d0-R;
        deltaend = d0+R;

        if rax(i) > alfastart && rax(i) < alfaend && dey(i) > deltastart && dey(i) < deltaend
            ra(i)=rax(i);
            de(i)=dey(i);
            ma(i)=MAfull(i);
        end
    end

    raFOV = ra(ra ~=0);
    deFOV = de(de ~=0);
    maFOV = ma(ma ~=0);
    [b1 b2] = size(raFOV);

    for i=1:b2
        a = raFOV(i);
        d = deFOV(i);
        xbar = cos(a)*cos(d);
        ybar = sin(a)*cos(d);
        zbar = sin(d);
        Xbar = [xbar;ybar;zbar];
        X = inv(M)*Xbar;
        x = X(1,1);
        y = X(2,1);
        z = X(3,1);
        x1(i)=f*x/z;
        y1(i)=f*y/z;
    end

    xtot = 2*tan((FOVx*pi/180)/2)*f;
    ytot = 2*tan((FOVy*pi/180)/2)*f;
    xpixel = l/xtot;
    ypixel = w/ytot;
    x1pixel = round(xpixel*x1);
    y1pixel = round(ypixel*y1);
    plot(x1pixel,y1pixel,'r*')

    Imstar0 = zeros(w,l);
    n = 1000;
    gmin = 0;
    gmax = 0;
    lmin = 1;
    lmax = l;
    wmin = 1;
    wmax = w;
    randomintensity = gmin + rand(1,n)*(gmax-gmin);
    lrand = round(lmin + rand(1,n)*(lmax-lmin));
    wrand = round(wmin + rand(1,n)*(wmax-wmin));
    for i=1:w
        for j=1:l
            Imstar0(i,j)=gmin + rand(1,1)*(gmax-gmin);
        end
    end

    [x1pixela x1pixelb] = size(x1pixel);

    for i=1:x1pixelb
        x1pix(i) = x1pixel(i)+l/2;
        y1pix(i) = -y1pixel(i)+w/2;
        x1disp0(i) = xpixel*x1(i)+l/2;
        y1disp0(i) = -ypixel*y1(i)+w/2;
        if x1pix(i) < 1 || x1pix(i) > l || y1pix(i) < 1 || y1pix(i) > w || x1disp0(i) < 1 || x1disp0(i) > l || y1disp0(i) < 1 || y1disp0(i) > w
            x1pix(i) = 0;
            y1pix(i) = 0;
            x1disp0(i) = 0;
            y1disp0(i) = 0;
            maFOV(i)=0;
        end
    end
    x1pixnew = x1pix(x1pix ~= 0 );
    y1pixnew = y1pix(y1pix ~= 0 );
    x1disp = x1disp0(x1disp0 ~= 0 );
    y1disp = y1disp0(y1disp0 ~=0 );
    maFOVt = maFOV(maFOV ~=0);
    [x1pixnewa x1pixnewb] = size(x1pixnew);
    matres = magmin;
    for i=1:x1pixnewb
        if maFOVt(i) > matres
            maFOVt0(i) = 0;
            x1pixnew0(i) = 0;
            y1pixnew0(i) = 0;
            x1disp0(i) = 0;
            y1disp0(i) = 0;
        else
            maFOVt0(i) = maFOVt(i);
            x1pixnew0(i) = x1pixnew(i);
            y1pixnew0(i) = y1pixnew(i);
            x1disp0(i) = x1disp(i);
            y1disp0(i) = y1disp(i);
        end
    end
    x1pixnew = x1pixnew0(x1pixnew0 ~=0);
    y1pixnew = y1pixnew0(y1pixnew0 ~=0);
    maFOVt = maFOVt0(maFOVt0 ~=0);
    x1disp = x1disp0(x1disp0 ~= 0 );
    y1disp = y1disp0(y1disp0 ~=0 );
    [x1pixnewa x1pixnewb] = size(x1pixnew);

    nmissmin = 30;
    nmissmax = 40;
    nmiss = round(nmissmin+rand(1,1)*(nmissmax-nmissmin));

    n1min = 1;
    n1max = x1pixnewb;
    numbmiss = round(n1min+rand(1,nmiss)*(n1max-n1min));

    [a b]=size(numbmiss);

    randmissstars = 0

    if randmissstars == 1


    for i=1:b
        x1pixnew(numbmiss(i))=0;
        y1pixnew(numbmiss(i))=0;
        x1disp(numbmiss(i))=0;
        y1disp(numbmiss(i))=0;
    end
    x1pixnew1 = x1pixnew(x1pixnew ~=0);
    y1pixnew1 = y1pixnew(y1pixnew ~=0);
    x1pixnew = x1pixnew1;
    y1pixnew = y1pixnew1;
    [x1pixnewa x1pixnewb] = size(x1pixnew);

    x1disp1 = x1disp(x1disp ~=0);
    y1disp1 = y1disp(y1disp ~=0);
    x1disp = x1disp1;
    y1disp = y1disp1;
    end
    for i=1:x1pixnewb
        rad = ROI;
        ystart = round(y1pixnew(i)-rad);
        yend = round(y1pixnew(i)+rad);

        if ystart <= 0
            ystart = 1;
        end

        if yend > w
            yend = w;
        end
            for u=ystart:yend
                xstart = round(x1pixnew(i)-sqrt(rad^2 - (u-y1pixnew(i))^2));
                xend = round(x1pixnew(i)+sqrt(rad^2 - (u-y1pixnew(i))^2));

                if xstart <= 0
                    xstart = 1;
                end

                if xend > l
                    xend = l;
                end

                for v=xstart:xend   
                    sigma = 1.5;
                    H = 1000*exp(-maFOVt(i)+1);
                    r = sqrt((x1pixnew(i)-v)^2 + (y1pixnew(i)-u)^2);
                    g1 = round(255*exp(-r));
                    g2 = H*exp(-((((u-y1disp(i))/w)*1024)^2+(((v-x1disp(i))/w)*1024)^2)/(2*sigma^2));
                    Imstar0(u,v)=Imstar0(u,v)+g2;
                end
            end
    end
end


