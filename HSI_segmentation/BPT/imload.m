function [im,imrgb,varargout] = imload(imname,varargin)

switch imname
    case 'beach'
        D = load('Data/beach');
        im = D.im;
        imrgb = im;
    case 'gtnp'
        D = load('Data/gtnp');
        im = D.im;
        imrgb = im;
    case 'highway'
        D = load('Data/highway');
        im = D.im;
        imrgb = im;
    case 'grenoble'
        D = load('Data/grenoble');
        im = D.im;
        imrgb = im;
    case 'monument'
        D = load('Data/monument');
        im = D.im;
        imrgb = im;
    case 'plane'
        D = load('Data/plane');
        im = D.im;
        imrgb = im;
    case 'stanford'
        D = load('Data/stanford');
        im = D.im;
        imrgb = im;
    case 'mountain'
        D = load('Data/mountain');
        im = D.im;
        imrgb = im;
    case 'elephant'
        D = load('Data/elephant');
        im = D.im;
        imrgb = im;
    case 'bus'
        D = load('Data/bus');
        im = D.im;
        imrgb = im;
    case 'food'
        D = load('Data/food');
        im = D.im;
        imrgb = im;
    case 'barn'
        D = load('Data/barn');
        im = D.im;
        imrgb = im;
    case 'toyimage'
        D = load('Data/toyimage');
        im = D.im;
        imrgb = im;
        if nargout==3
            lidar = D.lidar;
            varargout(1) = {lidar};
        end
    case 'cuprite'
        D = load('Data/cuprite');
        im = D.im;
        imrgb = D.imrgb;
%        im = cuprite(101:end,161:end,:);
%        v = [37,19,12];
    case 'pavia_univ'
        D = load('Data/pavia_univ');
        pavia = D.pavia;
        im = pavia;
        imrgb = D.imrgb;
%         v = [50 33 12];
    case 'moffett1'
        D = load('Data/moffett1');
        im = D.im;
        imrgb = D.imrgb;
    case 'moffett3'
        D = load('Data/moffett3');
        im = D.im;
        imrgb = D.imrgb;
    case 'washington_full'
        D = load('Data/washington_full');
        im = D.im;
        imrgb = D.imrgb;
    case 'washington_top'
        D = load('Data/washington_top');
        im = D.im;
        imrgb = D.imrgb;
    case 'washington_bottom'
        D = load('Data/washington_bottom');
        im = D.im;
        imrgb = D.imrgb;
    case 'salinas'
        D = load('Data/salinas');
        im = D.salinas;
        imrgb = D.imrgb;
%         v = [37 16 9];
    case 'forest'
        D = load('Data/forest');
        im = D.im;
        imrgb = rescale(im(:,:,[6 7 8]),1);
%         imrgb = D.imrgb;
    case 'nanawale1'
        D = load('Data/nanawale1');
        nanawale1 = D.nanawale1;
        im = nanawale1;
        v = [10 7 3];
        imrgb = rescale(im(:,:,v(1)),1);
        imrgb(:,:,2) = rescale(im(:,:,v(2)),1);
        imrgb(:,:,3) = rescale(im(:,:,v(3)),1);
        if nargout==3
            lidar = D.lidar1;
            varargout(1) = {lidar};
        end
    case 'nanawale2'
        D = load('Data/nanawale2');
        nanawale2 = D.nanawale2;
        im = nanawale2;
        v = [10 7 3];
        imrgb = rescale(im(:,:,v(1)),1);
        imrgb(:,:,2) = rescale(im(:,:,v(2)),1);
        imrgb(:,:,3) = rescale(im(:,:,v(3)),1);
        if nargout==3
            lidar = D.lidar2;
            varargout(1) = {lidar};
        end
    case 'nanawale_full'
        D = load('Data/nanawale_full');
        im = D.im;
        imrgb = D.imrgb;
        if nargout==3
            lidar = D.lidar;
            varargout(1) = {lidar};
        end
    case 'snowflake'
        D = load('Data/snowflake');
        snowflake = D.snowflake;
        im = snowflake;
        v = [1 2 3];
        imrgb = rescale(im(:,:,v(1)),1);
        imrgb(:,:,2) = rescale(im(:,:,v(2)),1);
        imrgb(:,:,3) = rescale(im(:,:,v(3)),1);
    case 'dfc2013_full'
        D = load('Data/dfc2013_full');
        im = D.dfc2013h;
        imrgb = D.imrgb;
%         v = [53 40 15];
        if nargout==3
            lidar = D.dfc2013l;
            varargout(1) = {lidar};
        end
    case 'dfc2013_left'
        D = load('Data/dfc2013_left');
        im = D.dfc2013h;
        imrgb = D.imrgb;
%         v = [53 40 15];
        if nargout==3
            lidar = D.dfc2013l;
            varargout(1) = {lidar};
        end
    case 'dfc2013_center'
        D = load('Data/dfc2013_center');
        im = D.dfc2013h;
        imrgb = D.imrgb;
%         v = [53 40 15];
        if nargout==3
            lidar = D.dfc2013l;
            varargout(1) = {lidar};
        end
    case 'dfc2013_right'
        D = load('Data/dfc2013_right');
        im = D.dfc2013h;
        imrgb = D.imrgb;
%         v = [53 40 15];
        if nargout==3
            lidar = D.dfc2013l;
            varargout(1) = {lidar};
        end
    case 'baden_full'
        D = load('Data/baden_full');
        im = D.im;
        imrgb = D.imrgb;
    case 'baden_topleft'
        D = load('Data/baden_topleft');
        im = D.im;
        imrgb = D.imrgb;
    case 'baden_topright'
        D = load('Data/baden_topright');
        im = D.im;
        imrgb = D.imrgb;        
    case 'baden_bottomleft'
        D = load('Data/baden_bottomleft');
        im = D.im;
        imrgb = D.imrgb;
    case 'baden_bottomright'
        D = load('Data/baden_bottomright');
        im = D.im;
        imrgb = D.imrgb;
    case 'plume'
        if nargin~=2
            error('wrong number of input argument. Please specify a frame number.')
        else
            framenumber = varargin{1};
            D = load('Data/Plume/vidrgb');
            imrgb = D.vidrgb(:,:,:,framenumber);
            framename = ['frameradiance_',int2str(framenumber)];
            F = load(strcat('Data/Plume/',framename));
            im = im2double(F.x);
        end
    otherwise
        error('there is no image with required name in folder ''Data''')
end
