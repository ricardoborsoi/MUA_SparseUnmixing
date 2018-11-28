function handle_IMout = drawborders(IMin,Seg,varargin)

borderwidth = 2;

if nargin == 2 % color is red by default and display in new figure
    color = [1 0 0];
    nofig = false;
elseif nargin == 3
    str = varargin{1};
    if strcmp(str,'nofig') % nofig option, color is red by default
        color = [1 0 0];
        nofig = true;
    elseif ischar(str) % color option, display in new figure by default
        nofig = false;
        color = rgb(str);
    elseif length(str)~=3
        error('Input color vector must have 3 components.')
    elseif length(str) == 3;
        nofig = false;
        color = str;
        color(color>1) = 1;
        color(color<0) = 0;
    else
        error('String argument must be either a color or ''nofig''.')
    end
elseif nargin == 4 % nofig option and color option
    str1 = varargin{1};
    str2 = varargin{2};
    if strcmp(str1,'nofig')
        nofig = true;
        strcol = str2;
    elseif strcmp(str2,'nofig')
        nofig = true;
        strcol = str1;
    else
        error('One input argument must be ''nofig''.')
    end 
    if ischar(strcol)
        color = rgb(strcol);
    else
        if length(strcol)~=3
            error('Input color vector must have 3 components.')
        else
            color = strcol;
            color(color>1) = 1;
            color(color<0) = 0;
        end
    end
else
    error('Wrong number of input arguments. Please check.')
end

labels = unique(Seg(:));
[lig col] = size(Seg);

if ~nofig
    handle_IMout = figure;
    imshow(IMin)
else
    handle_IMout = imshow(IMin);
end
hold on

if length(labels)==1 % unique region, draw border around the image
    
    plot(1:col,0.5*ones(size(1:col)),'Color',color,'LineWidth',borderwidth);
    plot(1:col,(lig+0.5)*ones(size(1:col)),'Color',color,'LineWidth',borderwidth);
    plot(0.5*ones(size(1:lig)),1:lig,'Color',color,'LineWidth',borderwidth);
    plot((col+0.5)*ones(size(1:lig)),1:lig,'Color',color,'LineWidth',borderwidth);
    
else
    for i=1:length(labels)
        
        Reg_i = Seg == labels(i);
        cont_i = imcontour(Reg_i,[0.5 0.5]);
        
        while ~isempty(cont_i)
            
            bdryLength = cont_i(2,1);
            plot(cont_i(1,2:bdryLength),cont_i(2,2:bdryLength),'Color',color,'LineWidth',borderwidth);
            cont_i(:,1:bdryLength+1) = [];
            
        end
        
    end
end

hold off
