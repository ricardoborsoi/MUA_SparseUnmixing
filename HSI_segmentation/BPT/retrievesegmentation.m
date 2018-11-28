function output = retrievesegmentation(prunedtree,initsegmap,string,varargin)

% segout = retrievesegmentation(prunedtree,initsegmap,stringoption)
% stringoption: 'regions' or 'borders'
% If stringoption = 'border', possibility of a 4th argument being the image
% on which the borders will be superimposed.

segmap = zeros(size(initsegmap));

PosNodes = [prunedtree.pruning];
node = prunedtree(logical(PosNodes));

for i=1:length(node)
    
    if isempty(node(i).nodeinfo.leaves)
        % the considered node is a leaf
        segmap(initsegmap==node(i).label) = node(i).label;
    else
        ind = ismember(initsegmap,node(i).nodeinfo.leaves);
        segmap(ind) = node(i).label;
    end
end

if strcmp(string,'regions')
    
    output = segmap;
    
elseif strcmp(string,'borders')
    
    reglbl = unique(segmap(:));
    borders = zeros(size(segmap));
    for i=1:length(reglbl)
        reg_i = double(segmap==reglbl(i));
        borders = borders + edge(reg_i,'Canny',0.95);
    end
    borders = borders|bwmorph(borders,'majority');
    borders = bwmorph(borders,'thin',Inf);
    borders = borders==0;
    
    if nargin == 3
        
        output = borders;
        
    elseif nargin == 4
        
        im = varargin{1};
        if size(im,1)==size(borders,1)&&size(im,2)==size(borders,2)
            
            if size(im,3)==3
                imgray = rgb2gray(im);
            else
                imgray = im(:,:,1);
            end
            white = imgray>=0.5;
            black = imcomplement(white);
            
            indwhite = borders==0&white;
            indblack = borders==0&black;
            indnwhite = repmat(indwhite,[1 1 size(im,3)]);
            indnblack = repmat(indblack,[1 1 size(im,3)]);
            im(indnwhite==1) = 0;
            im(indnblack==1) = 1;
            output = im;
            
        else
            error('The dimensions of the input image don''t match the segmentation map')
        end
        
    else
        error('Wrong number of input arguments. Please check')
    end
    
else
    error('Input string mismatches. Please choose between ''regions'' or ''borders''')
end
