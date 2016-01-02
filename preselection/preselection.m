function noduleCandidateMask = preselection(candidateMsak,xyzSpacing)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   rule-based pruing, to remove most vessles and other non-nodule region
%
%   input: 
%           candidateMsak: input segmentation binary mask
%           xyzSpacing: spacing information extracted from DICOM meta-information
%
%   output: nodule candidates mask
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Shiwen Shen
%   Date: 09/28/2014
%   Email: shiwenshen@ucla.edu
%   Copy rignt: medical imaging informatics group, UCLA


[x,y,z]=size(candidateMsak);
%nodule candidate mask
noduleCandidateMask=false(x,y,z);
%vessle mask
vessMask=zeros(x,y,z);
%temperal mask for each object
binaryMask=false(x,y,z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% threshold values
diameTMin=5;%mm
diameTMax=30;
areaTMin=(diameTMin/2)^2*pi;
areaTMax=(diameTMax/2)^2*pi;
volumTMin=3*(diameTMin/2)^3*pi/4;
volumTMax=3*(diameTMax/2)^3*pi/4;
% overlapT=0.5;
elongationTMax=4;
circulTMin=1/6;



    sizeInfor=[x,y,z];
    roiImage=candidateMsak;
    roiCC=bwconncomp(roiImage,18);
    %objects loop, for objects in each region
    for j=1:roiCC.NumObjects
        binaryMask=false(x,y,z);
        objectPosition=roiCC.PixelIdxList{j};
        [rowIn,coloIn,zIn]=ind2sub(sizeInfor,objectPosition);
        if length(objectPosition)<5
            continue;
        end
        binaryMask(objectPosition)=1;
%         stats = regionprops(binaryMask ,'BoundingBox');
        
        %diameter
        xLength=(max(coloIn)-min(coloIn)+1)*xyzSpacing(1);
        yLength=(max(rowIn)-min(rowIn)+1)*xyzSpacing(2);
        zLength=(max(zIn)-min(zIn)+1)*xyzSpacing(3);
        diameter=max([xLength,yLength,zLength]);
        
        if diameter<diameTMin
            continue;
        end
        
        if diameter>diameTMax
            continue;
        end
        
        %enlongation
        enlongation=max([xLength,yLength,zLength])/min([xLength,yLength,zLength]);
        if enlongation>elongationTMax
            continue;
        end
        
        %volume
        numPix=size(find(roiCC.PixelIdxList{j}),1);
        volume=numPix*xyzSpacing(1)*xyzSpacing(2)*xyzSpacing(3);
        
        if volume>volumTMax
            continue;
        end
        
        if volume<volumTMin
            continue;
        end
        
        %area
        [rowIn,coloIn,zIn]=ind2sub([x,y,z],objectPosition);
        midZ=uint8((max(zIn)+min(zIn))*0.5);
        bina2D=binaryMask(:,:,midZ);
        bina2DCC=bwconncomp(bina2D);
        numPixels = cellfun(@numel,bina2DCC.PixelIdxList);
        [largest1,idx1] = max(numPixels);
        bina2D= bina2D&0;
        bina2D(bina2DCC.PixelIdxList{idx1}) = 1;
%         
%         
%         
%         if bina2DCC.NumObjects>1
%             continue;
%         end
%         
        stats2D = regionprops(bina2D ,'Perimeter','Area');
        area=stats2D.Area;
        %circularity
        % if the lesion tend to round, the value is closer to 1
        roundDegree = 4*pi * stats2D.Area / (stats2D.Perimeter^2) ;     
        
        
        if roundDegree>circulTMin
            if area>areaTMin
                if area<areaTMax
                    noduleCandidateMask(objectPosition)=1;
                end
            end
        end
       
        
    end
