function featureAndPostion= featureExtractionCandidate(volume_image,xyzSpacing,noduleCandidateMask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   extract morphological and textural features for all nodule candidates
%   input: 
%           volume_image: original 3D image matrix
%           xyzSpacing: spacing information extracted from DICOM meta-information
%           noduleCandidateMask: 3D nodule candidates mask
%
%   output: detected features and position information for each nodule
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Shiwen Shen
%   Date: 09/28/2014
%   Email: shiwenshen@ucla.edu
%   Copy rignt: medical imaging informatics group, UCLA




[row,colomn,z]=size(volume_image);
sizeInfor=[row,colomn,z];


%image pixel normalization
meanValue=mean(volume_image(:));
stdValue=std(volume_image(:));
volume_image=(volume_image-meanValue)/stdValue;


roiCC=bwconncomp(noduleCandidateMask);
numberNodule=roiCC.NumObjects;
%final output array
featureAndPostion(numberNodule)=struct('feature','','postion','');
featureResult=zeros(1,27);

% loop to calculate features for each nodule
for i=1:numberNodule
    tempNoduleMask=false(sizeInfor);
    objectPosition=roiCC.PixelIdxList{i};
    tempNoduleMask(objectPosition)=1;
    
    %get the 2D midian slice for the nodule candidate
    [rowIn,coloIn,zIn]=ind2sub(sizeInfor,objectPosition);
     midZ=round((max(zIn)+min(zIn))*0.5);
     bina2D=tempNoduleMask(:,:,midZ);
     bina2DCC=bwconncomp(bina2D);
     numPixels = cellfun(@numel,bina2DCC.PixelIdxList);
     [largest1,idx1] = max(numPixels);
      bina2D= bina2D&0;
      bina2D(bina2DCC.PixelIdxList{idx1}) = 1;
     
     
     % calculate 2d geometric features
     featureResult(1:4)=GeometricFeature2D;
     % calculate 3d geometric features
     featureResult(5:12)=GeometricFeature3D;
     % calculate 2d intensity features
     featureResult(13:22)=intensityFeature2D;
     % calculate 3d intensity features
     featureResult(23:27)=intensityFeature3D;
     
     featureAndPostion(i).feature=featureResult;
     featureAndPostion(i).postion=objectPosition;
    
    
    
end
 

    function F1=GeometricFeature2D
    % calculate 2d geometric features
    % f1: Area
    % f2: Diameter
    % f3: Perimeter
    % f4: circularity
    % F1=[f1,f2,f3,f4]
    %bina2DCC=bwconncomp(bina2D);
    stats2D = regionprops(bina2D ,'Perimeter','Area','BoundingBox');   
    
    %area
    f1=stats2D.Area*xyzSpacing(1)*xyzSpacing(2);
    
    %Diameter
    xLength=(max(coloIn)-min(coloIn)+1)*xyzSpacing(1);
    yLength=(max(rowIn)-min(rowIn)+1)*xyzSpacing(2);
    f2=max([xLength,yLength]);
    
    %perimeter
    f3=imPerimeter(bina2D)*xyzSpacing(1);
    %circularity
    f4=4*pi*f1/(f3^2);
    
    F1=[f1,f2,f3,f4];
    
    
    end
    
    function F2=GeometricFeature3D
    % calculate 2d geometric features
    % f5: volume
    % f6: compactness; radius of  sphere with equivalent radius / root mean square distance
    %     >1 (sphere =1)
    % f7: x-y plane projection bounding box dimention rate: min(xlentgh,yLength)/xlentgh,yLength)
    % f8: 3d  bounding box dimention rate: min(xlentgh,yLength,zLength)/xlentgh,yLength,zLength)
    % f9: compactness 2: (Surface Area)^3/(volume)^2*36*pi
    % f10: Mean breadth
    % f11: Euler-Poincare Characteristic
    % f12: x-y plane projection compatness
    
    
    %volume
    numPix=length(objectPosition);
    f5=numPix*xyzSpacing(1)*xyzSpacing(2)*xyzSpacing(3);
    
    %compatness
    xLength=(max(coloIn)-min(coloIn)+1)*xyzSpacing(1);
    yLength=(max(rowIn)-min(rowIn)+1)*xyzSpacing(2);
    zLength=(max(zIn)-min(zIn)+1)*xyzSpacing(3);
    radius=max([xLength,yLength,zLength])*0.5;
    centerX=round((max(rowIn)+min(rowIn))*0.5);
    centerY=round((max(coloIn)+min(coloIn))*0.5);
    centerZ=midZ;
    rootMeanSquareDis=0;
    for kk=1:length(rowIn)
        dis=((rowIn(kk)-centerX)*xyzSpacing(1))^2+((coloIn(kk)-centerY)*xyzSpacing(2))^2+...
            ((zIn(kk)-centerZ)*xyzSpacing(3))^2;
        
        rootMeanSquareDis=rootMeanSquareDis+dis;
    end
    rootMeanSquareDis=rootMeanSquareDis/length(rowIn);
    rootMeanSquareDis=sqrt(rootMeanSquareDis);
    f6=radius/rootMeanSquareDis;
    
    %x-y plane projection bounding box dimention rate
    f7=min(xLength,yLength)/max(xLength,yLength);
    
    %3d  bounding box dimention rate
    f8=min([xLength,yLength,zLength])/max([xLength,yLength,zLength]);
    
    %compactness 2: (Surface Area)^3/(volume)^2*36*pi
    surfaceArea=imSurface(tempNoduleMask,[xyzSpacing(1),xyzSpacing(2),xyzSpacing(3)]);
    f9=surfaceArea^3/(f5^2*36*pi);
    
    
    
    %Mean breadth
    f10=imMeanBreadth(tempNoduleMask,[xyzSpacing(1),xyzSpacing(2),xyzSpacing(3)]);
    
    %Euler-Poincare Characteristic
    f11=imEuler3d(tempNoduleMask);
    
    
    %x-y plane projection compatness
    lengthArrow=length(rowIn);
    maskTem=false(size(tempNoduleMask,1),size(tempNoduleMask,2));
    for iii=1:lengthArrow
        maskTem(rowIn(iii),coloIn(iii))=1;
    end
    stats2D = regionprops(maskTem ,'Perimeter','Area','BoundingBox');
    areaT=stats2D.Area;
    perimeterT=stats2D.Perimeter;
    f12=4*pi*areaT/(perimeterT^2);
    
    F2=[f5,f6,f7,f8,f9,f10,f11,f12];
    end

    function F3=intensityFeature2D
    % calculate 2d intensity-based stasitical features
    % f13: minimum value inside
    % f14: mean contrast, (MeanInside-MeanOutside)/(MeanInside+MeanOutside)
    % f15: variance inside
    % f16: skewness inside
    % f17: kurtosis inside
    % f18-f22: moment 2,3,5,6,7, refer to: http://en.wikipedia.org/wiki/Image_moment
    
    medianSlice=volume_image(:,:,midZ);
    outBoudingSize=5;
    
    %minimum value inside
    f13=min(medianSlice(bina2D==1));
    
    % f14: mean contrast
    meanInside=mean(medianSlice(bina2D==1));
    indBin=find(bina2D==1);
    [row2D,col2D]=ind2sub([row,colomn],indBin);
    rowUp=min(row2D)-outBoudingSize;
    rowDown=max(row2D)+outBoudingSize;
    colLef=min(col2D)-outBoudingSize;
    colRig=max(col2D)+outBoudingSize;
    bina2DOut=false(size(bina2D));
    for ii=rowUp:rowDown
        for jj=colLef:colRig
            bina2DOut(ii,jj)=1;
        end
    end
    bina2DOut=bina2DOut&(~bina2D);
    meanOut=mean(medianSlice(bina2DOut==1));
    f14=(meanInside-meanOut)/(meanInside+meanOut);
    
    % variance inside
    f15=std(medianSlice(bina2D==1));
    
    % skewness, kurtosis inside
    f16=skewness(medianSlice(bina2D==1));
    f17=kurtosis(medianSlice(bina2D==1));
    
    %moment
%     M00=sum(medianSlice(bina2D==1));
%     
%     M10=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M10=M10+ii*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M01=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M01=M01+jj*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M11=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M11=M11+jj*ii*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M21=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M21=M21+ii^2*jj*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M12=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M12=M12+jj^2*ii*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M20=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M20=M20+ii^2*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M02=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M02=M02+jj^2*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M30=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M30=M30+ii^3*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     M03=0;
%     for ii=1:size(bina2D,1)
%         for jj=1:size(bina2D,2)
%             if bina2D(ii,jj)==1
%                 M03=M03+jj^3*medianSlice(ii,jj);
%             end
%         end
%     end
%     
%     xOverline=M10/M00;
%     yOverline=M01/M00;
%     u00=M00;
%     u11=M11-xOverline*M01;
%     u20=M20-xOverline*M10;
%     u02=M02-yOverline*M01;
%     u21=M21-2*xOverline*M11-yOverline*M20+2*(xOverline^2)*M01;
%     u12=M12-2*yOverline*M11-xOverline*M02+2*(yOverline^2)*M10;
%     u30=M30-3*xOverline*M20+2*(xOverline^2)*M10;
%     u03=M03-3*yOverline*M02+2*(yOverline^2)*M01;
%     gama11=u11/(u00^2);
%     gama20=u20/(u00^2);
%     gama02=u02/(u00^2);
%     gama21=u21/(u00^2.5);
%     gama12=u12/(u00^2.5);
%     gama30=u30/(u00^2.5);
%     gama03=u03/(u00^2.5);
%     I2=(gama20-gama02)^2+4*gama11^2;
%     I3=(gama30-3*gama12)^2+(3*gama21-gama03)^2;
%     I5=(gama30-3*gama12)*(gama30+gama12)*((gama30+gama12)^2-3*(gama21+gama03)^2)+...
%        (3*gama21-gama03)*(gama21+gama03)*(3*(gama30+gama12)^2-(gama21+gama03)^2);
%     I6=(gama20-gama02)*((gama30+gama12)^2-(gama21+gama03)^2)+4*gama11*(gama30+gama12)*(gama21+gama03);
%     I7=(3*gama21-gama03)*(gama30+gama12)*((gama30+gama12)^2-3*(gama21+gama03)^2)-...
%         (gama30-3*gama12)*(gama21+gama03)*(3*(gama30+gama12)^2-(gama21+gama03)^2);
    
     m00=0;
    m01=0;
    m10=0;
    m11=0;
    m12=0;
%     m21=0;
    [rowT,colT]=find(bina2D);
    for rr=1:length(rowT)
        for col=1:length(colT)
            m00=m00+medianSlice(rowT(rr),colT(col));
            m01=m01+(rowT(rr)^0)*(colT(col)^1)*medianSlice(rowT(rr),colT(col));
            m10=m10+(rowT(rr)^1)*(colT(col)^0)*medianSlice(rowT(rr),colT(col));
            m11=m11+(rowT(rr)^1)*(colT(col)^1)*medianSlice(rowT(rr),colT(col));
            m12=m12+(rowT(rr)^1)*(colT(col)^2)*medianSlice(rowT(rr),colT(col));
%             m21=m21+(rowT(rr)^2)*(colT(col)^1)*medianSlice(rr,col);
        end
    end
            
    
    f18=m01/m00;
    f19=m10/m00;
    f20=m11/m00;
    f21=m12/m00;
    f22=m00;
    F3=[f13,f14,f15,f16,f17,f18,f19,f20,f21,f22];
    
    end

    function F4=intensityFeature3D
    % calculate 3d intensity-based stasitical features
    % f23: minimum value inside
    % f24: mean contrast, (MeanInside-MeanOutside)/(MeanInside+MeanOutside)
    % f25: variance inside
    % f26: skewness inside
    % f27: kurtosis inside
    
    outBoudingSize=5;
    
    %minimum value inside
    f23=min(volume_image(tempNoduleMask==1));
    
    % f24: mean contrast
    meanInside=mean(volume_image(tempNoduleMask==1));
    rowUp=min(rowIn)-outBoudingSize;
    rowDown=max(rowIn)+outBoudingSize;
    colLef=min(coloIn)-outBoudingSize;
    colRig=max(coloIn)+outBoudingSize;
    zFr=min(zIn)-outBoudingSize;
    zBeh=max(zIn)+outBoudingSize;
    bina3DOut=false(size(volume_image));
    for ii=rowUp:rowDown
        for jj=colLef:colRig
            for kk=zFr:zBeh
                bina3DOut(ii,jj)=1;
            end
        end
    end
    bina3DOut=bina3DOut&(~tempNoduleMask);
    meanOut=mean(volume_image(bina3DOut==1));
    f24=(meanInside-meanOut)/(meanInside+meanOut);
    
    % variance inside
    f25=std(volume_image(tempNoduleMask==1));
    
    % skewness, kurtosis inside
    f26=skewness(volume_image(tempNoduleMask==1));
    f27=kurtosis(volume_image(tempNoduleMask==1));
    F4=[f23,f24,f25,f26,f27];
    
    end


end

