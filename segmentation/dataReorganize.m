function [volume_image,sliceLocationArray,xyzSpacing]=dataReorganize(directory)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   DATAREORGANIZE reorganize the dicom image in the directory
%   using correct order and store the sorted images in a 3-D array
%
%   Input:     
%              directory: the image diretory
%   Output: 
%              volume_image: the sorted image array
%              sliceLocationArray: the reordered image's relative location 
%                                  value
%              xyzSpacing: spacing between the center of adjacent pixels; 
%                          [row spacing,column spacing, z spacing]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Shiwen Shen
%   Date: 09/28/2014
%   Email: shiwenshen@ucla.edu
%   Copy rignt: medical imaging informatics group, UCLA


cd(directory);
d = ls('*.dcm');
m = size(d,1);

%get the spacing information
xyzSpacing=zeros(3,1);%row spacing, column spacing and z spacing
metadata = dicominfo(d(1,:));
[group, element] = dicomlookup('PixelSpacing');
xySpacing=metadata.(dicomlookup(group, element));
xyzSpacing(1:2)=xySpacing;
[group, element] = dicomlookup('SliceThickness');
zSpacing=metadata.(dicomlookup(group, element));
xyzSpacing(3)=zSpacing;


%get the reordered image data and position information
[group, element] = dicomlookup('InstanceNumber');
[group1, element1] = dicomlookup('SliceLocation');
sdata(m) = struct('imagename','','instance',0,'image','','SliceLocation','');

for i = 1:m
    metadata = dicominfo(d(i,:));
    position = metadata.(dicomlookup(group, element));
    imageArray=dicomread(d(i,:));
    SliceLocation= metadata.(dicomlookup(group1, element1));
    sdata(i) = struct('imagename',d(i,:),'instance',position,'image',imageArray,'SliceLocation',SliceLocation);
end

[sortedValue order] = sort([sdata(:).instance],'ascend');
sorted = sdata(order).';
[length,width]=size(imageArray);
volume_image=zeros(length,width,m);
sliceLocationArray=zeros(m,1);

for i = 1:m
    volume_image(:,:,i) = sorted(i).image;
    sliceLocationArray(i)=sorted(i).SliceLocation; 
end

end

