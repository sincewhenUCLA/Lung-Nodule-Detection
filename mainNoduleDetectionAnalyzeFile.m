%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           CDSC: lung nodule detection pipeline
% This is the main function to run the nodule detection
% pipeline demo. The input is low-dose (25% projections) CT
% reconstructed image using EM+TV method and the final output 
% is the detected nodule mask. This code is only for demo 
% purpose and some parts have been modified for the specific
% low-dose CT image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shiwen Shen
% Date: 09/28/2014
% Email: shiwenshen@ucla.edu
% Copy rignt: medical imaging informatics group, UCLA


clc;
clear;
close all;
wkdir = pwd;



    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%imgae read
info = analyze75info('PAT1_CT_RECON_FBP.hdr'); % Replace with path of Analyze 7.5 image file
xyzSpacing=[0.7;0.7;1.25];
volume_image=analyze75read('PAT1_CT_RECON_FBP');
volume_image = flipdim(volume_image,1);
volume_image=double(volume_image);

%%%%%%%%%%%%%%%%%%%%%%%
%get two phase segmentation result
% lp = 1e-12; % increase this value for less fine detail (maximum tested value = 2)
lp = 1e-13; % increase this value for less fine detail (maximum tested value = 2)
errb = [1e-1,5e-4];  
% ulab = [0.3, 0.6];
% ulab = [0.03, 0.25];
ulab = [0.1, 0.4];
%data clamp
upperBand=80;
lowwerBand=0;
volume_image(volume_image > upperBand) = upperBand;      
volume_image(volume_image < lowwerBand) = lowwerBand;
% scale data
volume_image = (volume_image-min(volume_image(:)))./(max(volume_image(:))-min(volume_image(:)));
[intialSegResult, timet] = CMF3D_Cutcv(volume_image, lp, errb, ulab);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    generate lung mask and candidates   %%%%%%
maskImageVolume= segmentationMask(intialSegResult, volume_image);
candidateMsak=maskImageVolume&intialSegResult;
candidateMsak(:,:,1:4)=0;
candidateMsak(:,:,end-4:end)=0;


%%%%%%%%%%%%%%%%%%%%%%%%%
% rule-based removal
noduleCandidateMask = preselection(candidateMsak,xyzSpacing);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%feature extraction
featureAndPostion= featureExtractionCandidate(volume_image,xyzSpacing,noduleCandidateMask);
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %classification
featureMask=[0,1,1,1,0,1,1,0,1,0,1,1,0,1,0,0,0,0,0,0,1,1,1,0,1,1,0];%select features for classifier

% %%kmeans classification preparation
load SelectedNegativesamples.mat
load SelectedPositivesamples.mat
load meanFeature.mat
load stdFeature.mat

meanP=mean(SelectedPositivesamples,1); %central points for postive training samples
meanN=mean(SelectedNegativesamples,1); %central pointes for negtive training samples 

meanP=meanP(featureMask==1);
meanN=meanN(featureMask==1);
meanFeature=meanFeature(featureMask==1);
stdFeature=stdFeature(featureMask==1);

finalNodule=false(size(volume_image)); %final output 
noduleNum=0;

%%%%%%%%%%%%%%%%%%%%%%
%nodule volume size
noduleVolume=[];

%%%%%%%%%%%%%%%%%%
%Z indexes of the nodule slices (which slices are detected with nodules) 
zIndexOfNoudle=[];

for jj=1:length(featureAndPostion)
    featureTemple=featureAndPostion(jj).feature;
    featureTemple=featureTemple(featureMask==1);
    featureTemple=(featureTemple-meanFeature)./stdFeature;
    temp1 = sum((featureTemple - meanP).^2);
    temp2 = sum((featureTemple - meanN).^2);
    
    if temp1<temp2
        finalNodule(featureAndPostion(jj).postion)=1;
        noduleNum=noduleNum+1;
        temNoduleVolume=length(featureAndPostion(jj).postion)*xyzSpacing(1)*xyzSpacing(2)*xyzSpacing(3);
        noduleVolume=[noduleVolume,temNoduleVolume];
    end
end
viewBinaryMask(finalNodule);



cd(wkdir);


    