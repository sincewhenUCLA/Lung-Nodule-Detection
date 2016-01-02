function viewBinaryMask(input3DArray)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   this function is used to view 3D image mask
%           Input:      3D image mask
%           Usage:      scroll mouse on image window to change slice
%                       or usage up arrow or down arrow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Author: Shiwen Shen
%   Date: 09/28/2014
%   Email: shiwenshen@ucla.edu
%   Copy rignt: medical imaging informatics group, UCLA




currenInd=1;
zth=size(input3DArray,3);
figure;
f= gcf;
imshow(input3DArray(:,:,currenInd),[]);
set(f,'KeyPressFcn',@(h_obj,evt) keymove(evt.Key));
set(f,'WindowScrollWheelFcn',@(h_obj,evt) keymove(evt.VerticalScrollCount));


function keymove(key)
    if strcmp(key,'uparrow') || sum(key)==-1 %If the uparrow is pressed or the mouse wheel is turned
        if ( currenInd<zth) 
            currenInd = currenInd+1;   
            imshow(input3DArray(:,:,currenInd),[]);
            currenInd
        end
    elseif strcmp(key,'downarrow') || sum(key)==1 %If the down arrow or mouse wheel is turned
        if (currenInd>1) 
            currenInd = currenInd-1;
            imshow(input3DArray(:,:,currenInd),[]);
            currenInd
        end
    end
end
end