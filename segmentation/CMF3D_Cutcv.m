function [ut, timet] = CMF3D_Cutcv(ur, lp, errb, ulab)

% We minimize wrt c_in and c_out, as well as lambda
% and now with adaptable errbound, having the effect that
% lambda is not solved to such a high accuracy until c1 and c2 start to
% converge.


%   Performing the continuous max-flow algorithm to solve the 
%   continuous min-cut problem in 3D
%    
%   Usage: [u, erriter, i, timet] = CMF3D_Cut;
%
%   Inputs: ur    -  data term (output of dataReorganize)
%           lp    -  length parameter. This parameters controls the level
%                       of fine detail in the result. Range: [1e-12, 1]
%           errb  -   [1x2] vector. The error bound for convergence.
%                    errb(1) error bound for convergence of u
%                    errb(2) error bound for convergence of c1 and c2
%           ulab  -  Initial estimates of background and foreground
%                       intensity
%
%
%   Outputs: 
%       - u: the final results u(x) in [0,1]. As the following paper,
%           the global binary result can be available by threshholding u
%           by any constant alpha in (0,1):
%
%           Nikolova, M.; Esedoglu, S.; Chan, T. F. 
%           Algorithms for Finding Global Minimizers of Image Segmentation and Denoising Models 
%           SIAM J. App. Math., 2006, 66, 1632-1648
%
%
%       - timet: gives the total computation time.
%       
%   Example:
%       >> [u, timet] = CMF3D_Cutcv(ur, lp, errb, ulab);

%   Please email noirinduggan@gmail.com for any questions, suggestions and bug reports
%
%   The Software is provided "as is", without warranty of any kind.
%
%




[rows, cols, heights]=size(ur);
szVol = rows*cols*heights;

% define the required parameters:

%
%   - cc: gives the step-size of the augmented Lagrangian method.
%       The optimal range of cc is [0.2, 3].

%
%   - numIter: the maximum iteration number.
%
%   - steps: the step-size for the graident-projection step to the
%       total-variation function. The optimal range of steps is [0.07,
%       0.12].
%



alpha =lp*ones(rows,cols,heights); % default  

cc = 0.35; % default = 0.35;
c_convergence = 3e-4;
steps = 0.11;
beta = 2;

% build up the priori L_2 data terms
Cs = abs(ur - ulab(1)).^beta;
Ct = abs(ur - ulab(2)).^beta;

% set the starting values
u = double((Cs-Ct) >= 0);
ps = min(Cs, Ct);
pt = ps;

pp1 = zeros(rows, cols+1, heights);
pp2 = zeros(rows+1, cols, heights);
pp3 = zeros(rows, cols, heights+1);
divp = - pp2(1:rows,:,:) + pp2(2:rows+1,:,:) + pp1(:,2:cols+1,:) ...
- pp1(:,1:cols,:) + pp3(:,:,2:heights+1) - pp3(:,:,1:heights);

erriter = zeros(300,1);
err_c = zeros(300,1);

tic

iter = 1;
err_c(1) = 1;
outerLoop=0;
innerloop=0;
while (err_c(iter) > errb(2))

    %  if cs and ct are beginning to converge, tighten the convergence limit
    % for lambda
    if err_c(iter) < c_convergence
        errb(1) = 5e-4;
    end

    i=1;
    erriter(1) = 1;
    outerLoop=outerLoop+1;
    printStr1=['outer loop number ',num2str(outerLoop)];
    disp(printStr1);
    while (erriter(i) > errb(1))
        innerloop=innerloop+1;
        printStr2=['inner loop number ',num2str(outerLoop)];
        disp(printStr2);

        i=i+1;

        % update the spatial flow field p = (pp1, pp2, pp3):
        %   the following steps are the gradient descent step with steps as the
        %   step-size.

        pts = divp - (ps - pt + u/cc);
        pp1(:,2:cols,:) = pp1(:,2:cols,:) + steps*(pts(:,2:cols,:) - pts(:,1:cols-1,:)); 
        pp2(2:rows,:,:) = pp2(2:rows,:,:) + steps*(pts(2:rows,:,:) - pts(1:rows-1,:,:));
        pp3(:,:,2:heights) = pp3(:,:,2:heights) + steps*(pts(:,:,2:heights) - pts(:,:,1:heights-1));

        % the following steps give the projection to make |p(x)| <= alpha(x)

        gk = sqrt((pp1(:,1:cols,:).^2 + pp1(:,2:cols+1,:).^2 + ...
        pp2(1:rows,:,:).^2 + pp2(2:rows+1,:,:).^2 + ...
        pp3(:,:,1:heights).^2 + pp3(:,:,2:heights+1).^2)*0.5);

        gk = double(gk <= alpha) + double(~(gk <= alpha)).*(gk ./ alpha);
        gk = 1 ./ gk;

        pp1(:,2:cols,:) = (0.5*(gk(:,2:cols,:) + gk(:,1:cols-1,:))).*pp1(:,2:cols,:); 
        pp2(2:rows,:,:) = (0.5*(gk(2:rows,:,:) + gk(1:rows-1,:,:))).*pp2(2:rows,:,:);
        pp3(:,:,2:heights) = (0.5*(gk(:,:,2:heights) + gk(:,:,1:heights-1))).*pp3(:,:,2:heights); 

        divp = - pp2(1:rows,:,:) + pp2(2:rows+1,:,:) + pp1(:,2:cols+1,:) ...
        - pp1(:,1:cols,:) + pp3(:,:,2:heights+1) - pp3(:,:,1:heights);

        % updata the source flow ps

        pts = divp - u/cc + pt + 1/cc;
        Cs = abs(ur - ulab(1)).^beta;
        ps = min(pts, Cs);

        % update the sink flow pt

        pts = - divp + ps + u/cc;
        Ct = abs(ur - ulab(2)).^beta;
        pt = min(pts, Ct);

        % update the multiplier u

        erru = cc*(divp + pt - ps);
        u = u - erru;

        % evaluate the avarage error

        erriter(i) = sum(sum(sum(abs(erru))))/szVol; 

        if (erriter(i) < errb(1))
            break;       
        end
    end


    ut = u > 0.5; 

    ulab_p = ulab;

    ulab(1) = sum(sum(sum( ur(:,:,:).*(1-ut(:,:,:)) ))./sum(sum( (1-ut(:,:,:)) )))/heights;
    ulab(2) = sum(sum(sum( ur(:,:,:).*ut(:,:,:) ))./sum(sum( ut(:,:,:) )))/heights;
    iter
    iter = iter + 1;
    err_c(iter) = abs(ulab_p(1) - ulab(1)) + abs(ulab_p(2) - ulab(2));
    fprintf('err_c = %f \r\n', err_c(iter,1));


end


toc
timet = toc

msg = sprintf('number of iterations = %u. \n', i);
disp(msg);

