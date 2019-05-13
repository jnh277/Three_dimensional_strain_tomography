function [Phi_yI, SLambda,lambdas,Phi_T] = beltrami_approx_setDomain(num_basis,A,B,C,D,E,F,...
    nu,entry,exit,X,Y,Z,p_dom,nsegs)
% [Phi_yI, SLambda,lambdas,Phi_T] = beltrami_approx_setDomain(num_basis,A,B,C,D,E,F,nu,entry,exit,X,Y,Z,dom_p,nsegs)
% num_basis: nominal basis function resolution in each direction
% A,B,C,D,E,F: structures containing parameters for each underlying basis
% function
% X,Y,Z: test points
% nu: poissons ratio
% entry: ray entry points (each point is a column/or each column contains the nsegs(i) points)
% exit: ray exit points (each point is a column)
% nsegs: number of segments each ray intercepted

m = num_basis;

if ~exist('nsegs','var')
    nsegs = [];
end

[ACF] = getCovFunc(A);
% [BCF] = getCovFunc(B);  % assumning all the same
% [CCF] = getCovFunc(C);
% [DCF] = getCovFunc(D);
% [ECF] = getCovFunc(E);
% [FCF] = getCovFunc(F);

% generate basis function points
[mm1,mm2,mm3] = meshgrid(1:m);   % grid of basis functions  
insideCircle = sqrt(mm1.^2 + mm2.^2 + mm3.^2)/m <=1+eps;    % points inside a circle
mm1 = mm1(insideCircle);
mm2 = mm2(insideCircle);
mm3 = mm3(insideCircle);
MM = [mm1'; mm2';mm3'];   % all the basis functions, basis functions to change across columns
[~,mm_adj] = size(MM);


if isempty(nsegs)
    L = sqrt(sum((entry-exit).^2));
    nhat = (exit-entry)./L;
else
    [r,n] = size(entry);   % calculate distance from first entry to last exit
    ind =sub2ind([r,n],nsegs*3-2,1:n);
    xdif = exit(ind) - entry(1,:);
    ind = ind + 1;
    ydif = exit(ind) - entry(2,:);
    ind = ind + 1;
    zdif = exit(ind) - entry(3,:);
    L = sqrt(xdif.^2 + ydif.^2 + zdif.^2);
    nhat = [xdif;ydif;zdif]./L;
end

%% calculate Phi_yI

if nargout < 4
    % Note can only do this and make one call to basisFuncCom_setDomain if
    % all the ones are using the same covariance function
    [Idfuncs,lambdasA,SA] = basisFuncCom_setDomain(MM,ACF,A.lx,A.ly,A.lz,A.sig_f,nhat,entry,exit,[1 1 1 1 1 1],[],[],[],p_dom,nsegs);
    IdydyA = Idfuncs(:,:,2);
    IdzdzA = Idfuncs(:,:,3);
    IdydzA = Idfuncs(:,:,6);
    
    lambdasB = lambdasA;
    lambdasC = lambdasA;
    lambdasD = lambdasA;
    lambdasE = lambdasA;
    lambdasF = lambdasA;
    

%     [Idfuncs,lambdasB,SB] = basisFuncCom(MM,m,BCF,B.lx,B.ly,B.lz,B.sig_f,nhat,entry,exit,[1 0 1 0 1 0],[],[],[],nsegs);
    IdxdxB = Idfuncs(:,:,1);
    IdzdzB = Idfuncs(:,:,3);
    IdxdzB = Idfuncs(:,:,5);
    

%     [Idfuncs,lambdasC,SC] = basisFuncCom(MM,m,CCF,C.lx,C.ly,C.lz,C.sig_f,nhat,entry,exit,[1 1 0 1 0 0],[],[],[],nsegs);
    IdxdxC = Idfuncs(:,:,1);
    IdydyC = Idfuncs(:,:,2);
    IdxdyC = Idfuncs(:,:,4);
    
%     [Idfuncs,lambdasD,SD] = basisFuncCom(MM,m,DCF,D.lx,D.ly,D.lz,D.sig_f,nhat,entry,exit,[0 0 1 1 1 1],[],[],[],nsegs);
    IdzdzD = Idfuncs(:,:,3);        
    IdxdyD = Idfuncs(:,:,4);
    IdxdzD = Idfuncs(:,:,5);
    IdydzD = Idfuncs(:,:,6);
    
%     [Idfuncs,lambdasE,SE] = basisFuncCom(MM,m,ECF,E.lx,E.ly,E.lz,E.sig_f,nhat,entry,exit,[0 1 0 1 1 1],[],[],[],nsegs);
    IdydyE = Idfuncs(:,:,2);       
    IdxdyE = Idfuncs(:,:,4);
    IdxdzE = Idfuncs(:,:,5);
    IdydzE = Idfuncs(:,:,6);
    
%     [Idfuncs,lambdasF,SF] = basisFuncCom(MM,m,FCF,F.lx,F.ly,F.lz,F.sig_f,nhat,entry,exit,[1 0 0 1 1 1],[],[],[],nsegs);
    IdxdxF = Idfuncs(:,:,1);       
    IdxdyF = Idfuncs(:,:,4);
    IdxdzF = Idfuncs(:,:,5);
    IdydzF = Idfuncs(:,:,6);
    


else

    % Note can only do this and make one call to basisFuncCom_setDomain if
    % all the ones are using the same covariance function
    [Idfuncs,lambdasA,SA,dfuncs] = basisFuncCom_setDomain(MM,ACF,A.lx,A.ly,A.lz,A.sig_f,nhat,entry,exit,[1 1 1 1 1 1],X,Y,Z,p_dom,nsegs);
    IdydyA = Idfuncs(:,:,2);
    IdzdzA = Idfuncs(:,:,3);
    IdydzA = Idfuncs(:,:,6);
    dydyA = dfuncs(:,:,2);
    dzdzA = dfuncs(:,:,3);
    dydzA = dfuncs(:,:,6);
    
    lambdasB = lambdasA;
    lambdasC = lambdasA;
    lambdasD = lambdasA;
    lambdasE = lambdasA;
    lambdasF = lambdasA;
    

%     [Idfuncs,lambdasB,SB,dfuncs] = basisFuncCom(MM,m,BCF,B.lx,B.ly,B.lz,B.sig_f,nhat,entry,exit,[1 0 1 0 1 0],X,Y,Z,nsegs);
    IdxdxB = Idfuncs(:,:,1);
    IdzdzB = Idfuncs(:,:,3);
    IdxdzB = Idfuncs(:,:,5);
    dxdxB = dfuncs(:,:,1);
    dzdzB = dfuncs(:,:,3);
    dxdzB = dfuncs(:,:,5);

%     [Idfuncs,lambdasC,SC,dfuncs] = basisFuncCom(MM,m,CCF,C.lx,C.ly,C.lz,C.sig_f,nhat,entry,exit,[1 1 0 1 0 0],X,Y,Z,nsegs);
    IdxdxC = Idfuncs(:,:,1);
    IdydyC = Idfuncs(:,:,2);
    IdxdyC = Idfuncs(:,:,4);
    dxdxC = dfuncs(:,:,1);
    dydyC = dfuncs(:,:,2);
    dxdyC = dfuncs(:,:,4);   
    
%     [Idfuncs,lambdasD,SD,dfuncs] = basisFuncCom(MM,m,DCF,D.lx,D.ly,D.lz,D.sig_f,nhat,entry,exit,[0 0 1 1 1 1],X,Y,Z,nsegs);
    IdzdzD = Idfuncs(:,:,3);        
    IdxdyD = Idfuncs(:,:,4);
    IdxdzD = Idfuncs(:,:,5);
    IdydzD = Idfuncs(:,:,6);
    dzdzD = dfuncs(:,:,3);        
    dxdyD = dfuncs(:,:,4);
    dxdzD = dfuncs(:,:,5);
    dydzD = dfuncs(:,:,6);
    
%     [Idfuncs,lambdasE,SE,dfuncs] = basisFuncCom(MM,m,ECF,E.lx,E.ly,E.lz,E.sig_f,nhat,entry,exit,[0 1 0 1 1 1],X,Y,Z,nsegs);
    IdydyE = Idfuncs(:,:,2);       
    IdxdyE = Idfuncs(:,:,4);
    IdxdzE = Idfuncs(:,:,5);
    IdydzE = Idfuncs(:,:,6);
    dydyE = dfuncs(:,:,2);       
    dxdyE = dfuncs(:,:,4);
    dxdzE = dfuncs(:,:,5);
    dydzE = dfuncs(:,:,6);
    
%     [Idfuncs,lambdasF,SF,dfuncs] = basisFuncCom(MM,m,FCF,F.lx,F.ly,F.lz,F.sig_f,nhat,entry,exit,[1 0 0 1 1 1],X,Y,Z,nsegs);
    IdxdxF = Idfuncs(:,:,1);       
    IdxdyF = Idfuncs(:,:,4);
    IdxdzF = Idfuncs(:,:,5);
    IdydzF = Idfuncs(:,:,6);
    dxdxF = dfuncs(:,:,1);       
    dxdyF = dfuncs(:,:,4);
    dxdzF = dfuncs(:,:,5);
    dydzF = dfuncs(:,:,6);
end

N1N1 = nhat(1,:)'.^2;
N2N2 = nhat(2,:)'.^2;
N3N3 = nhat(3,:)'.^2;
N1N2 = 2*nhat(1,:)'.*nhat(2,:)';
N1N3 = 2*nhat(1,:)'.*nhat(3,:)';
N2N3 = 2*nhat(2,:)'.*nhat(3,:)';

phi_yA = N1N1.*(-nu*IdzdzA-nu*IdydyA)+N2N2.*(IdzdzA-nu*IdydyA)+N3N3.*(IdydyA-nu*IdzdzA)...
    +N2N3.*(-(1+nu)*IdydzA); % relating measurements to the 'A' basis functions
phi_yB = N1N1.*(IdzdzB-nu*IdxdxB)+N2N2.*(-nu*IdzdzB-nu*IdxdxB)+N3N3.*(IdxdxB-nu*IdzdzB)...
    +N1N3.*(-(1+nu)*IdxdzB); % relating to the 'B' basis functions
phi_yC = N1N1.*(IdydyC-nu*IdxdxC)+N2N2.*(IdxdxC-nu*IdydyC)+N3N3.*(-nu*IdydyC-nu*IdxdxC)...
    +N1N2.*(-(1+nu)*IdxdyC); % relating to the 'C' basis functions
phi_yD = N1N1.*(2*nu*IdxdyD)+N2N2.*(2*nu*IdxdyD)+N3N3.*(-2*IdxdyD)...
    +N1N2.*(-(1+nu)*IdzdzD)+N1N3.*((1+nu)*IdxdzD)+N2N3.*((1+nu)*IdydzD);    % relating to the 'D' basis functions
phi_yE = N1N1.*(2*nu*IdxdzE)+N2N2.*(-2*IdxdzE)+N3N3.*(2*nu*IdxdzE)...
    +N1N2.*((1+nu)*IdydzE)+N1N3.*((1+nu)*IdxdyE)+N2N3.*(-(1+nu)*IdydyE); 
phi_yF = N1N1.*(-2*IdydzF)+N2N2.*(2*nu*IdydzF)+N3N3.*(2*nu*IdydzF)...
    +N1N2.*((1+nu)*IdxdzF)+N1N3.*(-(1+nu)*IdxdxF)+N2N3.*((1+nu)*IdxdyF); 

Phi_yI =[phi_yA,phi_yB,phi_yC,phi_yD,phi_yE,phi_yF]./L';

if nargout > 3
    np = length(X(:));
    % and combined for speed
    % interweave the Phis, so that we have [exx,eyy,ezz,exy,exz,ezz] for each
    % point
    Phi_T = NaN(np*6,mm_adj*6);
    Phi_T(1:6:end,:) = [-nu*dzdzA-nu*dydyA, dzdzB-nu*dxdxB, dydyC-nu*dxdxC, 2*nu*dxdyD, 2*nu*dxdzE, -2*dydzF]; % phi_exx;
    Phi_T(2:6:end,:) = [dzdzA-nu*dydyA, -nu*dzdzB-nu*dxdxB, dxdxC-nu*dydyC, 2*nu*dxdyD, -2*dxdzE, 2*nu*dydzF]; % phi_eyy;
    Phi_T(3:6:end,:) = [dydyA-nu*dzdzA, dxdxB-nu*dzdzB, -nu*dydyC-nu*dxdxC, -2*dxdyD, 2*nu*dxdzE, 2*nu*dydzF]; % phi_ezz;
    Phi_T(4:6:end,:) = [zeros(np,mm_adj), zeros(np,mm_adj), -(1+nu)*dxdyC, -(1+nu)*dzdzD, (1+nu)*dydzE, (1+nu)*dxdzF];  % phi_exy;
    Phi_T(5:6:end,:) = [zeros(np,mm_adj), -(1+nu)*dxdzB, zeros(np,mm_adj), (1+nu)*dxdzD, (1+nu)*dxdyE, -(1+nu)*dxdxF];  % phi_exz;
    Phi_T(6:6:end,:) = [-(1+nu)*dydzA, zeros(np,mm_adj), zeros(np,mm_adj), (1+nu)*dydzD, -(1+nu)*dydyE, (1+nu)*dxdyF];  % phi_eyz;
end




% SA was computed earlier
[SB] = calc_SLambda(B,lambdasB(1,:),lambdasB(2,:),lambdasB(3,:));
[SC] = calc_SLambda(C,lambdasC(1,:),lambdasC(2,:),lambdasC(3,:));
[SD] = calc_SLambda(D,lambdasD(1,:),lambdasD(2,:),lambdasD(3,:));
[SE] = calc_SLambda(E,lambdasE(1,:),lambdasE(2,:),lambdasE(3,:));
[SF] = calc_SLambda(F,lambdasF(1,:),lambdasF(2,:),lambdasF(3,:));

SLambda = [SA,SB,SC,SD,SE,SF];
lambdas = [lambdasA,lambdasB,lambdasC,lambdasD,lambdasE,lambdasF];


end



