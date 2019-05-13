function [Idfuncs,lambdas,SLambda,dfuncs] = basisFuncCom_setDomain(MM,covFunc,...
    lx,ly,lz,sig_f,nhat,entry,exit,req,X,Y,Z,p_dom,nsegs)
% [Idfuncs,lambdas,SLambda,dfuncs] = basisFuncCom_setDomain(MM,covFunc,lx,ly,lz,sig_f,nhat,entry,exit,req,X,Y,Z,dom_p,nsegs)
%   m: number of basis functions in each direction
%   MM: placing of basis functions to make a circular placement rather than
%   lx,ly,lz: length scales
%   sig_f: prior uncertainty
%   nhat: ray directions, each column corresponds to a ray
%   entry: entry intersects for a ray, each column corresponds to a ray
%   exit: exit intersects for a ray, each column corresponds to a ray
%   req: required ray integrals and derivatives
%   X,Y,Z: position of test points
%   nsegs: number of segments each ray intercepts

[~,n] = size(entry);        % number of ray measurements
[~,mm_adj] = size(MM);    

if ~exist('nsegs','var') || isempty(nsegs)
    nsegs = ones(1,n);
end

%% determine the problem domain and place the basis function lambdas

Lx = p_dom(1);            % Domain scales
Ly = p_dom(2);
Lz = p_dom(3);

lambda_x = MM(1,:)*pi/2/Lx;
lambda_y = MM(2,:)*pi/2/Ly;
lambda_z = MM(3,:)*pi/2/Lz;


n1 = nhat(1,:)'; n2 = nhat(2,:)'; n3 = nhat(3,:)';

n1lx = bsxfun(@times,n1,lambda_x);
n2ly = bsxfun(@times,n2,lambda_y);
n3lz = bsxfun(@times,n3,lambda_z);
Q1 = (n1lx-n2ly-n3lz);
Q2 = (n1lx+n2ly-n3lz);
Q3 = (n1lx-n2ly+n3lz);
Q4 = (n1lx+n2ly+n3lz);

I1 = ceil(find(abs(Q1(:)) <= eps^(1/2))/length(n1));
I2 = ceil(find(abs(Q2(:)) <= eps^(1/2))/length(n1));
I3 = ceil(find(abs(Q3(:)) <= eps^(1/2))/length(n1));
I4 = ceil(find(abs(Q4(:)) <= eps^(1/2))/length(n1));
I = unique([I1;I2;I3;I4]);

% attempt to avoid divide by zero (or very small number) by moving the
% basis functions
while (~isempty(I))
    lambda_x(I) = lambda_x(I) + eps^(1/2)*max(lambda_x)*(1-rand(1,length(I)));
    lambda_y(I) = lambda_y(I) + eps^(1/2)*max(lambda_y)*(1-rand(1,length(I)));
    lambda_z(I) = lambda_z(I) + eps^(1/2)*max(lambda_z)*(1-rand(1,length(I)));
    
    n1lx = bsxfun(@times,n1,lambda_x);
    n2ly = bsxfun(@times,n2,lambda_y);
    n3lz = bsxfun(@times,n3,lambda_z);
    Q1 = (n1lx-n2ly-n3lz);
    Q2 = (n1lx+n2ly-n3lz);
    Q3 = (n1lx-n2ly+n3lz);
    Q4 = (n1lx+n2ly+n3lz);
    
    I1 = ceil(find(Q1(:) == 0)/length(n1));
    I2 = ceil(find(Q2(:) == 0)/length(n1));
    I3 = ceil(find(Q3(:) == 0)/length(n1));
    I4 = ceil(find(Q4(:) == 0)/length(n1));
    I = unique([I1;I2;I3;I4]);
end

% put the lambdas together for output
lambdas = [lambda_x;lambda_y;lambda_z];
% compute the spectral intensities of the basis functions
SLambda = sig_f^2*(2*pi)^(3/2)*lx*ly*lz*exp(-0.5*(lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2));
% if strcmp('SE',covFunc)
%     SLambda = sig_f^2*(2*pi)^(3/2)*lx*ly*lz*exp(-0.5*(lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2));
% elseif strcmp('M5_2',covFunc)
%     SLambda = sig_f^2*8*pi^(3/2)*gamma(4)/gamma(5/2)*5^(2.5)*lx*ly*lz./(5+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(4);
% elseif strcmp('M3_2',covFunc)
%     SLambda = sig_f^2*8*pi^(3/2)*gamma(3)/gamma(3/2)*3^(1.5)*lx*ly*lz./(3+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(3);
% elseif strcmp('M1_2',covFunc)
%     SLambda = sig_f^2*8*pi^(3/2)*gamma(2)/gamma(1/2)*lx*ly*lz./(1+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(2);
% else
%     error('Invalid covariance function')
% end

% preallocate to zeros so multiseg stuff works  
Idfuncs = zeros(n,mm_adj,6);


%% Do the ray integrals of the function derivatives

for ss = 1:max(nsegs)
    segInd = find(nsegs >= ss);        % index of measurements that have at least this many segments
    
    x0 = entry(ss*3-2,segInd)'; y0 = entry(ss*3-1,segInd)'; z0 = entry(ss*3,segInd)';
    alpha_x = bsxfun(@times,Lx+x0,lambda_x);
    alpha_y = bsxfun(@times,Ly+y0,lambda_y);
    alpha_z = bsxfun(@times,Lz+z0,lambda_z);

    Gamma1s = cos(alpha_x-alpha_y-alpha_z);
    Gamma2s = cos(alpha_x+alpha_y-alpha_z);
    Gamma3s = cos(alpha_x-alpha_y+alpha_z);
    Gamma4s = cos(alpha_x+alpha_y+alpha_z);

    % exit points s = L
    xf = exit(ss*3-2,segInd)'; yf = exit(ss*3-1,segInd)'; zf = exit(ss*3,segInd)';

    alpha_x = bsxfun(@times,Lx+xf,lambda_x);
    alpha_y = bsxfun(@times,Ly+yf,lambda_y);
    alpha_z = bsxfun(@times,Lz+zf,lambda_z);

    Gamma1f = cos(alpha_x-alpha_y-alpha_z);
    Gamma2f = cos(alpha_x+alpha_y-alpha_z);
    Gamma3f = cos(alpha_x-alpha_y+alpha_z);
    Gamma4f = cos(alpha_x+alpha_y+alpha_z);

    % merging
    Gamma1 = (Gamma1f-Gamma1s)./Q1(segInd,:);
    Gamma2 = (Gamma2f-Gamma2s)./Q2(segInd,:);
    Gamma3 = (Gamma3f-Gamma3s)./Q3(segInd,:);
    Gamma4 = (Gamma4f-Gamma4s)./Q4(segInd,:);

    if req(1) % integral of d/dxdx
        Idfuncs(segInd,:,1) = Idfuncs(segInd,:,1) - (lambda_x.^2/sqrt(Lx*Ly*Lz)/4).*(Gamma1-Gamma2-Gamma3+Gamma4);  
    end
    if req(2) % integral of d/dydy
        Idfuncs(segInd,:,2) = Idfuncs(segInd,:,2) -(lambda_y.^2/sqrt(Lx*Ly*Lz)/4).*(Gamma1-Gamma2-Gamma3+Gamma4);
    end
    if req(3) % integral of d/dzdz 
        Idfuncs(segInd,:,3) = Idfuncs(segInd,:,3) -(lambda_z.^2/sqrt(Lx*Ly*Lz)/4).*(Gamma1-Gamma2-Gamma3+Gamma4);
    end
    if req(4) % integral of d/dxdy
        Idfuncs(segInd,:,4) = Idfuncs(segInd,:,4) +(lambda_x.*lambda_y/sqrt(Lx*Ly*Lz)/4).*(Gamma1+Gamma2-Gamma3-Gamma4);
    end
    if req(5)   % integral of d/dxdz
        Idfuncs(segInd,:,5) = Idfuncs(segInd,:,5) +(lambda_x.*lambda_z/sqrt(Lx*Ly*Lz)/4).*(Gamma1-Gamma2+Gamma3-Gamma4);
    end
    if req(6)   % integral of d/dydz
        Idfuncs(segInd,:,6) = Idfuncs(segInd,:,6) +(lambda_y.*lambda_z/sqrt(Lx*Ly*Lz)/4).*(-Gamma1-Gamma2-Gamma3-Gamma4);
    end
end
%% Compute the basis function derivates
if nargout > 3
    n_test = length(X(:));
    dfuncs = NaN(n_test,mm_adj,6);
    
    Bx = bsxfun(@times,X(:)+Lx,lambda_x);       % here X,Y,Z are for the test points
    By = bsxfun(@times,Y(:)+Ly,lambda_y);
    Bz = bsxfun(@times,Z(:)+Lz,lambda_z);

    sx = sin(Bx);       % avoid calling the sin function lots
    sy = sin(By);
    sz = sin(Bz);
    cx = cos(Bx);
    cy = cos(By);
    cz = cos(Bz);


    if req(1)   % d/dxdx
        dfuncs(:,:,1) =  -(lambda_x.^2/sqrt(Lx*Ly*Lz)) .* sx .* sy .* sz;
    end
    if req(2)   % d/dydy
        dfuncs(:,:,2) =  -(lambda_y.^2/sqrt(Lx*Ly*Lz)) .* sx .* sy .* sz;
    end
    if req(3)   % d/dzdz
        dfuncs(:,:,3) =  -(lambda_z.^2/sqrt(Lx*Ly*Lz)) .* sx .* sy .* sz;
    end
    if req(4)   % d/dxdy
        dfuncs(:,:,4) =  ((lambda_x.*lambda_y)/sqrt(Lx*Ly*Lz)) .* cx .* cy .* sz;
    end
    if req(5)   % d/dxdz
        dfuncs(:,:,5) =  ((lambda_x.*lambda_z)/sqrt(Lx*Ly*Lz)) .* cx .* sy .* cz;
    end
    if req(6)   % d/dydz
        dfuncs(:,:,6) =  ((lambda_y.*lambda_z)/sqrt(Lx*Ly*Lz)) .* sx .* cy .* cz;
    end    

end
    
end