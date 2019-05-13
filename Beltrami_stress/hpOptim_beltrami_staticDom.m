function [LogL,grad] = hpOptim_beltrami_staticDom(theta,Phi,lambdas,A,B,C,D,E,F,y,sig_m)
% [LogL,grad] = hpOptim_beltrami_staticDom(theta,Phi,lambdas,A,B,C,D,E,F,y,sig_m)

n = length(y);

if length(theta) ~= 24
    error('mode 3: theta = [A.sig_f,A.lx,A.ly,A.lz,B.sig_f,B.lx,B.ly,B.lz,...,F.sig_f,F.lx,F.ly,F.lz]')
end
A.sig_f = theta(1);
A.lx = theta(2);
A.ly = theta(3);
A.lz = theta(4);

B.sig_f = theta(5);
B.lx = theta(6);
B.ly = theta(7);
B.lz = theta(8);

C.sig_f = theta(9);
C.lx = theta(10);
C.ly = theta(11);
C.lz = theta(12);

D.sig_f = theta(13);
D.lx = theta(14);
D.ly = theta(15);
D.lz = theta(16);

E.sig_f = theta(17);
E.lx = theta(18);
E.ly = theta(19);
E.lz = theta(20);

F.sig_f = theta(21);
F.lx = theta(22);
F.ly = theta(23);
F.lz = theta(24);

l = theta([2:4,6:8,10:12,14:16,18:20,22:24]);

[~,nm] = size(lambdas);
nm = round(nm/6);
[SA] = calc_SLambda(A,lambdas(1,1:nm),lambdas(2,1:nm),lambdas(3,1:nm));
[SB] = calc_SLambda(B,lambdas(1,nm+1:2*nm),lambdas(2,nm+1:2*nm),lambdas(3,nm+1:2*nm));
[SC] = calc_SLambda(C,lambdas(1,2*nm+1:3*nm),lambdas(2,2*nm+1:3*nm),lambdas(3,2*nm+1:3*nm));
[SD] = calc_SLambda(D,lambdas(1,3*nm+1:4*nm),lambdas(2,3*nm+1:4*nm),lambdas(3,3*nm+1:4*nm));
[SE] = calc_SLambda(E,lambdas(1,4*nm+1:5*nm),lambdas(2,4*nm+1:5*nm),lambdas(3,4*nm+1:5*nm));
[SF] = calc_SLambda(F,lambdas(1,5*nm+1:6*nm),lambdas(2,5*nm+1:6*nm),lambdas(3,5*nm+1:6*nm));
SLambda = [SA,SB,SC,SD,SE,SF];



dif_stds = length(sig_m)>1;
[nrs, ncs] = size(sig_m);
if ncs ~= 1
    error('sig_m is wrong dimensions')
end

if dif_stds
    
    isnb = 1./sig_m.^2;
    
else
    isnb = ones(n,1)/sig_m^2;
    
end



[n,m] = size(Phi);          % total number of basies
nm = round(m/6);            % number of basies per basis function (A,B,C)
% solving using QR instead for numerical fuckign reasons
if nrs ~= n && dif_stds
    error('if differing sig_m s used, then must equal number of measurements')
end

if dif_stds
    Gamma = [Phi./sig_m;diag(1./sqrt(SLambda))];
else
    Gamma = [(ones(n,1)/sig_m).*Phi;diag(1./sqrt(SLambda))];
end
R = triu(qr(Gamma));
CZ = R(1:m,1:m);


logdetZ = 2*sum(log(abs(diag(CZ))));
if dif_stds
    logQ = sum(log(SLambda)) + sum(log(sig_m.^2)) + logdetZ;
else
    logQ = sum(log(SLambda)) + n*log(sig_m^2) + logdetZ;
end

optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
v = (linsolve(CZ,linsolve(CZ,Phi'*(isnb.*y),optsT),opts));
% yinvQy = y'*invSig_n*y - y'*invSig_n*Phi*(CZ\(CZ'\(Phi'*invSig_n*y)));  % slower
yinvQy = y'*(isnb.*y) - (y.*isnb)'*Phi*v;

LogL = logQ/2 + yinvQy/2 + n/2*log(2*pi);




if nargout == 2
    
    
    % GRADIENTS OF A
    inds = 1:nm;
    dlogQdsigf =2/A.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/A.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/A.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    %         dlogQdlx = sum(1/A.lx - A.lx*lambdas(1,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/A.lx - A.lx*lambdas(1,inds))./SLambda(inds)))));
    dlogQdlx = sum(1/A.lx - A.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/A.lx - A.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/A.lx - A.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    %         dlogQdly = sum(1/A.ly - A.ly*lambdas(2,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/A.ly - A.ly*lambdas(2,inds))./SLambda(inds)))));
    dlogQdly = sum(1/A.ly - A.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/A.ly - A.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/A.ly - A.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    %         dlogQdlz = sum(1/A.lz - A.lz*lambdas(3,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/A.lz - A.lz*lambdas(3,inds))./SLambda(inds)))));
    dlogQdlz = sum(1/A.lz - A.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/A.lz - A.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/A.lz - A.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradA = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    % GRADIENTS OF B
    inds = nm+1:2*nm;
    %         dlogQdsigf = 2/B.sig_f - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag(2/B.sig_f./SLambda(inds)))));
    dlogQdsigf = 2/B.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/B.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/B.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    %         dlogQdlx = sum(1/B.lx - B.lx*lambdas(1,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/B.lx - B.lx*lambdas(1,inds))./SLambda(inds)))));
    dlogQdlx = sum(1/B.lx - B.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/B.lx - B.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/B.lx - B.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    %         dlogQdly = sum(1/B.ly - B.ly*lambdas(2,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/B.ly - B.ly*lambdas(2,inds))./SLambda(inds)))));
    dlogQdly = sum(1/B.ly - B.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/B.ly - B.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/B.ly - B.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    %         dlogQdlz = sum(1/B.lz - B.lz*lambdas(3,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/B.lz - B.lz*lambdas(3,inds))./SLambda(inds)))));
    dlogQdlz = sum(1/B.lz - B.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/B.lz - B.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/B.lz - B.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradB = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    % GRADIENTS OF C
    inds = 2*nm+1:3*nm;
    %         dlogQdsigf = 2/C.sig_f - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag(2/C.sig_f./SLambda(inds)))));
    dlogQdsigf = 2/C.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/C.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/C.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    %         dlogQdlx = sum(1/C.lx - C.lx*lambdas(1,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds)))));
    dlogQdlx = sum(1/C.lx - C.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    %         dlogQdly = sum(1/C.ly - C.ly*lambdas(2,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds)))));
    dlogQdly = sum(1/C.ly - C.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    %         dlogQdlz = sum(1/C.lz - C.lz*lambdas(3,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/C.lz - C.lz*lambdas(3,inds))./SLambda(inds)))));
    dlogQdlz = sum(1/C.lz - C.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/C.lz - C.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/C.lz - C.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradC = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    % GRADIENTS OF D
    inds = 3*nm+1:4*nm;
    dlogQdsigf = 2/D.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/D.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/D.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    dlogQdlx = sum(1/D.lx - D.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/D.lx - D.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/D.lx - D.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    dlogQdly = sum(1/D.ly - D.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/D.ly - D.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/D.ly - D.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    dlogQdlz = sum(1/D.lz - D.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/D.lz - D.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/D.lz - D.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradD = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    % GRADIENTS OF E
    inds = 4*nm+1:5*nm;
    dlogQdsigf = 2/E.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/E.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/E.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    dlogQdlx = sum(1/E.lx - E.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/E.lx - E.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/E.lx - E.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    dlogQdly = sum(1/E.ly - E.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/E.ly - E.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/E.ly - E.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    dlogQdlz = sum(1/E.lz - E.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/E.lz - E.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/E.lz - E.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradE = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    % GRADIENTS OF F
    inds = 5*nm+1:6*nm;
    dlogQdsigf = 2/F.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/F.sig_f./SLambda(inds)),optsT),opts));
    dyinvQydsigf = -v(inds)'*diag(2/F.sig_f./SLambda(inds))*v(inds);
    grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
    
    dlogQdlx = sum(1/F.lx - F.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/F.lx - F.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlx = -v(inds)'*diag((1/F.lx - F.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
    grad_lx = dlogQdlx/2 + dyinvQydlx/2;
    
    dlogQdly = sum(1/F.ly - F.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/F.ly - F.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts));
    dyinvQydly = -v(inds)'*diag((1/F.ly - F.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
    grad_ly = dlogQdly/2 + dyinvQydly/2;
    
    dlogQdlz = sum(1/F.lz - F.lz*lambdas(3,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/F.lz - F.lz*lambdas(3,inds))./SLambda(inds)),optsT),opts));
    dyinvQydlz = -v(inds)'*diag((1/F.lz - F.lz*lambdas(3,inds))./SLambda(inds))*v(inds);
    grad_lz = dlogQdlz/2 + dyinvQydlz/2;
    
    gradF = [grad_sigf;grad_lx;grad_ly;grad_lz];
    
    grad = [gradA;gradB;gradC;gradD;gradE;gradF];
    
    
end



end