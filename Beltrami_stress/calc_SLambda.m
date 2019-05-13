function [SLambda] = calc_SLambda(basis,lambda_x,lambda_y,lambda_z)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[covFunc] = getCovFunc(basis);
lx = basis.lx;
ly = basis.ly;
lz = basis.lz;
sig_f = basis.sig_f;

if strcmp('SE',covFunc)
    SLambda = sig_f^2*(2*pi)^(3/2)*lx*ly*lz*exp(-0.5*(lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2));
elseif strcmp('M5_2',covFunc)
    SLambda = sig_f^2*8*pi^(1.5)*gamma(4)/gamma(5/2)*5^(2.5)*lx*ly*lz./(5+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(4);
elseif strcmp('M3_2',covFunc)
    SLambda = sig_f^2*8*pi^(1.5)*gamma(3)/gamma(3/2)*3^(1.5)*lx*ly*lz./(3+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(3);
elseif strcmp('M1_2',covFunc)
    SLambda = sig_f^2*8*pi^(1.5)*gamma(2)/gamma(1/2)*lx*ly*lz./(1+lambda_x.^2*lx^2+lambda_y.^2*ly^2+lambda_z.^2*lz^2).^(2);
else
    error('Invalid covariance function')
end
end

