function [covFunc] = getCovFunc(basisParams)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if isfield(basisParams,'covFunc')
    covFunc = basisParams.covFunc;
    if (~strcmp('SE',covFunc) && ~strcmp('M5_2',covFunc) && ~strcmp('M3_2',covFunc) && ~strcmp('M1_2',covFunc))
        error('Invalid Covariance Function')
    end
else
    covFunc = 'SE';
end
end

