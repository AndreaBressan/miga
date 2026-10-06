function S = Lvolume() 
%LVOLUME Construct a B-spline basis for an L-shaped volume.
%
%   S = Lvolume()
%   constructs a three-dimensional B-spline basis representing a volume
%   with L-shaped geometry.
%
%   Output arguments:
%       S - Vector-valued B-spline basis defining the parametrization of
%           the L-shaped volume.


knots = {[0 0 0 1 1 1];[0 0 0 0.5 1 1 1];[0 0 1 1]};
coefs = reshape([-1 -1 0 1 -0.5 -1 0 1 0 -1 0 1 -1 1 0 1 -0.707106781186548 ...
                   0.292893218813452 0 1 0 0 0 1 -1 1 0 1 -0.292893218813452 ...
                   0.707106781186548 0 1 0 0 0 1 1 1 0 1 1 0.5 0 1 1 0 ...
                   0 1 -1 -1 1 1 -0.5 -1 1 1 0 -1 1 1 -1 1 1 1 -0.707106781186548 ...
                   0.292893218813452 1 1 0 0 1 1 -1 1 1 1 -0.292893218813452 ...
                   0.707106781186548 1 1 0 0 1 1 1 1 1 1 1 0.5 1 1 1 0 ...
                   1 1], 4, 3, 4, 2);
B=[bSpline(2,knots{1});bSpline(2,knots{2});bSpline(1,knots{3})];    
tpB = tBasis(B); 
coefs=coefs(1:3,:,:,:);
coefs=permute(coefs,[2,3,4,1]);
coefs=reshape(coefs,[],1,3);
S=subBasis(tpB,coefs);
end

