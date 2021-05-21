function [V, J] = SigmoidFcn(theta,x)
% Evaluates the sigmoid function V = Ax/sqrt(1+Bx^2), where A = theta(1)
% and B = theta(2).
%
% Input arguments:
%   theta: Vector with values A and B, size 2 x 1
%   x: Vector of values where function is evaluated, size n x 1
%
% Output arguments:
%   V: Sigmoid function values, size n x 1
%   J: Jacobian of function for parameter optimisation, size n x 2
%
% W. D. Widanage 01/07/2014 (generous)

x = x(:);
A = theta(1);
B = theta(2);
V = A*x./sqrt(1+B*x.^2);

d1 = (1+B*x.^2).^(3/2);
J = [x./sqrt(1+B*x.^2), (-0.5*A*x.^3)./d1];