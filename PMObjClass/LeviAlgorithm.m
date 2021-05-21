function [Model] = LeviAlgorithm(G,w,nb,na, varargin)

p = inputParser;

% Create variable names and assign default values after checking the value
addRequired(p,'G', @isnumeric);
addRequired(p,'w', @isnumeric);
addRequired(p,'nb', @isnumeric);
addRequired(p,'na', @isnumeric);


% Optional parameteres
addParameter(p,'domain','s');
addParameter(p,'fs',1);


% Re-parse parObj
parse(p,G,w,nb,na,varargin{:})

G = p.Results.G;
w = p.Results.w;
nb = p.Results.nb;
na = p.Results.na;
fs = p.Results.fs;

if p.Results.domain == 's' 
    v = 1i*w*fs;
else
    v = exp(-1i*w);
end


F = length(w);

%Construct regressor matrix
for kk = 1:F
    Reg(kk,:) = [-G(kk)*(v(kk).^(na:-1:1)) v(kk).^(nb:-1:0)];
end

Regt = [real(Reg);imag(Reg)];
Gt = [real(G);imag(G)];

theta = Lls(Regt,Gt);

Model.A = theta(1:na);
Model.B = theta(na+1:end);