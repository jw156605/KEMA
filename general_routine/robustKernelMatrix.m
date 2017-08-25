% function K = robustKernelMatrix(ker,X,X2,sigma)
%
% Inputs: 
%	ker:    'lin','poly','rbf','sam'
%	X:	data matrix with training samples in columns and features in rows
%	X2:	data matrix with test samples in columnsand features in rows
%	sigma: width of the RBF kernel
% 	b:     bias in the linear and polynomial kernel 
%	d:     degree in the polynomial kernel
%
% Output:
%	K: kernel matrix

% With Fast Computation of the RBF kernel matrix
% To speed up the computation, we exploit a decomposition of the Euclidean distance (norm)
%
% Gustavo Camps-Valls, 2006
% Jordi (jordi@uv.es),
%   2007-11: if/then -> switch, and fixed RBF kernel
%   2010-04: RBF can be computed now also on vectors with only one feature (ie: scalars)

function K = robustKernelMatrix(ker,X,X2,sigma,b)

[X_ica] = fastica(X,'numOfIC',30);

switch ker
    case 'lin'
          K = X_ica' * X_ica;

    case 'poly'
          K = (X_ica' * X_ica + b).^d;

    case 'rbf'
        if size(X_ica,1) == 1
            n1sq = X_ica.^2;
        else
            n1sq = sum(X_ica.^2);
        end
        n1 = size(X_ica,2);

        D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X_ica'*X_ica;
        K = exp(-D/(2*sigma^2));
        
    otherwise
        error(['Unsupported kernel ' ker])
end

