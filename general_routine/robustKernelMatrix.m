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

if exist('X2','var');
	csvwrite('temp.csv',horzcat(X,X2));
	system('Rscript /broad/macosko/jwelch/CellIntegration/KEMA/general_routine/ica.R temp.csv ica.mm 30','-echo');
	[ica] = csvread('ica.mm',1,1);
	ica = ica';
	n1 = size(X,2);
	n2 = size(X2,2);
	X_ica = ica(:,1:n1);
	X2_ica = ica(:,(n1+1):(n1+n2));
else
	csvwrite('temp.csv',X);
	system('Rscript /broad/macosko/jwelch/CellIntegration/KEMA/general_routine/ica.R temp.csv ica.mm 30','-echo');
	[ica] = csvread('ica.mm');
	ica = ica';
	X_ica = fastica(X,'numOfIC',30);
	X2_ica = X_ica;
end

switch ker
    case 'lin'
          K = X_ica' * X2_ica;

    case 'poly'
          K = (X_ica' * X2_ica + b).^d;

    case 'rbf'
        if size(X_ica,1) == 1
            n1sq = X_ica.^2;
        else
            n1sq = sum(X_ica.^2);
        end
        n1 = size(X_ica,2);

        if exist('X2','var');
            if size(X2_ica,1) == 1
                n2sq = X2_ica.^2;
            else
              n2sq = sum(X2_ica.^2);
            end
            n2 = size(X2_ica,2);
            D = (ones(n2,1)*n1sq)' + ones(n1,1)*n2sq -2*X_ica'*X2_ica;
        else
            D = (ones(n1,1)*n1sq)' + ones(n1,1)*n1sq -2*X_ica'*X_ica;
        end
        K = exp(-D/(2*sigma^2));
        
    otherwise
        error(['Unsupported kernel ' ker])
end

