function nrm = scond(A,type)
%SCOND - Diagonal scaled condition number
%	nrm = scond(A,type) computes the condition number of
%	diagonally scaled A. The scaled A should have diagonals
% 	all ones.
%
%	Arguments:
%		A    : A symmetric matrix.
%		type : by default, is 2. One can pass any type as long as
%			it is recognized by MATLAB `cond`.
%
%	Outputs:
%		nrm  : the scaled condition number. 		

% Assert if only A is nonsymmetric!
% if ~issymmetric(A); error("Input should be symmetric!"); end

% parameters
n = size(A,1);
D = diag(diag(A).^(-1/2));

if nargin < 2
    nrm = cond(D*A*D);
else
    nrm = cond(D*A*D,type);
end

end
