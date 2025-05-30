function [V,D,NROT,NSWEEP,BOUND,SCOND] = mp_pjacobi(A, method)
%MP_PJACOBI - Mixed-precision preconditioned Jacobi algorithm
%	
%   Usage:
%      [V,D,NROT,NSWEEP,BOUND,SCOND] = MP_PJACOBI(A)
%   
%   Purpose:
%       MP_PJACOBI computes the spectral decomposition (SPD) of a real, 
%       N-by-N matrix A using the mixed-precision preconditioned Jacobi
%       algorithm. The SPD of A is written as 
%               A = V * diag(D) * V'
%       where V is an N-by-N orthogonal matrix, and D is an N dimensional
%       array whose entries are eigenvalues in descending order.
%       MP_PJACOBI can compute tiny eigenvalues with high relative accuracy
%       if A is symmetric positive definite. Even if A is ill--conditioned,
%       the MP_PJACOBI can still provide more digits of accuracy than eig. 
%       The algorithm ONLY works for DOUBLE PRECISION input, as the
%       algorithm requires a low precision spectral decomposition as a
%       preconditioner. We simulate the high precision using the Advanpix
%       Multiprecision Toolbox. 
%
%   Input:
%    - A is REAL matrix, dimension (N,N)
%       A is a SYMMETRIC, REAL matrix. The function will terminate if the
%       input is NONSYMMETRIC or Complex.
%	 - METHOD is a string
%		Default is "mp3", which indicates using high precision to
%		precondition the Jacobi algorithm.
%		The choice "mp2" uses working precision to precondition. 
%
%   Output:
%    - D is REAL vector, dimension (N)
%       If INFO = 0, D(i) is the ith largest eigenvalue of A.
%
%    - V is REAL matrix, dimension (N,N)
%       If INFO = 0 and EVEC = 1, V(:,i) is an normalized eigenvector 
%       corresponding to D(i). I.e. norm(V'*V - eye(N),2) = macheps, where 
%       macheps is the working precision unit roundoff. 
%       If INFO = 0 and EVEC = 0, V is an identity matrix. 
%   
%    - NROT is INTEGER 
%       The number of applied Jacobi rotation.
%
%    - NSWEEP is INTEGER 
%       The number of sweep the Jacobi rotation required to converge.
%       One sweep indicates the Jacobi algorithm has run through all A(i,j)
%       where i < j.
%
%   Author:
%       Zhengbo Zhou, Nov 2024, Manchester 
   
% Require eigenvectors 
EVEC = 1;

if nargin == 1 
    method = "mp3";
end

% Construct preconditioner at single precision
[Qlow,~] = eig(single(A));
[Qt,~] = qr(double(Qlow)); % HHQR
% [Qt,~] = mgs(double(Qlow)); % MGS
% [Qt,~] = ns(double(Qlow));
% temp1 = double(Qlow);
% I = eye(length(A));
% temp1 = 0.5*temp1*(3*I - temp1'*temp1);
% Qt = 0.5*temp1*(3*I - temp1'*temp1);

% Compute spectral decomposition
if method == "mp2"
    Atcomp = Qt'*A*Qt;
    Atcomp = (Atcomp + Atcomp')/2;
    [QJ,D,NROT,NSWEEP,INFO] = cjacobi(Atcomp,EVEC);
elseif method == "mp3"
    Athcomp = mp(Qt',34) * mp(A,34) * mp(Qt,34);
    Atcomp = double(Athcomp);
    Atcomp = (Atcomp + Atcomp')/2;
    [QJ,D,NROT,NSWEEP,INFO] = cjacobi(Atcomp,EVEC);
end


if (INFO == -1)
    warning("Routine cjacobi does not converge");
end

if nargout >= 5
    n = size(A,1);
    u = 2^(-53);
    p1 = sqrt(n); p2 = 7*n;
    SCOND = scond(Athcomp);
    BOUND = p2*u*SCOND+p1*u;
end

% Construct eigenvector
V = Qt * QJ;

% Sort eigenvalues in descending order 
[~,ind] = sort(D,'descend');
D = D(ind); 
V = V(:,ind);

end

