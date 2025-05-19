function [Q,R] = HouseQR(A)
%  Mixed precision Householder QR
% function [Q_fact,R] = HouseQR(A)
% Householder method for computing the QR factorization A = QR.
% A is mxn with m>=n.
% Q_fact is mxn and  houses the factored form representation of the 
%    orthogonal factor Q.
% R is nxn and upper triangular with nonnegative diagonal entries 
%    so A = QR.
% GVL4: Algorithm 5.2.1
A = single(A);
[m,n] = size(A);
for j=1:min(n,m-1)
    % Compute the jth Householder matrix Q{j}...
    [v,beta] = House(A(j:m,j));
    % Update...
    % A = Q{j}A
    A(j:m,j:n) = A(j:m,j:n) - (beta*v)*(v'*A(j:m,j:n));
    % Save the Householder vector...
    A(j+1:m,j) = v(2:m-j+1);
end
R = triu(A(1:n,1:n));
Q_fact = tril(A,-1);

%  Additional
R = double(R);
Q_fact = double(Q_fact);
Q = BackAccum(Q_fact);

end 

function Q = BackAccum(Q_fact,type)
% function Q = BackAccum(Q_fact,type)
% Explicit formation of an orthogonal matrix from its factored form
%   representation. Uses backward accumulation.
% Q_fact is mxn and Q_fact(j+1:m,j) is the essential part of the 
% jth Householder matrix Q_{j}.
% A call of the form BackAccum(Q_fact,'thin') produces a "thin", mxn Q
% with orthonormal columns, the first n columns of Q_{1}...Q_{n}.
% A call of the form BackAccum(Q_fact) produces an mxm Q
% that is the product of of Q_{1}...Q_{n}.
% GVL4: Section 5.2.2
[m,n] = size(Q_fact);
if nargin==2 && m>n && strcmp(type,'thin')
    Q = [];
    for j=n:-1:1
        v = [1;Q_fact(j+1:m,j)];
        beta = 2/(v'*v);
        k = m-j;
        Q = [1 zeros(1,n-j);zeros(m-j,1) Q];
        if norm(Q_fact(j+1:m,j))>0
            Q = Q - (beta*v)*(v'*Q);
        end
    end
else
    Q = eye(max(m-n,1),max(m-n,1));
    for j=min(n,m-1):-1:1
        v = [1;Q_fact(j+1:m,j)];
        beta = 2/(v'*v);
        k = m-j;
        Q = [1 zeros(1,k);zeros(k,1) Q];
        if norm(Q_fact(j+1:m,j))>0
            Q = Q - (beta*v)*(v'*Q);
        end
    end
end
end 

function [beta,x] = house(x)
%HOUSE Householder matrix related to x
% [beta,v] = HOUSE(x) retures beta and v which can be used 
% to construct P = I - beta * v * v', and 
% P * x = norm(x) * e1. 
% The Householder vector is scaled such that v(1) = 1.
%
% GOVA13-MC4 Algorithm 5.3.1
m = length(x); 
sigma = x(2:m)' * x(2:m);
if sigma == 0 && x(1) >= 0
	beta = 0;
elseif sigma == 0 && x(1) < 0
	beta = -2;
else
	mu = norm(x);
	if x(1) <= 0
		x(1) = x(1) - mu;
	else
		x(1) = -1 * sigma/(x(1) + mu);
	end
	beta = 2 * x(1) * x(1)/(x'*x);
	x = x/x(1);
end
end