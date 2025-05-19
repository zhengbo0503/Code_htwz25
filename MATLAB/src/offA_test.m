function offA_test(A)

n = size(A,1);
ul = 2^(-23);

fnrm = fnorm(A);

% Orthogonalization methods

[Ql,~] = eig(single(A)); Qd = double(Ql);
[Qt_qr,~] = qr(Qd);
[Qt_mgs,~] = mgs(Qd);
[Qt_cholqr,~] = cholqr(Qd);
Qt_ns = ns(Qd);

% tridiagonalization methods
Qt_tridiag = tridiagonal_method(A);

% collect errors
off_qr = off_error(A,fnrm,Qt_qr);
off_mgs = off_error(A,fnrm,Qt_mgs);
off_cholqr = off_error(A,fnrm,Qt_cholqr);
off_ns = off_error(A,fnrm,Qt_ns);
off_tridiag = off_error(A,fnrm,Qt_tridiag);

offA = [off_qr off_mgs off_cholqr off_ns off_tridiag];

bound = sqrt(n) * 5 * ul;



if sum(offA > bound) == 0
    fprintf("All good!\n");
else
    fprintf("off_qr  off_mgs  off_cholqr  off_ns  off_tridiag   bound\n");
    disp([offA bound]);
    fprintf("No good!\n");
end

end

function error = off_error(A,fnrm,Qt)
A = mp(A,34);
Qt = mp(Qt,34);
error = off(Qt'* A * Qt)/fnrm;
end