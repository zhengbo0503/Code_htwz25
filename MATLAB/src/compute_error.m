function err = compute_error(D_ref, D_comp, varargin)
if nargin > 2
    prec = varargin{1};
else
    prec = 71;
end
D_ref = mp(sort(D_ref), prec);
D_comp = mp(sort(D_comp), prec);
[n_ref,~] = size(D_ref);
[n_comp,~] = size(D_ref);
if (n_ref ~= n_comp)
    D_ref = D_ref';
end
err = abs((D_ref - D_comp)./D_ref);
end
