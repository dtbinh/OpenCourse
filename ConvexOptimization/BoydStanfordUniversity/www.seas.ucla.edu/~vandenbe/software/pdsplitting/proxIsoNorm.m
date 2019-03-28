function u = proxIsoNorm(mu,x)
% x is m rows by n columns by 2

nrms = sqrt(x(:,:,1).^2 + x(:,:,2).^2);
factors = max(0,1 - mu./nrms);  % Could be dividing by 0 here.
factors = cat(3,factors,factors);

u = factors.*x;



end