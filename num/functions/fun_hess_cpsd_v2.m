function H = fun_hess_cpsd_v2(theta,w,D,m)

% D:        mooring design parameters (2 x n matrix)
% theta:    calibration parameters

n = size(D,2);
r = length(w);

p = 1;
q = 2;

% System terms
Kp_s = -w.^2*(theta(1) + m) - 1i*w.*(theta(2) + D(2,p)) + D(1,p);
Kq   = -w.^2*(theta(1) + m) + 1i*w.*(theta(2) + D(2,q)) + D(1,q);
Kn   = -w.^2*(theta(1) + m) + 1i*w.*(theta(2) + D(2,n)) + D(1,n);

Hp = 1./Kp_s;
Hq = 1./Kq;
Hn = 1./Kn;

dZ1 = zeros(n-1,r);
dZ2 = zeros(n-1,r);
Wgt = zeros(n-1,r);

for i = 1:n-1

    Ki_s = -w.^2*(theta(1) + m) - 1i*w.*(theta(2) + D(2,i)) + D(1,i);
    Hi   = 1./Ki_s;

    % Jacobian terms
    dZ1(i,:) = (-w.^2.*(Kp_s + Kq).*Ki_s.*Kn + ...
                 w.^2.*(Ki_s + Kn).*Kp_s.*Kq) ./ (Ki_s.*Kn).^2;

    dZ2(i,:) = (1i*w.*(Kp_s - Kq).*Ki_s.*Kn - ...
                 1i*w.*(Ki_s - Kn).*Kp_s.*Kq) ./ (Ki_s.*Kn).^2;

    % Welch-aware diagonal weights (up to proportionality)
    Wgt(i,:) = (abs((Hi .* Hn) ./ (Hp .* Hq)).^2).^-1;

end

% Vectorize
v1 = dZ1(:);
v2 = dZ2(:);
wvec = Wgt(:);

% Weighted Hessian
H11 = real(sum(wvec .* conj(v1) .* v1));
H12 = real(sum(wvec .* conj(v1) .* v2));
H21 = real(sum(wvec .* conj(v2) .* v1));
H22 = real(sum(wvec .* conj(v2) .* v2));

H = [H11, H12; H21, H22];
