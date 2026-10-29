function [x,y,fsave] = fun_gen_resp_v2(w,f1,f2,pars,noise)

% System
% ------------------------------------------------------------------------
  
w_0 = 0.5 ;

A_ref = pars.A_ref ;    Dh_ref = pars.Dh_ref ;
A_dbl = pars.A_dbl ;    D_dbl = pars.Ddbl;

M = pars.M ; K = pars.K ;   D = pars.D ;

N = length(K) ;     r = length(w) ;

% Load preparation
% ------------------------------------------------------------------------

F_LF = zeros(N,r) ;

for i = 1:1:N

    f2e = awgn(f2(1,:),noise.loadSNR,'measured') ;

    if i == 1

        fsave = f2e ;
    end

    F_LF(i,:) = fftshift(fft(f2e)) ;

end

    F_WF = fftshift(fft(f1(1,:))) ;
    
    %F_WF =  0 ;

% response
% ------------------------------------------------------------------------

H_LF = zeros(N,r) ;
H_WF = zeros(N,r) ;

for i = 1:1:N

    H_LF(i,:) = (-w.^2*(M + A_ref) + ...
        1i*w*(Dh_ref + D(i)) + K(i)).^-1 ;   

    H_WF(i,:) = (-w.^2.*(M + A_dbl) + ...
        1i*w.*(D_dbl + D(i)) + K(i)).^-1 ;

end

list = abs(w)>w_0 ;

H = list.*H_WF + (~list).*H_LF ; H(:,1) = 0 ;

X = H.*(F_LF + F_WF);

x = zeros(size(X)) ;

for i = 1:1:N

    %X_shifted = [X(i,end/2 + 1:end) ny X(i,2:end/2)] ; 
    %x(i,:) = ifft(X_shifted)*r ;

    x(i,:) = ifft(ifftshift(X(i,:))) ;
end

% Noise on response
% ------------------------------------------------------------------------
%rng(1.0)

switch noise.opt 

    case '1'

        % (1) same noise level on all responses
        % ----------------------------------------------------------------
        
        y = awgn(x,noise.SNR,"measured") ;

    case '2'

        % (2) constant SNR based on individual responses
        % ----------------------------------------------------------------
        y = zeros(size(x)) ;

        for i = 1:1:N
        
                y(i,:) = awgn(x(i,:),noise.SNR,'measured') ;
        
        end 

end

