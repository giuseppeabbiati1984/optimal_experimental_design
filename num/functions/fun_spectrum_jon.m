function S_jon = fun_spectrum_jon( Tp, Hs, w, std)
% D. Stamenov (Aug 2023). JONSWAP spectrum for generating hydrodynamic loading

% Calculates the JONSWAP wave spectrum S
% Tp: Wave peak period [s]
% Hs: Significant wave height [m]
% w: Vector of wave frequencies (x-axis) [rad/s]
% g: Gravitational acceleration [m/s^2]
% std: not used but included for parallelism to the bell input
% OUTPUT
% S: Wave spectrum [m^2*s/rad] (length(w) x 1)


%###############COMPUTE JONSWAP############################################

% Formulation from: Zhengshun Cheng et al. (2018)

    wp = 2*pi/Tp ; % Wave peak frequency
    g = 9.81;         % gravitational acceleration

    % Calculate Gamma
    if Tp/sqrt(Hs) <= 3.6
        gamma = 5 ;
    elseif 3.6 < Tp/sqrt(Hs) && Tp/sqrt(Hs) < 5
        gamma = exp(5.75 - 1.15*Tp/sqrt(Hs)) ;
    else
        gamma = 1 ;
    end

    % Spectral parameter 1
    alpha = 5.061 * Hs^2/Tp^4 * (1 - 0.287 * log(gamma)) ;

    S_jon = zeros(length(w),1) ;

    for i=1:1:length(w)
        % Spectral parameter 2
        if w(i)<wp
            sigma = 0.07 ;
        else
            sigma = 0.09 ;
        end
        % Unidirectional wave spectrum (JONSWAP)
        S_jon(i) = (alpha*g^2)/w(i)^5 * exp(-1.25*(wp/w(i))^4) * gamma^(exp(-(w(i)-wp)^2/(2*sigma^2*wp^2))) ;
    end

end