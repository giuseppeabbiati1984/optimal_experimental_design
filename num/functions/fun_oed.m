function varargout = fun_oed(theta,w,D,m,optimality,option)

% This function evaluates P, which is the cost related to individual
% experimental designs

% Syntax
% ------------------------------------------------------------------------
% P = fun_cost_ad(theta,w,D,m,optimality,option)

% Input
% ------------------------------------------------------------------------
% w          : row vector of frequency steps                  (1xN double)
% M          : mass matrix                              (NdofxNdof double)
% theta      : calibration parameters                    (2x1xNmcs double)
% design     : design variable                         (NvarxNtest double)         
% F2         : LF loads matrix                     (NdofxN complex double) 
% optimality : 'A' or 'D'                                       (1x1 char)

% option                                                      (1x1 struct)
% ------------------------------------------------------------------------
% type       : 'debug', 'oed' ('oed' by default)                (1x1 char)
% noise      : Eps                                      ()
            
% geometry                                                    (1x1 struct)
% ------------------------------------------------------------------------
% sys        : 'diagonal' or 'nondiagonal'                      (1x1 char)
% type       : 'triangle' or 'rectangle'                        (1x1 char)
% dxe        : length of mooring line in x                    (1x1 double)
% dye        : length of mooring line in y                    (1x1 double)
% dxf        : floater dimension                              (1x1 double)
% dyf        : floater dimension                              (1x1 double)

% Output
% ------------------------------------------------------------------------
% P          : objective function                          (1xNmcs double)
% ------------------------------------------------------------------------

arguments    
    theta
    w
    D
    m
    optimality
    option.type = 'oed'
end

% ------------------------------------------------------------------------
%Nmcs  = numel(A_t(1,1,:)) ;
Nmcs = length(theta(1,1,:)) ;

P = zeros(Nmcs,1) ;

% ------------------------------------------------------------------------

for i = 1:1:Nmcs 
                               
        %H = fun_hess_cpsd(theta,w,D,m) ;
        H = fun_hess_cpsd_v2(theta(:,:,i),w,D,m) ;

    switch optimality 

        case 'A' % A-optimality
        
            P(i) = trace((H)^-1) ;

        case 'D' % D-optimality
            
            %P(i) = log(det((HessQ)^-1)) ; % (Lopez-Fidalgo, 2023) ch.2.4.3
            P(i) = -log(det((H))) ;    % (Lopez-Fidalgo, 2023) ch.2.4.3
            %P(i) = det((HessQ)^-1) ;

        case 'E' % E-optimality

            P(i) = min(eig(H))^-1 ; 

    end 

end

switch option.type

    case 'debug'

        varargout{1} = mean(P) ;
        varargout{2} = P ;
        %varargout{3} = X ;

    case 'oed'

        varargout{1} = mean(P) ;

end

end