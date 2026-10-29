function [A_t,Dh_t] = fun_gen_seq(a_ref,dh_ref,seqOptions) 

% This function generates a sequence of possible true values of unknown
% parameters. Points are generated using either a Sobol sequence or Latin
% Hybercube sampling (LHS). 

% Note that generated points are stored in diagonal matrices!

% Syntax
% ------------------------------------------------------------------------
% [A_t,Dh_t] = fun_gen_seq(a_ref,dh_ref,seqOptions)

% Input
% ------------------------------------------------------------------------
% a_ref      : reference true values, added mass           (Ndofx1 double)
% dh_ref     : reference true values, added damping        (Ndofx1 double)

% seqOptions 
% ------------------------------------------------------------------------
% Nset       : set dependency of coordinates (dof related)      (1x1 char)
% type       : set sampling sequence type                       (1x1 char)
% Nmcs       : number of quasi-monte carlo samples            (1x1 double)
% aBounds    : scale factors that defines bounds on a_t       (1x2 double)  
% dhBounds   : scale factors that defines bounds on dh_t      (1x2 double)

% Output
% ------------------------------------------------------------------------
% A_t        : Possible added mass matrices        (NdofxNdofxNmcs double)
% Dh_t       : Possible added damping matrices     (NdofxNdofxNmcs double)

% ------------------------------------------------------------------------

Ndof = length(a_ref) ;
Nmcs = seqOptions.Nmcs ;

% Specifing bounds on hyperspace
% ------------------------------------------------------------------------
if isfield(seqOptions,'aBounds')
    a_lb = seqOptions.aBounds(1) ;  a_ub = seqOptions.aBounds(2) ;

else 
    a_lb = 0.8 ;    a_ub = 1.2 ;

end

if isfield(seqOptions,'dhBounds')
    dh_lb = seqOptions.dhBounds(1) ;    dh_ub = seqOptions.dhBounds(2) ;

else
    dh_lb = 0.8 ;   dh_ub = 1.2 ;

end

% Generating sequence
% ------------------------------------------------------------------------

if Nmcs == 1 
    a_val = a_ref ;
    dh_val = dh_ref ;

    A_t = diag(a_val) ;
    Dh_t = diag(dh_val) ;

else

    switch seqOptions.Nset(1:3)

        case 'all'  % This makes sobol based on all unknown parameters

            switch seqOptions.type(1:3)

                case 'sob'      % Sobol sequence sampling 

                    p = sobolset(2*Ndof) ;                                                  
                    x = net(p,Nmcs) ;

                case 'lhs'      % Latin hypercube sampling 

                    x = lhsdesign(Nmcs,2*Ndof) ; 
            
            end

            a_min  =  a_lb*a_ref ;     a_max = a_ub*a_ref ; 
            dh_min = dh_lb*dh_ref ;   dh_max = dh_ub*dh_ref ; 
        
            a_val  = x(:,1)' .* (a_max  -  a_min) +  a_min ;
            dh_val = x(:,2)' .* (dh_max - dh_min) + dh_min ;

        case 'sep' % separate: This considers independency of some unknowns 

            switch seqOptions.type(1:3)

                case 'sob'      % Sobol sequence sampling 

                    p = sobolset(2) ;                                                  
                    x = net(p,Nmcs) ;

                case 'lhs'      % Latin hybercube sampling 

                    x = lhsdesign(Nmcs,2) ; 
            
            end

            a_min  =  a_lb*a_ref ;     a_max = a_ub*a_ref ; 
            dh_min = dh_lb*dh_ref ;   dh_max = dh_ub*dh_ref ; 
        
            a_val  = x(:,1)' .* (a_max  -  a_min) +  a_min ;
            dh_val = x(:,2)' .* (dh_max - dh_min) + dh_min ;
    end

    A_t = zeros(Ndof,Ndof,Nmcs) ;
    Dh_t = zeros(Ndof,Ndof,Nmcs) ;

    for i = 1:1:Nmcs

        A_t(:,:,i) = diag(a_val(:,i)) ;
        Dh_t(:,:,i) = diag(dh_val(:,i)) ;

    end
end

% Visualization
% % ------------------------------------------------------------------------
% size = 0.25 ;
% 
% figure(1) ; clf
% sgtitle(['Type of sampling: ' seqOptions.type ', ' seqOptions.Nset])
% subplot(1,3,1); plot(a_val(1,:),dh_val(1,:),'o',LineWidth=1.2);
% pbaspect([1 1 1]); grid on
% xlabel('$a_{11}$'); ylabel('$d_{h,11}$'); title('Surge')
% %xlim([(1-size)*a_ref(1) (1+size)*a_ref(1)])
% %ylim([(1-size)*dh_ref(1) (1+size)*dh_ref(1)])
% xline(a_lb*a_ref(1),'--') ; xline(a_ub*a_ref(1),'--')
% yline(dh_lb*dh_ref(1),'--') ; yline(dh_ub*dh_ref(1),'--')
% 
% subplot(1,3,2); plot(a_val(2,:),dh_val(2,:),'o',LineWidth=1.2); 
% pbaspect([1 1 1]); grid on
% xlabel('a'); ylabel('$d_h$'); title('Sway')
% %xlim([(1-size)*a_ref(2) (1+size)*a_ref(2)])
% %ylim([(1-size)*dh_ref(2) (1+size)*dh_ref(2)])
% 
% subplot(1,3,3); plot(a_val(3,:),dh_val(3,:),'o',LineWidth=1.2); 
% pbaspect([1 1 1]); grid on
% xlabel('a'); ylabel('$d_h$'); title('Yaw')
% %xlim([(1-size)*a_ref(3) (1+size)*a_ref(3)])
% %ylim([(1-size)*dh_ref(3) (1+size)*dh_ref(3)])
% % ------------------------------------------------------------------------

end