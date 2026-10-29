clc ; clear ; close all

% ========================================================================
%                               INPUTS
% ========================================================================

mydir = fileparts(pwd) ;
addpath([mydir '\num\functions'])

% Floater parameters 
% ------------------------------------------------------------------------

m = 1.4176e4  ;                     % mass                    [t]

% km = [53 90 53]  ;                   % mooring stiffness (3)   [kN/m]
% dm = [218 0 0]  ;                    % mooring damping (3)     [kNs/m]

% km = [53 90 53]  ;                   % mooring stiffness (3)   [kN/m]
% dm = [0  0 218]  ;                   % mooring damping (3)     [kNs/m]

% km = [90 53 90]  ;                   % mooring stiffness (3)   [kN/m]
% dm = [284 218 0]  ;                  % mooring damping (3)     [kNs/m]

% km = [90 53 90]  ;                   % mooring stiffness (3)   [kN/m]
% dm = [0 218 284]  ;                  % mooring damping (3)     [kNs/m]

km = [90 90 90]  ;                     % mooring stiffness (3)   [kN/m]
dm = [284 0   0]  ;                    % mooring damping (3)     [kNs/m]

%km = [150 40 150]  ;                  % mooring stiffness (3)   [kN/m]
%dm = [ 0   0  0]  ;                   % mooring damping (3)     [kNs/m]

% km = [150 40 40]  ;                   % mooring stiffness (3)   [kN/m]
% dm = [ 0   0  0]  ;                   % mooring damping (3)     [kNs/m]

n = length(km) ;

D_trial = zeros(2,n) ;

    D_trial(1,:) = km ;
    D_trial(2,:) = dm ;

theta = [1.0e4 , 700] ; 

dw = 0.001 ;
w = 0:dw:0.25 ;

% Check hessian function
% ------------------------------------------------------------------------
H = fun_hess_cpsd_v2(theta,w,D_trial,m) ;

% Check Z-functions
% ------------------------------------------------------------------------
Kp_s = -w.^2*(theta(1) + m) - 1i*w.*(theta(2)*0.5 + D_trial(2,1)) ...
    + D_trial(1,1);
Kq   = -w.^2*(theta(1) + m) + 1i*w.*(theta(2)*0.5 + D_trial(2,2)) ...
    + D_trial(1,2);
Kn   = -w.^2*(theta(1) + m) + 1i*w.*(theta(2)*0.5 + D_trial(2,n)) ...
    + D_trial(1,n);

Z = zeros(2,length(w)) ;

for i = 1:n-1

    Ki_s = -w.^2*(theta(1) + m) - 1i*w.*(theta(2)*0.5 + D_trial(2,i)) ...
        + D_trial(1,i);
    
       Z(i,:) = Kp_s.*Kq./(Ki_s.*Kn) ;

end

figure(1) ; clf ; 
subplot(2,1,1) ; hold all ; box on
plot(w,real(Z(1,:)))
plot(w,real(Z(2,:)))
ylim([-1 4])

subplot(2,1,2) ; hold all ; box on
plot(w,imag(Z(1,:)))
plot(w,imag(Z(2,:)))
ylim([-4 1])

figure(2) ; clf 
subplot(2,1,1) ; hold all ; box on
plot(w,abs(Z(1,:)))
plot(w,abs(Z(2,:)))

subplot(2,1,2) ; hold all ; box on
plot(w,phase(Z(1,:)))
plot(w,phase(Z(2,:)))

%%
% ========================================================================
%                     OPTIMAL EXPERIMENTAL DESIGN
% ========================================================================

% Initial guess for optimal control variables 
D_init = 1 + randn(2,n)/50 ;  

% Lower and upper bound for control variables                               
D_lb = ones(size(D_init)) ; 
D_lb(1,:) = 0.4*D_lb(1,:) ; D_lb(2,:) = 0*D_lb(1,:) ;
D_ub = ones(size(D_init)) ; 
D_ub(1,:) = 1.5*D_ub(1,:) ; D_ub(2,:) = 3.0*D_ub(2,:) ;

% Objective function                                                        
Pcost = @(D) fun_oed(theta_in,w,D.*ref,m,optimality) ;

opt_algorithm = 'ga'

tic 
switch opt_algorithm

    case 'gd'                       % Gradient descent (gd)
    % --------------------------------------------------------------------

        options = optimoptions('fmincon','display','iter') ; 
        D_opt   = fmincon(Pcost,D_init,[],[],[],[],...
                  D_lb,D_ub,[],options) ; 

    case 'gs'                       % Global search (gs)
    % --------------------------------------------------------------------

        rng default                 % For reproducibility

        options = optimoptions(@fmincon,'Algorithm','sqp',...
            'display','iter');

        options.OptimalityTolerance    = 1e-5 ;
        options.FunctionTolerance      = 1e-5 ; 
        options.StepTolerance          = 1e-5 ;
        options.ConstraintTolerance    = 1e-5 ;
        options.MaxIterations          = 10 ; 
        options.MaxFunctionEvaluations = 100 ;

        problem = createOptimProblem('fmincon','objective',Pcost,'x0',...
            D_init,'lb',D_lb,'ub',D_ub,'options',options);
        gs = GlobalSearch;
        [x,f] = run(gs,problem)

end
toc

% Note that km = [40 150 40] and km = [150 40 150] result in the same
% cost.

%%
% ========================================================================
%                      INVESTIGATE DESIGN SPACE
% ========================================================================

% Define reference design
%D_ref = ones(2,3);  % or any nominal values you want as the fixed back
D_ref = D_trial/100 ;

par = 'stiffness' ;   % 'stiffness', 'damping' , 'mixed'

% Choose two indices to vary
 j1 = 1;  % first parameter
 j2 = 2;  % second parameter


switch par

    case 'stiffness'
           i1 = 1;
           i2 = 1;

           d1 = linspace(1*0.1, 10*1, 60);
           d2 = linspace(1*0.1, 10*1, 60);

    case 'damping'
            i1 = 2;
            i2 = 2;

           d1 = linspace(0*1, 3.0*1, 60);
           d2 = linspace(0*1, 3.0*1, 60);

    case 'mixed'
            i1 = 1;
            i2 = 2;

           d1 = linspace(0.4*1, 1.5*1, 60);
           d2 = linspace(0*1, 3.0*1, 60);
end

[D1,D2] = meshgrid(d1,d2);

% Evaluate cost on the grid
Z = zeros(size(D1));
for k = 1:numel(D1)
    D = D_ref;
    D(i1,j1) = D1(k);
    D(i2,j2) = D2(k);
    Z(k) = Pcost(D);
end

% Plot contour
% ------------------------------------------------------------------------
figure(1) ; clf 
contourf(D1,D2,Z,30); colorbar
xlabel(sprintf('D(%d,%d)', i1,j1))
ylabel(sprintf('D(%d,%d)', i2,j2))
title('Cost contours for two design parameters')
hold on
plot(D_ref(i1,j1), D_ref(i2,j2), 'rx', 'LineWidth', 2) % mark ref point
hold off

% Surface plot
% ------------------------------------------------------------------------
figure(2) ; clf 
surf(D1,D2,Z)
shading interp
colorbar
xlabel(sprintf('D(%d,%d)', i1,j1))
ylabel(sprintf('D(%d,%d)', i2,j2))
zlabel('Cost')
%title('Surface plot of cost function')
hold on
% plot3(D_ref(i1,j1), D_ref(i2,j2), Pcost(D_ref), 'rx',...
%     'MarkerSize', 10, 'LineWidth', 2) % reference point
hold off
