clc ; clear ; close all

% ========================================================================
%                             INPUT AND SETUP
% ========================================================================

inputs_08_01_2026       % AsRun-ED


color1 = '#0072BD' ;        % Blue
color2 = '#A2142F' ;        % Red
color3 = '#EDB120' ;        % Yellow

color12 = '#7E2F8E' ;       % purple
color13 = '#77AC30' ;       % green
color23 = '#D95319' ;       % orange

%%
% ========================================================================
%                  SIMULATION OF RESPONSE MEASUREMENTS
% ========================================================================

% System and response
% ------------------------------------------------------------------------

rng(1.0)

pars.A_ref = A_ref ;    pars.Dh_ref = Dh_ref ;
pars.A_dbl = A_dbl ;    pars.Ddbl = D_dbl ;

pars.M = M ; pars.K = K ;  pars.D = D ;

[x,y_raw] = fun_gen_resp_v2(w,f1,f2,pars,noise) ;

%%
figure(1) ; clf 
set(gcf,'Position',[600 100 1400 1000])
subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
plot(t,y_raw(1,:),'-k')
plot(t,x(1,:),'-','Color',color1)
legend('Noise polluted response, $y_1$','Actual response, $x_1$')
xlabel('Time, t [s]')
ylabel('Displacement [m]')

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on 
plot(t,y_raw(2,:),'-k')
plot(t,x(2,:),'-','Color',color2)
legend('Noise polluted response, $y_2$','Actual response, $x_2$')
xlabel('Time, t [s]')
ylabel('Displacement [m]')

subplot(3,1,3,'fontsize',fontS) ; hold all ; box on 
plot(t,y_raw(3,:),'-k')
plot(t,x(3,:),'-','Color',color3)
legend('Noise polluted response, $y_3$','Actual response, $x_3$')
xlabel('Time, t [s]')
ylabel('Displacement [m]')

figure(2) ; clf 
set(gcf,'Position',[600 100 1400 1000])

subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
histogram(y_raw(1,:),'FaceColor',color1)
%histogram(x(1,:),'FaceColor',color1)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-10,40])

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on 
histogram(y_raw(2,:),'FaceColor',color2)
%histogram(x(2,:),'FaceColor',color2)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-10,40])

subplot(3,1,3,'fontsize',fontS) ; hold all ; box on 
histogram(y_raw(3,:),'FaceColor',color3)
%histogram(x(3,:),'FaceColor',color3)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-10,40])

%%
% ========================================================================
%                        CHECK SPECTRAL ESTIMATES
% ========================================================================

%y0 = x - mean(x,2) ;
y0 = y_raw - mean(y_raw,2) ;

fs = 1/dt ;
Gmat = cell(N,N) ;

%nfft = 2^11 ;
nfft = 2000 ;
Nw = nfft ;
%Nw = [] ;

for i = 1:1:N

    for j = 1:1:N

        %[G,f] = cpsd(y0(i,:),y0(j,:),[],[],nfft,fs) ; 
        [G,f] = cpsd(y0(i,:),y0(j,:),Nw,[],nfft,fs) ;

        Gmat{i,j} = G'/(2*pi) ;

    end
end

w2 = (2*pi*f)' ;

figure(3) ; clf ; hold all ; box on 
set(gcf,'Position',[800 300 1200 400])

x_start = 0;
x_end = 0.15;
y_min = 1e-2;
y_max = 1e4;

% Fade parameters
n_fade = 20;  % number of fading steps
fade_color = [0.8 0.9 1];  % light blue

for i = 1:n_fade
    % Linearly interpolate x position
    x1 = x_start + (i-1)*(x_end - x_start)/n_fade;
    x2 = x_start + i*(x_end - x_start)/n_fade;
    
    % Decrease alpha from left to right
    alpha = 1 - (i-1)/n_fade;
    

    fill([x1 x2 x2 x1], [y_min y_min y_max y_max], ...
         fade_color, ...
         'EdgeColor', 'none', ...
         'FaceAlpha', alpha, ...
         'HandleVisibility', 'off');  % not shown in legend
end

%subplot(2,1,1,'FontSize',16) ; hold all ; box on 
% plot(w2,abs(Gmat{1,2}).^2./(Gmat{1,1}.*Gmat{2,2}),'Color','#7E2F8E',...
%     'LineWidth',1.5) 
% plot(w2,abs(Gmat{1,3}).^2./(Gmat{1,1}.*Gmat{3,3}),'Color',"#77AC30",...
%     'LineWidth',1.5) 
% plot(w2,abs(Gmat{2,3}).^2./(Gmat{2,2}.*Gmat{3,3}),'Color',"#D95319",...
%     'LineWidth',1.5) 
% ylim([0,1.5])
% %xlim([0,0.1])
% xlim([0,1])
% %xlim([0,0.25])
% %xlim([0,2.0])
% xlabel('Angular frequency, $\omega$ [rad/s]')
% ylabel('Coherence, $\gamma^2_{ij}$')
% legend('$\gamma^2_{12}$','$\gamma^2_{13}$','$\gamma^2_{23}$')
% title('\textbf{Numerical}')

% subplot(2,1,2,'yscale','linear','FontSize',16) ; hold all ; box on 
plot(w2,Gmat{1,1},...
    'LineWidth',1.5)
plot(w2,Gmat{2,2},'Color','#A2142F',...
    'LineWidth',1.5)
plot(w2,Gmat{3,3},...
    'LineWidth',1.5)
%xlim([0,0.1])
xlim([0,1])
%xlim([0,0.25])
legend('$\hat{G}_{\tilde{x}_1\tilde{x}_1}$',...
    '$\hat{G}_{\tilde{x}_2\tilde{x}_2}$','$\hat{G}_{\tilde{x}_3\tilde{x}_3}$')
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('PSD, $\hat{G}_{\tilde{x}_i\tilde{x}_i}(\omega)$ [m$^2$s/rad]')
yscale('log')
ylim([1e-2,1e4])

fontsize(fontS,'points')

% figure(5) ; clf ; hold all ; box on 
% set(gcf,'Position',[800 300 1200 400])
% axis equal
% plot(real(Gmat{1,2}./Gmat{2,3}),real(Gmat{1,3}./Gmat{2,3}),'*')
%plot(real(Gmat{2,3}./Gmat{1,3}),real(Gmat{2,3}./Gmat{1,2}),'*')
%plot(imag(Gmat{1,2}./Gmat{2,3}),imag(Gmat{1,3}./Gmat{2,3}),'*')

%%
% ========================================================================
%            CUT-OFF SENSITIVITY (VARIANCE-BASED ESTIMATION)
% ========================================================================

w_cut = 0.01:0.01:0.25 ; 

lb = -4.0*[A_ref;Dh_ref] ;
ub = 10.0*[A_ref;Dh_ref] ;

rng(1.0)
theta_init = lb + (ub - lb) .* rand(2,1) ;

theta_opt_mcs2 = zeros(2,length(w_cut)) ;

for i = 1:1:length(w_cut)

    list = abs(w)<w_cut(i) ;

    Y = zeros(size(x)) ;   

    for j = 1:1:N

        Y(j,:) = fftshift(fft(y_raw(j,:))).*list ;
        Yin{j} = Y(j,:) ;

    end
    fun_freq = @(theta) fun_obj_VAR(theta,M,K,D,w,Yin) ; 

    theta_opt_mcs2(:,i) = fmincon(fun_freq,theta_init,[],[],[],[],lb,ub,[]) ;
end

%%
% ========================================================================
%              CUT-OFF SENSITIVITY (CPSD-BASED ESTIMATION)
% ========================================================================

theta_opt_mcs3 = zeros(2,length(w_cut)) ;

lb = 0.0*[A_ref;Dh_ref] ;
ub = 4.0*[A_ref;Dh_ref] ;

rng(1.0)
theta_init = lb + (ub - lb) .* rand(2,1) ;

Qoptions.idx = [1 2] ;
Qoptions.diff = 'NLSE' ;

for i = 1:1:length(w_cut)

    Qoptions.wcut = w_cut(i) ;

    fun_cpsd = @(theta) fun_obj_CPSD(theta,M,K,D,Gmat,w2,Qoptions) ;

    theta_opt_mcs3(:,i) = fmincon(fun_cpsd,theta_init,[],[],[],[],lb,ub,[]) ;

    disp(num2str(i))

end

%%
figure(4) ; clf ; hold all ; box on %; grid on
% ------------------------------------------------------------------------
%set(gcf,'Position',[1600 300 1200 400])
set(gcf,'Position',[1600 300 800 400])
plot(w_cut,theta_opt_mcs2(1,:),'x',...
    'MarkerEdgeColor',"#D95319",'LineWidth',1.5,'MarkerSize',10)

plot(w_cut,theta_opt_mcs3(1,:),'o',...
    'MarkerEdgeColor',"#0072BD",'LineWidth',1.5)

ylim([0,2e4])
fontsize(fontS,'points')
xlabel('$w_{cut}$ [rad/s]')
ylabel('Added mass estimate, $\hat{a}$ [t]')
legend('S2021 Estimator','C2025 Estimator')
%legend('Estimates')

figure(5) ; clf ; hold all ; box on %; grid on
% ------------------------------------------------------------------------
set(gcf,'Position',[1600 300 800 400])
plot(w_cut,theta_opt_mcs2(2,:),'x',...
    'MarkerEdgeColor',"#D95319",'LineWidth',1.5,'MarkerSize',10)

plot(w_cut,theta_opt_mcs3(2,:),'o',...
    'MarkerEdgeColor',"#0072BD",'LineWidth',1.5)

ylim([0,1e3])
fontsize(fontS,'points')
xlabel('$w_{cut}$ [rad/s]')
ylabel('Damping estimate, $\hat{d}_h$ [kNs/m]')
legend('S2021 Estimator','C2025 Estimator')

%%
% ========================================================================
%           MONTE CARLO SAMPLING (VARIANCE-BASED ESTIMATION 0.1)
% ========================================================================

w_cut = 0.1 ;

N_samples = 1e3 ;
theta_opt_mcs1 = zeros(2,N_samples) ;

lb = -4.0*[A_ref;Dh_ref] ;
ub = 4.0*[A_ref;Dh_ref] ;

rng(1.0)
theta_init = lb + (ub - lb) .* rand(2,1) ;
list = abs(w)<w_cut ;

for i = 1:1:N_samples

    [x,y_raw] = fun_gen_resp_v2(w,f1,f2,pars,noise) ;

    Y = zeros(size(x)) ;   

    for j = 1:1:N

        Y(j,:) = fftshift(fft(y_raw(j,:))).*list ;                          % est. with output noise
        Yin{j} = Y(j,:) ;

    end

    fun_freq = @(theta) fun_obj_VAR(theta,M,K,D,w,Yin) ; 

    theta_opt_mcs1(:,i) = fmincon(fun_freq,theta_init,[],[],[],[],lb,ub,[]) ;

    disp(num2str(i))

end

%%
% ========================================================================
%         MONTE CARLO SAMPLING (VARIANCE-BASED ESTIMATION 0.2)
% ========================================================================

w_cut = 0.2 ;

N_samples = 1e3 ;
theta_opt_mcs2 = zeros(2,N_samples) ;

lb = -4.0*[A_ref;Dh_ref] ;
ub = 4.0*[A_ref;Dh_ref] ;

rng(1.0)
theta_init = lb + (ub - lb) .* rand(2,1) ;
list = abs(w)<w_cut ;

for i = 1:1:N_samples
    
    [x,y_raw] = fun_gen_resp_v2(w,f1,f2,pars,noise) ;

    Y = zeros(size(x)) ;   

for j = 1:1:N

        Y(j,:) = fftshift(fft(y_raw(j,:))).*list ;                          % est. with output noise
        Yin{j} = Y(j,:) ;

end

    fun_freq = @(theta) fun_obj_VAR(theta,M,K,D,w,Yin) ; 

    theta_opt_mcs2(:,i) = fmincon(fun_freq,theta_init,[],[],[],[],lb,ub,[]) ;

    disp(num2str(i))

end

%%
% ========================================================================
%               MONTE CARLO SAMPLING (CPSD-BASED ESTIMATION)
% ========================================================================

N_samples = 1e3  ;
theta_opt_mcs3 = zeros(2,N_samples) ;

lb = 0.0*[A_ref;Dh_ref] ;
ub = 4.0*[A_ref;Dh_ref] ;

rng(1.0)
theta_init = lb + (ub - lb) .* rand(2,1) ;

Qoptions.idx = [1 2] ;  % {p,q}={1,2}
Qoptions.diff = 'NLSE' ;

for i = 1:1:N_samples
    
    [x,y_raw] = fun_gen_resp_v2(w,f1,f2,pars,noise) ;

     y0 = y_raw - mean(y_raw,2) ;
   
        fs = 1/dt ;
        Gmat = cell(N,N) ;

        for k = 1:1:N
        
            for j = 1:1:N
        
                [G,f] = cpsd(y0(k,:),y0(j,:),Nw,[],nfft,fs) ; 
            
                Gmat{k,j} = G'/(2*pi) ;
        
            end
        end

        w2 = (2*pi*f)' ;

        fun_cpsd = @(theta) fun_obj_CPSD(theta,M,K,D,Gmat,w2,Qoptions) ;
        
    theta_opt_mcs3(:,i) = fmincon(fun_cpsd,theta_init,[],[],[],[],lb,ub,[]) ;
    
    disp(num2str(i))

    Gsave{i} = Gmat ;

end

%%
% ------------------------------------------------------------------------
% Plots
% ------------------------------------------------------------------------

[l1,h1] = ksdensity(theta_opt_mcs1(1,:)) ;
[l2,h2] = ksdensity(theta_opt_mcs2(1,:)) ;
[l3,h3] = ksdensity(theta_opt_mcs3(1,:)) ;

figure(6) ; clf ; hold all ; box on ; grid on
% ------------------------------------------------------------------------
xline(A_ref,'k','LineWidth',2)
plot(h1,l1,'LineWidth',1.5,'Color',color1)
plot(h2,l2,'--','LineWidth',1.5,'Color',color2)
plot(h3,l3,':','LineWidth',1.5,'Color',color3)

xlabel('LF added mass $\hat{a}$ [t]','fontSize',16)
%xlim([0.5*A_ref 1.5*A_ref])
xlim([0,14000])
%xlim([4500 ; 11000])
legend('True value','Variance-based estimation, $\omega_0=0.1$',...
    'Variance-based estimation, $\omega_0=0.2$',...
        'CPSD-based estimation','Location','northwest')
fontsize(fontS,'points')
set(gcf,'Position',[800 300 1200 800])
title(['$T=$ ' num2str(Tmax) ' s' '   and   ',...
    '$n=$ ' num2str(N)])

[l1,h1] = ksdensity(theta_opt_mcs1(2,:)) ;
[l2,h2] = ksdensity(theta_opt_mcs2(2,:)) ;
[l3,h3] = ksdensity(theta_opt_mcs3(2,:)) ;


figure(7) ; clf ; hold all ; box on ; grid on
% ------------------------------------------------------------------------
xline(Dh_ref,'k','LineWidth',2)
plot(h1,l1,'LineWidth',1.5,'Color',color1)
plot(h2,l2,'--','LineWidth',1.5,'Color',color2)
plot(h3,l3,':','LineWidth',1.5,'Color',color3)
xlabel('Wave drift damping $\hat{d}_h$ [kNs/m]','fontSize',16)
xlim([0,1000])
%xlim([0.5*Dh_ref 1.5*Dh_ref])
%xlim([340 ; 470])
legend('True value','Variance-based estimation, $\omega_0=0.1$',...
    'Variance-based estimation, $\omega_0=0.2$',...
        'CPSD-based estimation','Location','northwest')
fontsize(fontS,'points')
set(gcf,'Position',[800 300 1200 800])
title(['$T=$ ' num2str(Tmax) ' s' '   and   ',...
    '$n=$ ' num2str(N)])

