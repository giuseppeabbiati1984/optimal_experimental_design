clc ; clear ; close all 

% ========================================================================
%                               INPUTS
% ========================================================================

input_12_09_2025

%%
% ========================================================================
%                            ALIGNING SIGNALS
% ========================================================================

s_idx = 2e4 ;                 % idx for cutting off initial part

[t_resampled,y_resampled] = fun_align_resp(y_raw,fs_raw,s_idx) ;

dt_resampled = t_resampled(2) - t_resampled(1) ;
fs_resampled = 1/dt_resampled ;

figure(2) ; clf 
set(gcf,'Position',[800 300 1600 1000])
subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
plot(t_raw{1},y_raw{1},'Color',[.7 .7 .7])
plot(t_resampled,y_resampled(1,:),'Color',color1)
ylabel('$\tilde{x}_1$ [m]')
%xlabel('Time, $t$ [s]')

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on 
plot(t_raw{2},y_raw{2},'Color',[.7 .7 .7])
plot(t_resampled,y_resampled(2,:),'Color',color2)
ylabel('$\tilde{x}_2$ [m]')
%xlabel('Time, $t$ [s]')

subplot(3,1,3,'fontsize',fontS) ; hold all ; box on 
plot(t_raw{3},y_raw{3},'Color',[.7 .7 .7])
plot(t_resampled,y_resampled(3,:),'Color',color3)
ylabel('$\tilde{x}_3$ [m]')
xlabel('Time, $t$ [s]')

figure(3) ; clf 
set(gcf,'Position',[600 100 1600 1000])
subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
histogram(y_resampled(1,:),'FaceColor',color1)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-40,10])

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on 
histogram(y_resampled(2,:),'FaceColor',color2)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-40,10])

subplot(3,1,3,'fontsize',fontS) ; hold all ; box on 
histogram(y_resampled(3,:),'FaceColor',color3)
ylabel('Number of samples')
xlabel('Displacement [m]')
xlim([-40,10])

mean_y = mean(y_resampled,2) ;
y0_resampled = y_resampled - mean_y ;
[R13,lags] = xcorr(y0_resampled(1,:),y0_resampled(3,:)) ;
[R12, ~] = xcorr(y0_resampled(1,:),y0_resampled(2,:)) ;
[R23, ~] = xcorr(y0_resampled(2,:),y0_resampled(3,:)) ;

idx13 = find(R13==max(R13)) ;   lag013 = lags(idx13)/fs_resampled ;
idx12 = find(R12==max(R12)) ;   lag012 = lags(idx12)/fs_resampled ;
idx23 = find(R23==max(R23)) ;   lag023 = lags(idx23)/fs_resampled ;

figure(4) ; clf
set(gcf,'Position',[800 300 1600 1000])
subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
plot(lags,R12,'Color',color12)
xlabel('Time lag, $\tau$ [s]')
ylabel('$\hat{R}_{y_1y_2}(\tau)$')

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on
plot(lags,R13,'Color',color13)
xlabel('Time lag, $\tau$ [s]')
ylabel('$\hat{R}_{y_1y_3}(\tau)$')

subplot(3,1,3,'fontsize',fontS) ; hold all ; box on
plot(lags,R23,'Color',color23)
xlabel('Time lag, $\tau$ [s]')
ylabel('$\hat{R}_{y_2y_3}(\tau)$')

%%
% ========================================================================
%                        DOWNSAMPLED SIGNALS
% ========================================================================

DSopt.fac = 400 ;                           % LF only
%DSopt.fac = 100 ;                          % LF and WF

DSopt.filt_type = 'fir' ;                   %  'butter' or 'fir'

[t,y] = fun_down_sample(t_resampled,y_resampled,fs_resampled,DSopt) ;

y0 = y - mean(y,2) ;
dt = t(2) - t(1) ;
Tmax = (t(end)-t(1)) + dt ;
r = length(t) ;
w = linspace(-r/2,r/2-1,r)*2*pi/Tmax ;

figure(6) ; clf 
set(gcf,'Position',[800 300 1200 1000])
subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
plot(t_raw{1},y_raw{1},'Color',[.7 .7 .7])
plot(t,y(1,:),'Color',color1,'LineWidth',1.5)
%xlabel('Time, $t$ [s]')
%ylabel('$y_1(t)$ [m]')
ylabel('$\tilde{x}_1(t)$ [m]')

subplot(3,1,2,'fontsize',fontS) ; hold all ; box on 
plot(t_raw{2},y_raw{2},'Color',[.7 .7 .7])
plot(t,y(2,:),'Color',color2,'lineWidth',1.5)
%xlabel('Time, $t$ [s]')
%ylabel('$y_2(t)$ [m]')
ylabel('$\tilde{x}_2(t)$ [m]')

subplot(3,1,3,'fontSize',fontS) ; hold all ; box on 
plot(t_raw{3},y_raw{3},'Color',[.7 .7 .7])
plot(t,y(3,:),'Color',color3,'LineWidth',1.5)
xlabel('Time, $t$ [s]')
%ylabel('$y_3(t)$ [m]')
ylabel('$\tilde{x}_3(t)$ [m]')

%%
% ========================================================================
%                   MEASUREMENT CROSS SPECTRAL DENSITY 
% ========================================================================

dt = t(2) - t(1) ;
fs = 1/dt ;

Gmat = cell(N,N) ;

%Npoints = [] ;
Npoints = 2^8 ;             % For estimation
%Npoints = 2^10 ;           % For visualization

for i = 1:1:N

    for j = 1:1:N

        %[G,f] = cpsd(x0(i,:),x0(j,:),[],3e3,[],fs) ; 
        [G,f] = cpsd(y0(i,:),y0(j,:),[],[],Npoints,fs) ;

        Gmat{i,j} = G'/(2*pi) ;

    end
end

w2 = (2*pi*f)' ;

figure(7) ; clf ; hold all ; box on
%set(gcf,'Position',[800 300 1200 800])
set(gcf,'Position',[800 300 1200 400])
%subplot(2,1,1,'FontSize',16) ; hold all ; box on 

% ------------------------------------------------------------------------
x_start = 0;
x_end = 0.15;
%y_min = 0 ;   % linear scale PSD
%y_max = 1e3 ; % linear scale PSD
%y_min = 1e-2; % log scale PSD
%y_max = 1e4;  % log scale PSD
y_min = 0;    % coherence  
y_max = 1.5;  % coherence 

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

% Coherence:
% ------------------------------------------------------------------------
plot(w2,abs(Gmat{1,2}).^2./(Gmat{1,1}.*Gmat{2,2}),'Color','#7E2F8E',...
    'LineWidth',1.5) 
plot(w2,abs(Gmat{1,3}).^2./(Gmat{1,1}.*Gmat{3,3}),'Color',"#77AC30",...
    'LineWidth',1.5) 
plot(w2,abs(Gmat{2,3}).^2./(Gmat{2,2}.*Gmat{3,3}),'Color',"#D95319",...
    'LineWidth',1.5) 
fontsize(fontS,'points')
ylim([0,1.5])
%xlim([0,0.1])
xlim([0,1])
title('\textbf{Experimental}')
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('Coherence, $\gamma^2_{ij}$ [-]')
legend('$\gamma^2_{12}$','$\gamma^2_{13}$','$\gamma^2_{23}$')

% PSD:
% ------------------------------------------------------------------------
% plot(w2,Gmat{1,1},'Color',color1,'LineWidth',1.5)
% plot(w2,Gmat{2,2},'Color',color2,'LineWidth',1.5)
% plot(w2,Gmat{3,3},'Color',color3,'LineWidth',1.5)
% fontsize(fontS,'points')
% yscale('log')
% %xlim([0,0.1])
% ylim([1e-2,1e4])
% %ylim([0,1e3])
% legend('$\hat{G}_{\tilde{x}_1\tilde{x}_1}(\omega)$',...
%     '$\hat{G}_{\tilde{x}_2\tilde{x}_2}(\omega)$',...
%     '$\hat{G}_{\tilde{x}_3\tilde{x}_3}(\omega)$')
% xlabel('Angular frequency, $\omega$ [rad/s]')
% ylabel('PSD, $G_{\tilde{x}_i\tilde{x}_i}(\omega)$ [m$^2$s/rad]')

%%
% ========================================================================
%                    CROSS SPECTRAL DENSITY PLOTS
% ========================================================================

figure(8) ; clf 
set(gcf,'Position',[200 300 1400 800])
subplot(2,1,1,'fontsize',fontS,'yscale','linear') ; hold all ; box on 
plot(w2,real(Gmat{1,2}),'-','Color','#7E2F8E','LineWidth',1.5)
plot(w2,real(Gmat{1,3}),'-','Color',"#77AC30",'LineWidth',1.5)
plot(w2,real(Gmat{2,3}),'-','Color',"#D95319",'LineWidth',1.5)
legend('Re$(G_{y_1y_2})$','Re$(G_{y_1y_3})$','Re$(G_{y_2y_3})$')
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('Real part of measurement CPSD')
xlim([0,0.2])

subplot(2,1,2,'fontsize',fontS,'yscale','linear') ; hold all ; box on 
plot(w2,imag(Gmat{1,2}),'-','Color','#7E2F8E','LineWidth',1.5)
plot(w2,imag(Gmat{1,3}),'-','Color',"#77AC30",'LineWidth',1.5)
plot(w2,imag(Gmat{2,3}),'-','Color',"#D95319",'LineWidth',1.5)
legend('Im$(G_{y_1y_2})$','Im$(G_{y_1y_3})$','Im$(G_{y_2y_3})$')
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('Imaginary part of measurement CPSD')
xlim([0,0.2])

figure(9) ; clf 
set(gcf,'Position',[200 300 1400 800])
subplot(2,1,1,'fontsize',fontS,'yscale','linear') ; hold all ; box on 
plot(w2,abs(Gmat{1,2}),'-','Color','#7E2F8E','LineWidth',1.5)
plot(w2,abs(Gmat{1,3}),'-','Color',"#77AC30",'LineWidth',1.5)
plot(w2,abs(Gmat{2,3}),'-','Color',"#D95319",'LineWidth',1.5)
legend('$|G_{y_1y_2}|$','$|G_{y_1y_3}|$','$|G_{y_2y_3}|$')
xlabel('Time, t [s]')
ylabel('Absolute value of CPSD')
xlim([0,0.2])

subplot(2,1,2,'fontsize',fontS,'yscale','linear') ; hold all ; box on 
plot(w2,phase(Gmat{1,2}),'-','Color','#7E2F8E','LineWidth',1.5)
plot(w2,phase(Gmat{1,3}),'-','Color',"#77AC30",'LineWidth',1.5)
plot(w2,phase(Gmat{2,3}),'-','Color',"#D95319",'LineWidth',1.5)
legend('$\theta_{y_1y_2}$','$\theta_{y_1y_3}$','$\theta_{y_2y_3}$')
xlabel('Time, $t$ [s]')
ylabel('phase of CPSD')
xlim([0,0.2])

%%
% ========================================================================
%                       PARAMETER ESTIMATION (1)
% ========================================================================

rng(1.0)

lb = 0.25*[1 ; 1] ;     ub = 5.0*[1 ; 1] ;
theta_init = lb + (ub - lb) .* rand(2,1) ;

%options.diff = 'mean' ;
options.idx = [1,2] ;
options.wcut = w2(end) ;

Q1 = @(theta) fun_obj_CPSD([A_ref;Dh_ref].*theta,M,Km,Dm,Gmat,w2,options) ;

% The iterative output is disabled by this line
fminoptions = optimoptions('fmincon','display','none') ;    
%fminoptions.OptimalityTolerance = 1e-6 ;
%fminoptions.FunctionTolerance = 1e-9 ;
%fminoptions.StepTolerance = 1e-8 ;

theta_opt = fmincon(Q1,theta_init,[],[],[],[],lb,ub,[],fminoptions) ;


%filename2 = 'theta_opt_paper.mat'; 
%save(filename2,'theta_opt'); 

%%
% ========================================================================
%                       PARAMETER ESTIMATION (2)
% ========================================================================

w_cut = 0.1 ;
%w_cut = 0.2 ;
%w_cut = 0.25 ;
%w_cut = w2(end) ;
list = abs(w)<w_cut ;

Y = zeros(N,length(y)) ;
Yin = cell(N,1) ;

for i = 1:1:N
    Y(i,:) = fftshift(fft(y(i,:))) ;
    Yin{i} = Y(i,:).*list ;
end

rng(1.0)
lb = 0.0*[1 ; 1] ;     ub = 5.0*[1 ; 1] ;
theta_init = lb + (ub - lb) .* rand(2,1) ;

Q2 = @(theta) fun_obj_VAR([A_ref;Dh_ref].*theta,M,Km,Dm,w,Yin) ;

% The iterative output is disabled by this line
fminoptions = optimoptions('fmincon','display','none') ;    
%options.OptimalityTolerance = 1e-20 ;
%options.StepTolerance = 1e2 ;

theta_opt_2 = fmincon(Q2,theta_init,[],[],[],[],lb,ub,[],fminoptions) ;

%%
% ========================================================================
%                   CROSS SPECTRAL DENSITY SCALING
% ========================================================================

p = options.idx(1) ;
q = options.idx(2) ;
% creating transfer functions 
% ------------------------------------------------------------------------
H_hat = cell(N,1) ;

for i = 1:1:N

    H_hat{i} =...
        (-w2.^2*(theta_opt(1)*A_ref + M) + ...
        1i*w2*(theta_opt(2)*Dh_ref + Dm(i)) + Km(i)).^-1 ;

        % H_hat{i} =...
        % (-w2.^2*(theta_opt_2(1)*A_ref + M) + ...
        % 1i*w2*(theta_opt_2(2)*Dh_ref + Dm(i)) + Km(i)).^-1 ;


end

% create cross spectra estimates
% ------------------------------------------------------------------------
%fac = cell(N,N) ;

Gpq_hat = cell(N,N) ;

for i = 1:1:N

    for j = (i+1):1:N

        fac = conj(H_hat{p}).*H_hat{q}./(conj(H_hat{i}).*H_hat{j}) ;

        Gpq_hat{i,j} = fac.*Gmat{i,j} ;

    end

end

figure(10) ; clf
set(gcf,'Position',[800 300 1400 800])
subplot(2,1,1,'fontsize',fontS) ; hold all ; box on 
plot(w2,abs(Gpq_hat{1,2}),'color',color12,'LineWidth',1.5)
plot(w2,abs(Gpq_hat{1,3}),'color',color13,'LineWidth',1.5)
plot(w2,abs(Gpq_hat{2,3}),'color',color23,'LineWidth',1.5)
xlim([0,0.1])
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('$|(\hat{G}_{pq}(\omega))_{ij}|$')
legend('$|(\hat{G}_{pq}(\omega))_{12}|$',...
    '$|(\hat{G}_{pq}(\omega))_{13}|$','$|(\hat{G}_{pq}(\omega))_{23}|$')
 
subplot(2,1,2,'fontsize',fontS) ; hold all ; box on 
plot(w2,phase(Gpq_hat{1,2}),'color',color12,'LineWidth',1.5)
plot(w2,phase(Gpq_hat{1,3}),'color',color13,'LineWidth',1.5)
plot(w2,phase(Gpq_hat{2,3}),'color',color23,'LineWidth',1.5)
xlim([0,0.1])
xlabel('Angular frequency, $\omega$ [rad/s]')
ylabel('$|(\hat{\theta}_{pq})_{ij}|$')
legend('$|(\hat{\theta}_{pq})_{12}|$',...
    '$|(\hat{\theta}_{pq})_{13}|$','$|(\hat{\theta}_{pq})_{23}|$')

%%
% ========================================================================
%                      EXCITATION LOAD ESTIMATION
% ========================================================================

yd = zeros(size(y)) ;  Yd = zeros(size(y)) ; 
ydd = zeros(size(y)) ; Ydd = zeros(size(y)) ;  

%w_cut = 0.20 ;
w_cut = w2(end) ;

list = abs(w)<w_cut ;

Y = Y.*list ;

for i = 1:1:N
 
    Yd(i,:) = 1i*w.*Y(i,:) ;            % Yd(i,1) = 0 ;
    Ydd(i,:) = -w.^2.*Y(i,:) ;          % Ydd(i,1) = 0 ;

    y(i,:) = ifft(ifftshift(Y(i,:))) ;
    yd(i,:) = ifft(ifftshift(Yd(i,:))) ;
    ydd(i,:) = ifft(ifftshift(Ydd(i,:))) ;

end
 
f_hat = (M + A_ref*theta_opt(1)).*ydd + ...
    (Dh_ref.*theta_opt(2) + Dm').*yd + Km'.*y ;
 
% f_hat = (M + A_ref*theta_opt_2(1)).*ydd + ...
%      (Dh_ref*theta_opt_2(2) + Dm').*yd + Km'.*y ;

%filename2 = 'force_recon_paper.mat' ; 
%save(filename2,'theta_opt','t','f_hat'); 
%%
% ========================================================================
%                           LOAD SPECTRA
% ========================================================================

f_mean = mean(f_hat) ;
e = f_mean - f_hat ;

figure(11) ; clf 
set(gcf,'Position',[800 300 1600 1000])
subplot(3,1,1,'fontsize',fontS); hold all ; box on 
plot(w,abs(Y(1,:)),'color',color1)
plot(w,abs(Y(2,:)),'color',color2)
plot(w,abs(Y(3,:)),'color',color3)
xlabel('Angular frequnecy, $\omega$ [rad/s]')
ylabel('$|Y(\omega)|$')
legend('$|Y_1(\omega)|$','$|Y_2(\omega)|$','$|Y_3(\omega)|$')

subplot(3,1,2,'fontsize',fontS); hold all ; box on 
plot(w,abs(Yd(1,:)),'color',color1)
plot(w,abs(Yd(2,:)),'color',color2)
plot(w,abs(Yd(3,:)),'color',color3)
xlabel('Angular frequnecy, $\omega$ [rad/s]')
ylabel('$|i\omega Y(\omega)|$')
legend('$|i\omega Y_1(\omega)|$','$|i\omega Y_2(\omega)|$',...
    '$|i\omega Y_3(\omega)|$')

subplot(3,1,3,'fontsize',fontS); hold all ; box on 
plot(w,abs(Ydd(1,:)),'color',color1)
plot(w,abs(Ydd(2,:)),'color',color2)
plot(w,abs(Ydd(3,:)),'color',color3)
xlabel('Angular frequnecy, $\omega$ [rad/s]')
ylabel('$|-\omega^2Y(\omega)|$')
legend('$|-\omega^2Y_1(\omega)|$','$|-\omega^2Y_2(\omega)|$',...
    '$|-\omega^2Y_3(\omega)|$')

figure(12) ; clf 
set(gcf,'Position',[800 300 1400 800])
subplot(2,1,1,'fontsize',fontS) ; hold all ; box on 
plot(t,f_hat(1,:),'color',color1,'LineWidth',1.5)
plot(t,f_hat(2,:),'color',color2,'LineWidth',1.5)
plot(t,f_hat(3,:),'color',color3,'LineWidth',1.5)
plot(t,f_mean(1,:),'k')
xlabel('Time, $t$ [s]')
ylabel('Force estimate, $\hat{f}(t)$ [N]')
legend('$\hat{f}_1(t)$','$\hat{f}_2(t)$','$\hat{f}_3(t)$',...
    '$\bar{f}(t)$','Location','southeast')
xlim([3000,11000])

subplot(2,1,2,'fontsize',fontS) ; hold all ; box on 
plot(t,e(1,:),'color',color1)
plot(t,e(2,:),'color',color2)
plot(t,e(3,:),'color',color3)
xlabel('Time, $t$ [s]')
ylabel('Error, $e(t)$')
legend('$e_1(t)$','$e_2(t)$','$e_3(t)$')

figure(13) ; clf ; hold all ; box on
set(gcf,'Position',[600 100 1200 400])
%subplot(3,1,1,'fontsize',fontS) ; hold all ; box on 
%histogram(e(1,:),'FaceColor',color1)
[ks1,h1] = ksdensity(e(1,50:end-50)) ;
plot(h1/1000,ks1,'-','LineWidth',1.5,'Color',color1)
[ks2,h2] = ksdensity(e(2,50:end-50)) ;
plot(h2/1000,ks2,'-','LineWidth',1.5,'Color',color2)
[ks3,h3] = ksdensity(e(3,50:end-50)) ;
plot(h3/1000,ks3,'-','LineWidth',1.5,'Color',color3)
fontsize(fontS,'points')
ylabel('Density')
xlabel('Estimated error, $\hat{v}(t)$ [kN]')
xlim([-1.5e5/1000,1.5e5/1000])
ylim([0,2e-5])
legend('$\hat{v}_1(t)$','$\hat{v}_2(t)$','$\hat{v}_3(t)$')

Npoints = 2^10;

f0_mean = f_mean - mean(f_mean) ;

f0 = f_hat - mean(f_hat,2) ;

[Gff_mean,ff] = cpsd(f0_mean/1000,f0_mean/1000,[],[],Npoints,fs) ;
wf = ff*2*pi ;

Gff = zeros(N,length(Gff_mean))' ;
Gee = zeros(N,length(Gff_mean))' ;

for i = 1:1:N

    Gff(:,i) = cpsd(f0(i,:)/1000,f0(i,:)/1000,[],[],Npoints,fs) ;
    Gee(:,i) = cpsd(e(i,:)/1000,e(i,:)/1000,[],[],Npoints,fs) ;

end

Gff_mean = Gff_mean'/(2*pi) ;
Gff = Gff'./(2*pi) ;
Gee = Gee'./(2*pi) ;

%filename2 = 'PSDs_cut_off_cpsd.mat'; 
%save(filename2,'Gee','Gff','Gff_mean','wf'); 

figure(14) ; clf ; hold all ; box on 
set(gcf,'Position',[600 100 1200 600])
plot(wf,Gff(1,:),'--','Color',color1,'LineWidth',1.5)
plot(wf,Gff(2,:),'--','Color',color2,'LineWidth',1.5)
plot(wf,Gff(3,:),'--','Color',color3,'LineWidth',1.5)
plot(wf,Gff_mean,'k','LineWidth',1.5)
xlim([0,0.2])
ylim([0,1.0e6])
xlabel('Angular frequency, $\omega$ [rad/s]')
fontsize(fontS,'points')
ylabel('LF load PSD, $G_{ff}(\omega)$ [kNs/rad]')
legend('$\hat{G}_{f_1f_1}(\omega)$',...
    '$\hat{G}_{f_2f_2}(\omega)$',...,
    '$\hat{G}_{f_3f_3}(\omega)$','$\hat{G}_{ff}(\omega)$')

figure(15) ; clf ; hold all ; box on 
set(gcf,'Position',[600 100 1200 600])
plot(wf,Gee(1,:),'--','Color',color1,'LineWidth',1.5)
plot(wf,Gee(2,:),'--','Color',color2,'LineWidth',1.5)
plot(wf,Gee(3,:),'--','Color',color3,'LineWidth',1.5)
plot(wf,mean(Gee),'Color','k','LineWidth',1.5)
xlim([0,0.2])
ylim([0,6e4])
xlabel('Angular frequency, $\omega$ [rad/s]')
fontsize(fontS,'points')
ylabel('Load error PSD, $G_{vv}(\omega)$ [kNs/rad]')
legend('$\hat{G}_{v_1v_1}(\omega)$',...
    '$\hat{G}_{v_2v_2}(\omega)$',...,
    '$\hat{G}_{v_3v_3}(\omega)$',...,
    '$\hat{G}_{vv}(\omega)$') 

figure(16) ; clf ; hold all ; box on 
set(gcf,'Position',[600 100 1200 600])
plot(wf,Gff_mean,'k','LineWidth',1.5)
plot(wf,mean(Gee),'Color',[0.7 0.7 0.7],'LineWidth',1.5)
xlim([0,0.2])
ylim([0,1.0e6])
xlabel('Angular frequency, $\omega$ [rad/s]')
fontsize(fontS,'points')
ylabel('LF load PSD, $G_{ff}(\omega)$ [kNs/rad]')
yscale('linear')


%%
% ========================================================================
%                           CUT-OFF SENSITIVITY
% ========================================================================

Window = 0.01:0.01:0.24 ;

for j = 1:1:length(Window) 

% CPSD-BASDED PARAMETER ESTIMATION
% ------------------------------------------------------------------------
options.wcut = Window(j) ;

Q1 = @(theta) fun_obj_CPSD([A_ref;Dh_ref].*theta,M,Km,Dm,Gmat,w2,options) ;


% The iterative output is disabled by this line
fminoptions = optimoptions('fmincon','display','none') ;    
%options.OptimalityTolerance = 1e-20 ;
%options.StepTolerance = 1e2 ;

theta_opt_win(:,j) = fmincon(Q1,theta_init,[],[],[],[],lb,ub,...
    [],fminoptions) ;

disp(num2str(j))

% VARIANCE-BASED PARAMETER ESTIMATION
% ------------------------------------------------------------------------

list = abs(w)<Window(j) ;

for i = 1:1:3

    Yin{i} = Y(i,:).*list ;

end

Q2 = @(theta) fun_obj_VAR([A_ref;Dh_ref].*theta,M,Km,Dm,w,Yin) ;

theta_opt_2 = fmincon(Q2,theta_init,[],[],[],[],lb,ub,[],fminoptions) ;

% Optimization
theta_opt_win2(:,j) = fmincon(Q2,theta_init,[],[],[],[],...
    lb,ub,[],fminoptions) ; 


end

figure(17) ; clf ; hold all ; box on %; grid on
% ------------------------------------------------------------------------
%set(gcf,'Position',[1600 300 1200 400])
set(gcf,'Position',[1600 300 800 400])
 plot(Window,theta_opt_win2(1,:)*A_ref/1000,'x',...
     'MarkerEdgeColor',"#D95319",'LineWidth',1.5,'MarkerSize',10)

plot(Window,theta_opt_win(1,:)*A_ref/1000,'o',...
    'MarkerEdgeColor',"#0072BD",'LineWidth',1.5)

ylim([0,3e4])
fontsize(fontS,'points')
xlabel('$w_{cut}$ [rad/s]')
ylabel('Added mass estimate, $\hat{a}$ [t]')
legend('S2021 Estimator','C2025 Estimator')
%legend('Estimates')

figure(18) ; clf ; hold all ; box on %; grid on
% ------------------------------------------------------------------------
set(gcf,'Position',[1600 300 800 400])
plot(Window,theta_opt_win2(2,:)*Dh_ref/1000,'x',...
     'MarkerEdgeColor',"#D95319",'LineWidth',1.5,'MarkerSize',10)

plot(Window,theta_opt_win(2,:)*Dh_ref/1000,'o',...
    'MarkerEdgeColor',"#0072BD",'LineWidth',1.5)

ylim([0,1e3])
fontsize(fontS,'points')
xlabel('$w_{cut}$ [rad/s]')
ylabel('Damping estimate, $\hat{d}_h$ [kNs/m]')
legend('S2021 Estimator','C2025 Estimator')
%legend('Estimates')

%filename2 = 'convergence_exp_paper.mat' ; 
%save(filename2,'A_ref','Dh_ref','theta_opt_win','theta_opt_win2','Window'); 