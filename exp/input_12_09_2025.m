% Set text format 
% ------------------------------------------------------------------------
set(0,'defaulttextinterpreter','latex') ;
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');

% ========================================================================
%                               LOAD DATA
% ========================================================================

mydir = fileparts(pwd) ;
addpath([mydir '\exp\data'])
addpath([mydir '\exp\functions'])

testNum = ["4090", "4072", "4060"] ; % I_tot = {1,2,3}

fontS = 16 ;
color1 = '#0072BD' ;        % Blue
color2 = '#A2142F' ;        % Red
color3 = '#EDB120' ;        % Yellow

color12 = '#7E2F8E' ;        % purple
color13 = '#77AC30' ;        % green
color23 = '#D95319' ;        % orange

N = numel(testNum) ;

%%
% ========================================================================
%                  PREPARING PARAMETERS FOR OPTIMIZATION 
% ========================================================================

% Mooring parameters for {4090 4072 4060}
% ------------------------------------------------------------------------
Km_m = [32.4 55.0 32.4] ;
Km = 40^2*Km_m*(1.025) ;

Dm_m = [21.0 0 0] ;
Dm = 40^(2.5)*Dm_m*(1.025) ;

M =  1.4176e7 ;	
A_ref = 1e7 ;
Dh_ref = 1e6 ;

%%
% ========================================================================
%                               RAW DATA
% ========================================================================

load('expData.mat') 

figure(1) ; clf ; hold all ; box on 
set(gcf,'Position',[800 300 1200 600])
plot(t_raw{1},y_raw{1},'color',color1)
plot(t_raw{2},y_raw{2},'color',color2)
plot(t_raw{3},y_raw{3},'Color',color3)
fontsize(fontS,'points')
ylabel('Raw displacement measurements, $y$ [m]')
xlabel('Time, $t$ [s]')
legend('$y_1$','$y_2$','$y_3$')

fs_raw = [1/(t_raw{1}(2) - t_raw{1}(1))...
              1/(t_raw{2}(2) - t_raw{2}(1))...
              1/(t_raw{3}(2) - t_raw{3}(1))] ;

fs_max = max(fs_raw) ;