% This input file contains the numerical reference of the three
% experimental tests (AsRun-ED) used in the paper.

% Set text format 
% ------------------------------------------------------------------------
set(0,'defaulttextinterpreter','latex') ;
set(0,'defaultAxesTickLabelInterpreter','latex'); 
set(0,'defaultLegendInterpreter','latex');

fontS = 16 ;

% Directory
% ------------------------------------------------------------------------
mydir = fileparts(pwd) ;

addpath([mydir '\num\functions'])
addpath([mydir '\num\hydrostar'])
addpath([mydir '\num\loads'])

% For preliminary analysis and OED comparison
filename = "loads_T_max_12000_JS_Tp_12_Hs_11_Head_1.mat" ;

% For consistency check 
% filename = "loads_T_max_4000_JS_Tp_12_Hs_11_Head_1.mat" ;
% filename = "loads_T_max_8000_JS_Tp_12_Hs_11_Head_1.mat" ;
% filename = "loads_T_max_16000_JS_Tp_12_Hs_11_Head_1.mat" ;
% filename = "loads_T_max_32000_JS_Tp_12_Hs_11_Head_1.mat" ;

load(filename)
clear F2 F1 

% Wave spectrum                 
% ------------------------------------------------------------------------
r = length(t) ;
dt = t(2) - t(1) ;        % time increment [s]
Tmax = t(end) + dt;       % length of signal time domain [s] 
wnyq = pi/dt ;            % Nyquist frequency [rad/s]
dw = 2*pi/Tmax ;          % frequency increment [rad/s]
w = dw:dw:wnyq ;          % frequency vector for one sided spectrum
STD = 0.7 ;

numExtrac = regexp(filename, '\d+\.\d+|\d+', 'match');
numExtrac = str2double(numExtrac);

Z = fun_spectrum_jon(numExtrac(2), numExtrac(3), w, STD)'; 
Z = [fliplr(conj(Z)) , Z] ;
clear w 

% Preparing frequency-domain wave loads                 
% ------------------------------------------------------------------------

w = (-r/2:r/2-1)*2*pi/Tmax ;
w_0 = 0.25 ;

A_ref = 1.0e4 ;         % added mass "true value"
Dh_ref = 700 ;          % damping "true value"

CM = read_hs_CA_CM('CM_11.rao') ;
CA = read_hs_CA_CM('CA_11.rao') ;

idx = find(w == 0) ;

A_WF = interp1(CM.freq,CM.val(1,:),w(idx:end),'linear','extrap') ;
A_WF = A_WF/1e3 ;
D_WF = interp1(CA.freq,CA.val(1,:),w(idx:end),'linear','extrap') ;
D_WF = D_WF/1e3 ;

A_dbl = [A_WF(end) , fliplr(A_WF(2:end)) , A_WF] ;
D_dbl = [D_WF(end) , fliplr(D_WF(2:end)) , D_WF] ;

% Floater parameters 
% ------------------------------------------------------------------------

M = 1.4176e4  ;                      % mass                    [t]

% AsRun-ED
K = [53.1 90.2 53.1]  ;              % mooring stiffness (3)   [kN/m]       
D = [217.8 0 0]  ;                   % mooring damping (3)     [kNs/m]

% Poor-ED
%K = [70.0 80.0 90.0]  ;              % mooring stiffness (3)   [kN/m]
%D = [100.0 0 0]  ;                   % mooring damping (3)     [kNs/m]

% SpaceFilling-ED
%K = [40.0 80.0 60.0]  ;              % mooring stiffness (3)   [kN/m]
%D = [0 150.0 250.0]  ;               % mooring damping (3)     [kNs/m]

% Optimal-ED
%K = [150.0 40.0 150.0]  ;            % mooring stiffness (3)   [kN/m]
%D = [0 0 0]  ;                       % mooring damping (3)     [kNs/m]

N = length(K) ;

% Noise
% ------------------------------------------------------------------------

noise.loadSNR = 10 ;
noise.SNR = 35 ;        noise.opt = '1' ;       % '1' or '2' (dB)

% figure(100) ; clf ; hold on ; box on 
% set(gcf,'Position',[600 100 1000 800])
% subplot(2,1,1,'fontsize',16) ; hold on ; box on 
% plot(w,A_dbl,'LineWidth',1)
% xlim([0,2])
% ylabel('Added mass [kg]')
% xlabel('$\omega$ [rad/s]')
% 
% subplot(2,1,2,'fontsize',16) ; hold on ; box on 
% plot(w,D_dbl,'LineWidth',1)
% xlim([0,2])
% ylabel('Radiation damping [Ns/m]')
% xlabel('$\omega$ [rad/s]')