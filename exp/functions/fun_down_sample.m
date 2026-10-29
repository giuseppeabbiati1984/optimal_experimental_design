function [t,y] = fun_down_sample(t_resampled,y_resampled,fs_resampled,DSopt)

% ========================================================================
%                        DOWNSAMPLED SIGNALS
% ========================================================================

N = height(y_resampled) ;
fontS = 16 ;

fac = DSopt.fac ;

fs_down = fs_resampled/fac ;

f_cut = fs_down / 2;                     % Cutoff frequency (Hz)

filt_type = DSopt.filt_type ;

y_filt = zeros(size(y_resampled)) ;

switch filt_type 

    case 'butter' 
        [a, b] = butter(6, f_cut /(fs_resampled / 2));                      

        figure(5) ; clf ; hold on ; box on 
        set(gcf,'Position',[800 300 1400 800])
        fontsize(fontS,'points')
        freqz(a,b,[],fs_resampled)
        %xlim([0,1])
        
        % Apply the low-pass filter to the signal
        
        for i = 1:1:N
        
            % Zero-phase filtering to avoid phase distortion
              y_filt(i,:) = filtfilt(a, b, y_resampled(i,:));  
        
        end 

    case 'fir'                                    

        ord = 6000;                                                         % Filter order (adjust for sharper roll-off)
        b = fir1(ord, f_cut/(fs_resampled / 2), 'low');                     % Design the FIR low-pass filter
        
        figure(5) ; clf ; hold on ; box on 
        set(gcf,'Position',[800 300 1600 1000])
        fontsize(fontS,'points')
        freqz(b,1,[],fs_resampled)
        %xlim([0,0.1])

        for i = 1:1:N
            y_filt(i,:) = filtfilt(b, 1, y_resampled(i,:)) ;

        end
end

y = y_filt(:,1:fac:end) ;
t = t_resampled(:,1:fac:end) ;