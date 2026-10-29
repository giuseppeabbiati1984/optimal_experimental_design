function [t_resampled,y_resampled] = fun_align_resp(y_raw,fs_raw,s_idx)

% ========================================================================
%                            ALIGNING SIGNALS
% ========================================================================
N = length(y_raw) ;

fs_resampled = min(fs_raw);   % Target sampling frequency

% Get rational approximation and resample (now low-pass filtering)
% ------------------------------------------------------------------------

y_temp = cell(N,1) ;
Npoints = zeros(N,1) ;

for i = 1:1:N

    % Resample signals
    [p1, q1] = rat(fs_resampled / fs_raw(i));  
    y_temp{i} = resample(y_raw{i}, p1, q1);

    Npoints(i) = length(y_temp{i}) ;
end

% Check if the lengths match and trim if necessary
[min_length,idx] = min(Npoints);

t_resampled = (0:length(y_temp{idx})-1) / fs_resampled;

y_resampled = zeros(N,min_length) ;

for i = 1:1:N
    
    y_resampled(i,:) = y_temp{i}(1:min_length)';

end

t_resampled = t_resampled(s_idx:end-10) ;
y_resampled = y_resampled(:,s_idx:end-10) ;