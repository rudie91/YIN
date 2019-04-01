function [pitch, ap_pwr, RMS] = myYin(wavData, fs, tHop, tW, f0Min, f0Max, dp_th)
% My YIN implementation of Step 1-5.
% The first estimate is performed at time 0 (i.e. the first integration
% window is from the first sample to time tW, which means the window is
% centered at tW/2) and the following estimates are spaced by tHop.
%
% Input
%   - wavData       : input single-channel audio wave (a column vector)
%   - fs            : sampling rate of wavData
%   - tHop          : time interval between two adjacent estimates in second (default 0.01)
%   - tW            : integration window size in second (default 0.025)
%   - f0Min         : lowest possible F0 in Hz (default 40)
%   - f0Max         : highest possible F0 in Hz (default 400)
%   - dp_th         : the threshold of dips of d_prime(default 0.1)
% Output
%   - pitch         : estimated pitches (a row vector)
%   - ap_pwr        : the corresponding d_prime, which is approximatedly
%                       the aperiodic power over the total signal power. 
%                       It can be used as the salience of this pitch
%                       estimate.
%   - rms           : the RMS value of the signal in the integration
%                       window, which can also be used to determine if
%                       there is a pitch in the signal or not.
%
% Author: Yoon mo Yang
% Created: Sep 16th 2018
% Last modified: Sep 19th 2018

% default parameters for speech
if nargin<7 dp_th=0.1; end
if nargin<6 f0Max=400; end
if nargin<5 f0Min=40; end
if nargin<4 tW=0.025; end
if nargin<3 tHop=0.01; end

% Start your implementation here
W = floor(tW * fs); % Window length in samples 
R = floor(fs * tHop); % Hop size in samples
nframes = floor((length(wavData)-W)/R); % number of frames based on input, W, R

for i = 1:nframes
    fprintf("The %i th frame is being processed now\n",i);
    d_t = zeros(W+1,1); % Initialize difference eqn d_t
    d_prime = zeros(W+1,1); % Initialize normalized difference eqn d_prime
    RMS(i) = rms(wavData((i-1)*R+1:(i-1)*R+W)); % RMS value of windowed signal
    for tau = 1:W+1 
        % Step 2: difference eqn
        for j = (i-1)*R+1:(i-1)*R+W
            if j+(tau-1) > length(wavData) % for the last frame
                break;
            else
                d_t(tau) = (wavData(j) - wavData(j+(tau-1)))^2 + d_t(tau); % culumulative sum of diff eqn
            end
        end
        % Step 3: Normalization      
        d_prime(tau) = (tau-1)*d_t(tau)/sum(d_t(1:tau)); % normalized diff eqn
    end
    % Step 4: Absolute Threshold
    tau_min = min(find(d_prime < dp_th)); % find the tau value that gives the d_prime less than dp_th
    if isempty(tau_min) % if there is no tau value that satisfies the condition above
        [d_min,tau_min] = min(d_prime); % global minimum
        ap_pwr(i) = d_prime(tau_min); % ap_pwr 
        pitch(i) = fs/(tau_min); % pitch with the global minimum
% 1 part (c): Before modifying the frame w/o pitch
%         pitch(i) = fs/(tau_min);
%         if pitch(i) < f0Min || f0Max < pitch(i) % Check the freq range 
%             pitch(i) = 0;
%         end

% After modifying the frame w/o pitch by setting the pitch to 0 (w/
% thresholds for ap_pwr and RMS)
        if (ap_pwr(i) >= 0.15 && RMS(i) <= 0.05) || (pitch(i) < f0Min || f0Max < pitch(i))
            pitch(i) = 0; 
        end
    else
        % Step 5: Parabolic interpolation
        a = tau_min-1; % the previous bin
        b = tau_min; % the current bin
        c = tau_min+1;  % the future bin
        alpha = d_prime(a); % the previous value (alpha)
        beta = d_prime(b); % the current value (beta)
        gamma = d_prime(c); % the next value (gamma)
        
        % update tau_min using  parabolic interpolation (int_min:interpolated minimum)
        int_min = floor(b + 0.5*((alpha-beta)*(c-b)^2-(gamma-beta)*(b-a)^2)/((alpha-beta)*(c-b)-(gamma-beta)*(b-a)));
        pitch(i) = fs/(int_min); % update the pitch with interpolated tau_min
        ap_pwr(i) = d_prime(int_min); % update the ap_pwr with interpolated tau_min
        if pitch(i) < f0Min || f0Max < pitch(i) % Check the freq range 
            pitch(i) = 0;
        end
    end
end % for loop (nframe)
end % end of function
