% Oliver Gordon & Dominic Coy, 2017
% Is Drumming Fractal? - Main Code

%% Setup/Settings
%clear variables;
%close all;

% Settings
sample_start = 1;                                   % Sample to start from
sample_end = 2;                                     % Sample to end with                                
technique = 'Double Handed';                        % Double Handed OR 
                                                    % Single Handed OR 
                                                    % Metronomic

bpm = 88;                                           % BPM of piece

filter_order = 100;                                 % FIR filter order
max_pass_freq = 16000;                              % Max bandpass freq
min_pass_freq = 600;                                % Min bandpass freq
hard_freq = 14000;                                  % Max hardpass freq

basic_thresh_peak_no = 100;                         % Simple onset drift
basic_thresh_peak = 15;                             % Simple onset amp

thresh_onset = 0.034;                               % Max onset drift
thresh_amp = 0.1;                                   % Max onset amp
s_min = 5;                                          % Min window size
s_max = 100;                                        % Max window size

db1_end = 20;                                       % Max window size H1
db2_start = 30;                                     % Min window size H2
db_thresh = 0.1;                                    % Straight/Barrelled?

db = 1;
exp_onset = 60/bpm/4;

% Clear sound if sound still playing
if exist('player','var') == 1
    stop(player)
end

% Ensure 2016b or newer because of syntax change
if verLessThan ('matlab','9.1.0.441655') ~= 0
    error('Circshift changed syntax after 2016a. Use 2016b or newer.')
end

% Read in Rasenen data for checking |Onset Times|Amplitudes|
 all_data = dlmread('Docs\Check File.txt');

% Set up matrix to store exponents. Rows = samples, Columns = double
% barrelled test (11 if forward db, 12 if reverse, 13 if straight),
% exponent 1, error, exponent 2, error, exponent straight, error
final_result = NaN(7,(sample_end-sample_start+1));

% Set up matrices to store onsets and power variation. Columns = sample no.
final_onset_difs = NaN((sample_end-sample_start+1),3000);
final_power_difs = NaN((sample_end-sample_start+1),3000);

% Repeat entire process for all samples
for sample_no = sample_start:sample_end
%% Read In Files

% Read in raw audio file and make mono
Sample_Name = [technique,' ',num2str(sample_no)];   
[audio_data,audio_freq] = audioread(['Samples\',Sample_Name,'.wav']);
[~,channels] = size(audio_data);
if channels == 2
    audio_data = (audio_data(:,1) + audio_data(:,2)) / 2;
end

%% Pre-Processing

% 100th order FIR bandpass for 0.6-16kHz to remove low and high noise
[b,a] = fir1(filter_order,[min_pass_freq/(audio_freq/2),...
    max_pass_freq/(audio_freq/2)],'bandpass');
audio_data_FIR = filter(b, a, audio_data);

% Calculate spectogram
[~,spec_freq,spec_time,power] = ...
    spectrogram(audio_data_FIR, 64, [], [], audio_freq, 'yaxis');

% Hard lowpass for <14kHz so only well defined power changes will be used
% to calculate onset times
hard_ind = round((hard_freq/22500)*length(spec_freq));
F_hard = spec_freq((1:hard_ind),:);

% Calculate corresponding power
power = 10*log10(abs(power)+eps);
power = power(1:hard_ind,:);
power_tot = mean(power);

%% Onset Times - Basic Algorithm

% Define time axes
time = 1/audio_freq : 1/audio_freq : length(audio_data)/audio_freq;
time_power = linspace(1/audio_freq,length(audio_data)/audio_freq, ...
    length(power_tot));

Npeaksmax = round(length(audio_data_FIR)/audio_freq/60*bpm*4);

% Find locations of peaks & delete peaks incorrectly appearing during onset
    locs = find(power_tot - (circshift(power_tot,-1)) < -basic_thresh_peak)+1;
    peak_power = power_tot(locs);

% Calculate corresponding times between peaks
onset_time = time_power(locs);
onset_dif = diff(onset_time);
mean_onset_dif = mean(onset_dif);

%% Onset Times - Stronger Algorithm

% Calculate array number corresponding to window time
exp_onset_samples = round(exp_onset/(time_power(2)-time_power(1)));
window_onset_samples = round(thresh_onset/(time_power(2)-time_power(1)));

% Find first peak using basic algorithm to begin more advanced peak
% detection - code runs in <0.02 seconds, so don't bother adapting for just
% first peak
window_mid = locs(5);

% Set up loop
window_no = 1;
locs1 = NaN(1,Npeaksmax);
locs1(1) = window_mid;

while window_mid < length(time_power)-window_onset_samples
    % Define window range to search in
    window_range = window_mid-window_onset_samples ...
        : window_mid+window_onset_samples;
    
    % Detect potential onset peaks and record the highest, or carry on if
    % none found
    [~,peak_loc] = findpeaks(power_tot(window_range),'sortstr','descend');
    if isempty(peak_loc) == 0
        locs1(window_no) = window_mid-window_onset_samples+peak_loc(1)-1;
    end
    
    % Move to next windowing region
    window_mid = locs1(window_no)+exp_onset_samples;
    window_no = window_no+1;
end

% Calculate time and amplitude of peaks for each time axis
locs1 = locs1(1:window_no-2);
peak_power1 = power_tot(locs1);
onset_time1 = time_power(locs1);
locs2 = round(onset_time1./(time(2)-time(1)));
onset_time2 = time(locs2-(filter_order/2));
onset_amp2 = audio_data(locs2-(filter_order/2));

% Calculate difference between onset times, and remove %erroneous results
delta_t = diff(all_data(:,1))';             % Test onset differences
delta_p = diff(all_data(:,2))';
%delta_t = step(dsp.ColoredNoise(...         % Noise of chosen B.
%     'InverseFrequencyPower',0.5,...        % a = (B+1)/2
%    'SamplesPerFrame',5000))';
%delta_t = exp_onset.*ones(5000,1);          % Metronome

%delta_t = diff(onset_time1);                 % Actual data
%delta_p = diff(peak_power1);

% Remove errenous onset differences
delta_t(exp_onset-delta_t>thresh_onset | ...
    exp_onset-delta_t<-thresh_onset)=[];
delta_p(delta_t(exp_onset-delta_t>thresh_onset | ...
    exp_onset-delta_t<-thresh_onset)) = [];

% Add random offset for error
random_offset = (0.01+0.01).*rand(1,length(delta_t)) - 0.01;
%delta_t = delta_t+random_offset;

mean_t = mean(delta_t);

% Store power difference and onset time difference
final_onset_difs(sample_no,1:length(delta_t)) = delta_t;
final_power_difs(sample_no,1:length(delta_p)) = delta_p;

%% DFA

% Pre-allocate array of F_k
F_k_s = zeros(floor(length(delta_t)/s_max),s_max);

for s = s_min:s_max
    
    % Pre-allocate/reset inner loops
    delta_y = zeros(1,length(delta_t));
    fitrange = 1:s:length(delta_y)-s;
    delta_y_s = zeros(1,fitrange(end)+s-1);
    F_k_substep = zeros(1,length(s));
    
    % Sum onset difference over i from i=1 to i=N data points
    for i = 1:length(delta_t)
        % Sum over j from j=1 to j=i
        delta_y(i) = sum(delta_t(1:i)-mean_t);
    end
    
    % Omit data when not enough present to form new window of size s
    delta_y = delta_y(1:length(delta_y_s));
    
    % Perform linear fits for onset differences in windows of size s
    for fitx = 1:s:length(delta_y)
        polyco = polyfit_s(1:s , delta_y(fitx:fitx+s-1) , 1);
        delta_y_s(fitx:fitx+s-1) = ...
            (polyco(1).*(1:s))+polyco(2);
    end
    
    % Calculate RMS fluctuations of linear detrending
    for k_range = 0:((length(delta_y_s)/s)-1)
        sumstart = (k_range*s)+1;
        sumend = (k_range*s)+s;
        
        F_k_substep(k_range+1) = sum((delta_y(sumstart:sumend) - ...
            delta_y_s(sumstart:sumend)).^2);
    end
    
    F_k_s(1:length(F_k_substep),s) = sqrt((1/s).*F_k_substep); 
end

% Calculate mean (ignoring leftover values = 0) and define axis for cftool
F_s = (sum(F_k_s,1)./sum(F_k_s~=0,1));
s_axis = (1:s);

%% Find Hurst Exponent

% Remove NaN data and configure axis for use in curve fit
F_s(isnan(F_s))=[];
F_s = log(F_s);
s_axis = log(s_axis(end-length(F_s)+1:end));

% Curve fit first line of window sizes to first order
[pcoeff1,pstruct1] = polyfit((s_axis(1:db1_end)),(F_s(1:db1_end)),1);
[yfit1,perror1] = polyval(pcoeff1,log(s_axis(1:db1_end)),pstruct1);
perror1 = mean(perror1);

% Curve fit second line of window sizes to first order
[pcoeff2,pstruct2] = polyfit((s_axis(db2_start:end)),(F_s(db2_start:end)),1);
[yfit2,perror2] = polyval(pcoeff2,log(s_axis(db2_start:end)),pstruct2);
perror2 = mean(perror2);

% See if double barrelled by checking gradients+errors against threshold
% and store results
if (pcoeff2(1)-perror2) - (pcoeff1(1)+perror1) > db_thresh
    % Forward bias double barrel
    
    % Print values
    disp('Forward Double Barrelled')
    disp(['Hurst Exponent 1 = ',num2str(pcoeff1(1)),...
        char(177),num2str(sqrt(perror1.^2+0.02^2))])
    disp(['Hurst Exponent 2 = ',num2str(pcoeff2(1)),...
        char(177),num2str(sqrt(perror2.^2+0.02^2))])
    
    % Store values
    final_result(1,sample_no) = 11;
    final_result([2,4],sample_no) = [pcoeff1(1);pcoeff2(1)];
    final_result([3,5],sample_no) = [perror1;perror2];
elseif (pcoeff1(1)-perror1) - (pcoeff2(1)+perror2) > db_thresh
    % Reverse bias double barrel
    
    % Print values
    disp('Reverse Double Barrelled')
    disp(['Hurst Exponent 1 = ',num2str(pcoeff1(1)),...
        char(177),num2str(perror1)])
    disp(['Hurst Exponent 2 = ',num2str(pcoeff2(1)),...
        char(177),num2str(perror2)])
    
    % Store values
    final_result(1,sample_no) = 12;
    final_result([2,4],sample_no) = [pcoeff1(1);pcoeff2(1)];
    final_result([3,5],sample_no) = [perror1;perror2];
else
    % Straight line
    disp('Straight')
    
    % Fit straight line
    [pcoeff3,pstruct3] = polyfit((s_axis),(F_s),1);
    [yfit3,perror3] = polyval(pcoeff2,log(s_axis),pstruct2);
    perror3 = mean(perror3);  
    disp(['Hurst Exponent = ',num2str(pcoeff3(1)),...
        char(177),num2str(perror3)])
    
    % Store value
    final_result(1,sample_no) = 13;
    final_result([6,7],sample_no) = [pcoeff3(1);perror3];
    db=0;
end

%% Plotting

% Clear all axis
%arrayfun(@cla,findall(0,'type','axes'))

% ------
% Figure 1
% ------
% Plot graph of gradient Hurst Coefficient
Fig1 = figure(1);
hold on
plot(s_axis,F_s,'kx')

% Plot Hurst exponents
if db == 0
    plot(s_axis,((s_axis.*pcoeff3(1))+pcoeff3(2)))
    db = 1;
else
    plot(s_axis(1:db1_end),...
        ((s_axis(1:db1_end).*pcoeff1(1))+pcoeff1(2)))
    plot(s_axis(db2_start:end),...
        ((s_axis(db2_start:end).*pcoeff2(1))+pcoeff2(2)))
end
hold off
xlim([min(s_axis),max(s_axis)])
ylim([min(F_s),max(F_s)])
xlabel('Log Window Number, s')
ylabel('Log F(s)')
title(['Graph to Calculate Hurst Exponent, ',Sample_Name])
saveas(gcf,['Figures\',technique,' ',num2str(sample_no),...
'Hurst Exponent'],'png')

% ------
% Figure 2
% ------
% Plot spectogram
Fig2 = figure(2);
imagesc(spec_time, F_hard, power);
image_spec = gca;
image_spec.YDir = 'normal';
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title(['Spectogram, ',Sample_Name])
cbar = colorbar;
cbar.Label.String = 'Power/Freq (dB/Hz)';
axis tight;
%saveas(gcf,['Figures\',technique,' ',num2str(sample_no),...
%'Spectogram'],'epsc')

% ------
% Figure 3
% ------
% Plot original waveform
Fig3 = figure(3);
plot_orig = subplot(3,1,1);
hold on
plot(time,audio_data);
plot(onset_time2,onset_amp2,'kx');
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0,time(end)]);
hold off

% Plot filtered waveform
plot_FIR = subplot(3,1,2);
plot(time,audio_data_FIR);
xlabel('Time (s)');
ylabel('Amplitude');
xlim([0,time(end)]);

% Plot power
plot_power = subplot(3,1,3);
hold on
plot(time_power,power_tot);
plot(onset_time1,peak_power1,'kx');
xlabel('Time (s)');
ylabel('PSD (dB/Hz)');
xlim([0,time(end)]);
hold off
%saveas(gcf,['Figures\',technique,' ',num2str(sample_no),...
%'Waveforms'],'epsc')

% ------
% Figure 4
% ------
% Plot onset time difference
Fig4 = figure(4);
hold on
plot_delta_y = subplot(2,1,1);
xlim([0,length(delta_y)]);

% Plot windowed onset time difference fits
hold on
plot(delta_y_s,'k-','LineWidth',1)
plot(delta_y,'LineWidth',0.02)
xlabel('Onset Number')
ylabel('Drift From Mean (s)')
xlim([0,length(delta_y)]);

% Plot vertical lines of width window length
for linedif = 1:s:length(delta_y)
    ax = gca;
    line([linedif,linedif],ax.YLim,...
        'Color','r','LineStyle','--','LineWidth',0.01)
end

% Plot onset time variation
plot_onsetdif = subplot(2,1,2);
hold all
plot(exp_onset-delta_t);
xlabel('Onset Number');
ylabel('Onset Variation (s)');
title(['Onset Time Variation, ',Sample_Name])
%saveas(gcf,['Figures\',technique,' ',num2str(sample_no),' Drift'],'epsc')

%% "Fun" Stuff

% Link axes tile figure windows (Created by Peter J. Acklam, 1998)
linkaxes([plot_orig,plot_FIR,plot_power,image_spec],'x')
linkaxes([plot_delta_y,plot_onsetdif],'x')
tilefigs

% Output time /percentage of onsets found.
fprintf(['Theoretical onset difference = ',num2str(60/bpm/4),' secs. \n'])
fprintf(['Calculated onset difference = ',num2str(mean_t),' secs.\n'])
disp([num2str(length(locs1)/(Npeaksmax-4)*100),'% of peaks found'])

end

% Play
%player = audioplayer(audio_data, 1*audio_freq);
%player = audioplayer(audio_data_FIR, 1*audio_freq);
%play(player);

% Save final exponents and differences to file
%dlmwrite(['Docs\',num2str(technique),' Exponents.txt'],final_result)
%dlmwrite(['Docs\',num2str(technique),' Variation.txt'],...
%    [final_onset_difs;final_power_difs])

beep
disp('Done!')