% Oliver Gordon & Dominic Coy, 2017
% Is Drumming Fractal? - Analysis Code

clear variables
close all

%% Importing
technique_all = {'Double Handed','Single Handed','Metronomic'};
no_samples = 22;
no_techniques = 2;

% Preallocate storage arrays
variation_all = NaN(no_samples*2,3000*no_techniques);
exponents_all = NaN(7,no_samples*no_techniques);

% Read in stored values
experience_all = xlsread('Docs\People.xlsx');
for technique = 1:2           
    
    % Read in experiences
    experience = repmat(experience_all(:,4),no_techniques,1);
    
    % Read in variations and separate out
    variation_all(:,((((technique-1)*3000)+1):(technique*3000))) = ...
        dlmread(['Docs\',char(technique_all(technique)),' Variation.txt']);
    variation_onset = variation_all(1:no_samples,:);
    variation_power = variation_all(no_samples+1:end,:);
    
    % Read in exponents and separate out
    exponents_all(:,((((technique-1)*no_samples)+1):...
        (technique*no_samples))) = ...
        dlmread(['Docs\',char(technique_all(technique)),' Exponents.txt']);
    exponents_class = exponents_all(1,:);
    exponents_double = exponents_all([2,4],:);
    exponents_double_err = sqrt((exponents_all([3,5],:).^2)+(0.02^2));
    exponents_single = exponents_all(6,:);
    exponents_single_err = sqrt((exponents_all(7,:)).^2+(0.02^2));
    
end
%% Analysis

% Find index of values that are forward db, reverse db, straight
index_forward = find(exponents_class == 11);
index_backward = find(exponents_class == 12);
index_straight = find(exponents_class == 13);

% Find index of values in each experience range
index_0_5 = find(experience<=5);
index_5_10 = find(experience>5 & experience<=10);
index_10_15 = find(experience>10 & experience<=15);
index_15_max = find(experience>15);

% Find percentage of results forward, reverse, straight for experience
% ranges
perc_forward_0_5 = sum(sum(index_0_5 == index_forward))/...
    (no_samples.*no_techniques)*100;
perc_forward_5_10 = sum(sum(index_5_10 == index_forward))/...
    (no_samples.*no_techniques)*100;
perc_forward_10_15 = sum(sum(index_10_15 == index_forward))/...
    (no_samples.*no_techniques)*100;
perc_forward_15_end = sum(sum(index_15_max == index_forward))/...
    (no_samples.*no_techniques)*100;

perc_backward_0_5 = sum(sum(index_0_5 == index_backward))/...
    (no_samples.*no_techniques)*100;
perc_backward_5_10 = sum(sum(index_5_10 == index_backward))/...
    (no_samples.*no_techniques)*100;
perc_backward_10_15 = sum(sum(index_10_15 == index_backward))/...
    (no_samples.*no_techniques)*100;
perc_backward_15_end = sum(sum(index_15_max == index_backward))/...
    (no_samples.*no_techniques)*100;

perc_straight_0_5 = sum(sum(index_0_5 == index_straight))/...
    (no_samples.*no_techniques)*100;
perc_straight_5_10 = sum(sum(index_5_10 == index_straight))/...
    (no_samples.*no_techniques)*100;
perc_straight_10_15 = sum(sum(index_10_15 == index_straight))/...
    (no_samples.*no_techniques)*100;
perc_straight_15_end = sum(sum(index_15_max == index_straight))/...
    (no_samples.*no_techniques)*100;

% Store this as matrix
perc_age = [perc_straight_0_5, perc_forward_0_5, perc_backward_0_5 ; ...
    perc_straight_5_10, perc_forward_5_10, perc_backward_5_10 ; ...
    perc_straight_10_15, perc_forward_10_15, perc_backward_10_15 ; ...
    perc_straight_15_end, perc_forward_15_end, perc_backward_15_end];

% Find percentage of results that are forward, reverse, straight
perc_straight_tot = sum(perc_age(:,1));
perc_forward_tot = sum(perc_age(:,2));
perc_backward_tot = sum(perc_age(:,3));

% Find percentage of results outside 1 s.d. of 0.5 for straight 
mean = 0.5;
ins_sing = exponents_single((exponents_single > ...
    mean-exponents_single_err) & ...
    (exponents_single < mean+exponents_single_err));
index_ins_sing = find(ismember(exponents_single,ins_sing));

% Find percentage of straight 
% results outside 1 s.d of 0.5, for different experiences
ins_0_5 = length(find(ismember(index_ins_sing,index_0_5)))...
    /length(index_0_5)*100;
ins_5_10 = length(find(ismember(index_ins_sing,index_5_10)))...
    /length(index_5_10)*100;
ins_10_15 = length(find(ismember(index_ins_sing,index_10_15)))...
    /length(index_10_15)*100;
ins_15_end = length(find(ismember(index_ins_sing,index_15_max)))...
    /length(index_15_max)*100;

% Find total percentage of straight results outside 1 s.d. 
disp([num2str((length(ins_sing))/...
    no_samples/no_techniques*100), '% of results within 1 s.d of '...,
    num2str(mean)]);
disp([num2str(ins_0_5),'% of 0-5 years within 1 s.d. of ',...
    num2str(mean)])
disp([num2str(ins_5_10),'% of 5-10 years within 1 s.d. of ',...
    num2str(mean)])
disp([num2str(ins_10_15),'% of 10-15 years within 1 s.d. of ',...
    num2str(mean)])
disp([num2str(ins_15_end),'% of 15+ years within 1 s.d. of ',...
    num2str(mean)])

% Find standard deviation of straight results for different experiences
disp(['std of 0-5 group = ',...
    num2str(nanstd(exponents_single(index_0_5)))]);
disp(['std of 5-10 group = ',...
    num2str(nanstd(exponents_single(index_5_10)))]);
disp(['std of 10-15 group = ',...
    num2str(nanstd(exponents_single(index_10_15)))]);
disp(['std of 15+ group = ',...
    num2str(nanstd(exponents_single(index_15_max)))]);

% Find means of onsets and power for both axis
mean_variation_onset = nanmean(variation_onset);
mean_variation_power = nanmean(variation_power);
row_mean_variation_onset = nanmean(variation_onset,2);
row_mean_variation_power = nanmean(variation_power,2);

% Find drift away from mean drift for both power and onset
bpm = 88;
exp_onset = 60/bpm/4;

drift_mean_onset = mean_variation_onset - exp_onset;
drift_mean_power = mean_variation_power - exp_onset;

%% Plotting

% ------
% Figure 1
% ------
% Plot drift from mean onset time
figure(1)
plot(cumsum(drift_mean_onset(~isnan(drift_mean_onset)))); 
xlabel('Onset Number');
ylabel('Drift From Mean (s)');
title('Graph to Show the Drift Away From the Mean for Onset time');

% ------
% Figure 2
% ------
% Plot drift from mean onset power
figure(2)
plot(cumsum(drift_mean_power(~isnan(drift_mean_power)))); 
xlabel('Onset Number');
ylabel('Drift Away From the Mean');
title('Graph to Show the Drift Away From the Mean for Onset Power');

% ------
% Figure 3
% ------
% Plot histogram of hurst exponent, with different subfigures for different
% experience ranges
figure(3);
max_val = max(experience);
experience_onset = [experience, exponents_all(6, :)'];
filterLims(experience_onset, 0, 5, 1);
filterLims(experience_onset, 6, 10, 2) ;
filterLims(experience_onset, 11, 15, 3);
filterLims(experience_onset, 16, max_val, 4) ;

% -------------
% Figures 4 - 8
% -------------
nbins = 15;
sample_no = 1:14;

% Histogram of hurst exponent
figure(4)
histogram(exponents_single)
xlabel('Hurst Exponent');
ylabel('Number');
title('Histogram to Show Hurst Exponents'); 
disp(['std = ', num2str(nanstd(exponents_single))])
disp(['range = ', num2str(max(exponents_single)-min(exponents_single))])
disp(['mean = ', num2str(nanmean(exponents_single))])

% Histogram of all onset time variation
figure(5)
histogram(variation_onset(sample_no,:));
xlabel('Onset Variation (s)');
ylabel('Number');
title('Histogram to Show Onset Time Variations');

% Histogram of all onset power variation
figure(6)
histogram(variation_power(sample_no,:));
xlabel('Power Variation');
ylabel('Number');
title('Histogram to Show Power Variations');
xlim([min(min(variation_power(sample_no,:))), ...
    max(max(variation_power(sample_no,:)))])

% Histogram of mean onset time variation
figure(7)
histogram(row_mean_variation_onset,nbins);
xlabel('Mean Onset Variation (s)');
ylabel('Number');
title('Histogram to Show the Mean Onset Variations');

% Histogram of mean power variation
figure(8)
histogram(row_mean_variation_power,nbins);
xlabel('Mean Power Variation');
ylabel('Number');
title('Histogram to Show the Mean Power Variations');