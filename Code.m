% Load spike data
data_P = load('Pe170417_spikes.mat'); % Load data for subject P
data_W = load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for a specific unit
% Example for subject P
unit_spike_times_P = data_P.ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial
% Example for subject W
unit_spike_times_W = data_W.ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial

% Define bin size for histogram
bin_size = 0.01; % 10 ms bins

% Create spike train (histogram) for subject P
edges_P = 0:bin_size:max(unit_spike_times_P); % Define the edges of bins
spike_train_P = histcounts(unit_spike_times_P, edges_P); % Calculate the histogram

% Create spike train (histogram) for subject W
edges_W = 0:bin_size:max(unit_spike_times_W); % Define the edges of bins
spike_train_W = histcounts(unit_spike_times_W, edges_W); % Calculate the histogram

% Plot spike train histogram for subject P
figure;
bar(edges_P(1:end-1), spike_train_P, 'histc'); % Plot histogram
title('Spike Train Histogram - Subject P');
xlabel('Time (s)');
ylabel('Spike Count');

% Plot spike train histogram for subject W
figure;
bar(edges_W(1:end-1), spike_train_W, 'histc'); % Plot histogram
title('Spike Train Histogram - Subject W');
xlabel('Time (s)');
ylabel('Spike Count');


% Number of trials and spike times
num_trials = 100;
num_spikes_per_trial = 10;

% Generate random spike times
spike_times = rand(num_trials, num_spikes_per_trial) * 360;

% Create a figure for the raster plot
figure;
hold on;

% Plot spike times for each trial
for trial = 1:num_trials
    for spike = 1:num_spikes_per_trial
        line([spike_times(trial, spike), spike_times(trial, spike)], ...
             [trial - 0.4, trial + 0.4], 'Color', 'k');
    end
end

% Set plot limits and labels
xlim([0, 360]);
ylim([0, num_trials + 1]);
xlabel('Time (ms)');
ylabel('Trial');
title('Raster Plot of Spike Times');

hold off;

% Load spike data
data_P = load('Pe170417_spikes.mat'); % Load data for subject P
data_W = load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for specific units
% Example for subject P
unit_spike_times_P = data_P.ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial for subject P
% Example for subject W
unit_spike_times_W = data_W.ex.EVENTS{2, 1, 1}; % Extract spike times for the second unit, first stimulus, first trial for subject W

% Define total recording time for each subject (assume it is the max spike time)
total_time_P = max(unit_spike_times_P);
total_time_W = max(unit_spike_times_W);

% Calculate the number of spikes for each subject
num_spikes_P = length(unit_spike_times_P);
num_spikes_W = length(unit_spike_times_W);

% Calculate the firing rate for each subject
firing_rate_P = num_spikes_P / total_time_P; % Firing rate for subject P (spikes per second)
firing_rate_W = num_spikes_W / total_time_W; % Firing rate for subject W (spikes per second)

% Display the results
fprintf('Firing Rate for Subject P: %.2f spikes/second\n', firing_rate_P);
fprintf('Firing Rate for Subject W: %.2f spikes/second\n', firing_rate_W);

% Plot the firing rates
figure;
bar([firing_rate_P, firing_rate_W]);
set(gca, 'XTickLabel', {'Subject P', 'Subject W'});
title('Firing Rates of Subjects');
xlabel('Subjects');
ylabel('Firing Rate (spikes/second)');

% Load spike data
data_P = load('Pe170417_spikes.mat'); % Load data for subject P
data_W = load('Wi170428_spikes.mat'); % Load data for subject W

% Initialize variables
num_stimuli = 9; % Number of stimuli
bin_size = 0.01; % 10 ms bins

% Example for subject P
tuning_curve_P = zeros(1, num_stimuli);
for stim = 1:num_stimuli
    % Extract spike times for each stimulus
    spike_times_P = data_P.ex.EVENTS{1, stim, 1}; % First unit, different stimuli, first trial for subject P
    % Calculate firing rate for each stimulus
    firing_rate_P = length(spike_times_P) / (max(spike_times_P) - min(spike_times_P));
    tuning_curve_P(stim) = firing_rate_P;
end

% Example for subject W
tuning_curve_W = zeros(1, num_stimuli);
for stim = 1:num_stimuli
    % Extract spike times for each stimulus
    spike_times_W = data_W.ex.EVENTS{1, stim, 1}; % First unit, different stimuli, first trial for subject W
    % Calculate firing rate for each stimulus
    firing_rate_W = length(spike_times_W) / (max(spike_times_W) - min(spike_times_W));
    tuning_curve_W(stim) = firing_rate_W;
end

% Plot the tuning curves
figure;
subplot(2, 1, 1);
bar(tuning_curve_P);
title('Tuning Curve - Subject P');
xlabel('Stimulus');
ylabel('Firing Rate (spikes/second)');

subplot(2, 1, 2);
bar(tuning_curve_W);
title('Tuning Curve - Subject W');
xlabel('Stimulus');
ylabel('Firing Rate (spikes/second)');

% Load spike data
load('Pe170417_spikes.mat'); % Load data for subject P

% Number of stimuli (example: 8 directions)
num_stimuli = 8;
directions = linspace(0, 2*pi, num_stimuli + 1);
directions = directions(1:end-1); % Remove the last point to avoid duplication

% Initialize array to store firing rates
firing_rates = zeros(1, num_stimuli);

% Calculate firing rate for each stimulus
for stim = 1:num_stimuli
    % Extract spike times for each stimulus
    spike_times = ex.EVENTS{1, stim, 1}; % Example for the first unit, different stimuli
    % Calculate firing rate for each stimulus
    firing_rate = length(spike_times) / (max(spike_times) - min(spike_times));
    firing_rates(stim) = firing_rate;
end

% Plot the tuning curve as a polar plot
figure;
polarplot(directions, firing_rates, '-o');
title('Tuning Curve');

% Load spike data
load('Pe170417_spikes.mat'); % Load data for subject P

% Extract spike times for specific unit
unit_spike_times_P = ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial for subject P

% Calculate ISI for subject P
isi_P = diff(unit_spike_times_P); % ISI for subject P

% Plot the ISI histogram for subject P
figure;
histogram(isi_P, 'BinWidth', 0.01); % Histogram for ISI of subject P
title('ISI Histogram - Subject P');
xlabel('Inter-Spike Interval (s)');
ylabel('Frequency');

% Load spike data
load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for specific unit
unit_spike_times_W = ex.EVENTS{2, 1, 1}; % Extract spike times for the second unit, first stimulus, first trial for subject W

% Calculate ISI for subject W
isi_W = diff(unit_spike_times_W); % ISI for subject W

% Plot the ISI histogram for subject W
figure;
histogram(isi_W, 'BinWidth', 0.01); % Histogram for ISI of subject W
title('ISI Histogram - Subject W');
xlabel('Inter-Spike Interval (s)');
ylabel('Frequency');

% Load spike data
load('Pe170417_spikes.mat'); % Load data for subject P

% Extract spike times for specific unit
unit_spike_times_P = ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial for subject P

% Calculate ISI for subject P
isi_P = diff(unit_spike_times_P); % ISI for subject P

% Calculate CV for subject P
cv_P = std(isi_P) / mean(isi_P); % CV for subject P

% Display the CV result
fprintf('CV for Subject P: %.2f\n', cv_P);

% Plot the ISI histogram for subject P and display CV
figure;
histogram(isi_P, 'BinWidth', 0.01); % Histogram for ISI of subject P
title(['ISI Histogram - Subject P (CV = ', num2str(cv_P, '%.2f'), ')']);
xlabel('Inter-Spike Interval (s)');
ylabel('Frequency');


% Load spike data
load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for specific unit
unit_spike_times_W = ex.EVENTS{2, 1, 1}; % Extract spike times for the second unit, first stimulus, first trial for subject W

% Calculate ISI for subject W
isi_W = diff(unit_spike_times_W); % ISI for subject W

% Calculate CV for subject W
cv_W = std(isi_W) / mean(isi_W); % CV for subject W

% Display the CV result
fprintf('CV for Subject W: %.2f\n', cv_W);

% Plot the ISI histogram for subject W and display CV
figure;
histogram(isi_W, 'BinWidth', 0.01); % Histogram for ISI of subject W
title(['ISI Histogram - Subject W (CV = ', num2str(cv_W, '%.2f'), ')']);
xlabel('Inter-Spike Interval (s)');
ylabel('Frequency');


% Load spike data
load('Pe170417_spikes.mat'); % Load data for subject P

% Extract spike times for specific unit
unit_spike_times_P = ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial for subject P

% Define time periods (in seconds)
pre_decision_period = [0, 1]; % Example: 0 to 1 second before decision
during_decision_period = [1, 2]; % Example: 1 to 2 seconds during decision
post_decision_period = [2, 3]; % Example: 2 to 3 seconds after decision

% Extract spike times for each period
spikes_pre_decision_P = unit_spike_times_P(unit_spike_times_P >= pre_decision_period(1) & unit_spike_times_P < pre_decision_period(2));
spikes_during_decision_P = unit_spike_times_P(unit_spike_times_P >= during_decision_period(1) & unit_spike_times_P < during_decision_period(2));
spikes_post_decision_P = unit_spike_times_P(unit_spike_times_P >= post_decision_period(1) & unit_spike_times_P < post_decision_period(2));

% Calculate firing rates for each period
firing_rate_pre_P = length(spikes_pre_decision_P) / diff(pre_decision_period);
firing_rate_during_P = length(spikes_during_decision_P) / diff(during_decision_period);
firing_rate_post_P = length(spikes_post_decision_P) / diff(post_decision_period);

% Display the results
fprintf('Firing Rate for Subject P (Pre-decision): %.2f spikes/second\n', firing_rate_pre_P);
fprintf('Firing Rate for Subject P (During-decision): %.2f spikes/second\n', firing_rate_during_P);
fprintf('Firing Rate for Subject P (Post-decision): %.2f spikes/second\n', firing_rate_post_P);

% Plot the firing rates
figure;
bar([firing_rate_pre_P, firing_rate_during_P, firing_rate_post_P]);
set(gca, 'XTickLabel', {'Pre-decision', 'During-decision', 'Post-decision'});
title('Firing Rates in Different Periods - Subject P');
xlabel('Periods');
ylabel('Firing Rate (spikes/second)');

% Load spike data
load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for specific unit
unit_spike_times_W = ex.EVENTS{2, 1, 1}; % Extract spike times for the second unit, first stimulus, first trial for subject W

% Define time periods (in seconds)
pre_decision_period = [0, 1]; % Example: 0 to 1 second before decision
during_decision_period = [1, 2]; % Example: 1 to 2 seconds during decision
post_decision_period = [2, 3]; % Example: 2 to 3 seconds after decision

% Extract spike times for each period
spikes_pre_decision_W = unit_spike_times_W(unit_spike_times_W >= pre_decision_period(1) & unit_spike_times_W < pre_decision_period(2));
spikes_during_decision_W = unit_spike_times_W(unit_spike_times_W >= during_decision_period(1) & unit_spike_times_W < during_decision_period(2));
spikes_post_decision_W = unit_spike_times_W(unit_spike_times_W >= post_decision_period(1) & unit_spike_times_W < post_decision_period(2));

% Calculate firing rates for each period
firing_rate_pre_W = length(spikes_pre_decision_W) / diff(pre_decision_period);
firing_rate_during_W = length(spikes_during_decision_W) / diff(during_decision_period);
firing_rate_post_W = length(spikes_post_decision_W) / diff(post_decision_period);

% Display the results
fprintf('Firing Rate for Subject W (Pre-decision): %.2f spikes/second\n', firing_rate_pre_W);
fprintf('Firing Rate for Subject W (During-decision): %.2f spikes/second\n', firing_rate_during_W);
fprintf('Firing Rate for Subject W (Post-decision): %.2f spikes/second\n', firing_rate_post_W);

% Plot the firing rates
figure;
bar([firing_rate_pre_W, firing_rate_during_W, firing_rate_post_W]);
set(gca, 'XTickLabel', {'Pre-decision', 'During-decision', 'Post-decision'});
title('Firing Rates in Different Periods - Subject W');
xlabel('Periods');
ylabel('Firing Rate (spikes/second)');

% Load spike data
load('Pe170417_spikes.mat'); % Load data for subject P

% Extract spike times for specific unit
unit_spike_times_P = ex.EVENTS{1, 1, 1}; % Extract spike times for the first unit, first stimulus, first trial for subject P

% Define time periods (in seconds)
pre_decision_period = [0, 1]; % Example: 0 to 1 second before decision
during_decision_period = [1, 2]; % Example: 1 to 2 seconds during decision
post_decision_period = [2, 3]; % Example: 2 to 3 seconds after decision

% Extract spike times for each period
spikes_pre_decision_P = unit_spike_times_P(unit_spike_times_P >= pre_decision_period(1) & unit_spike_times_P < pre_decision_period(2));
spikes_during_decision_P = unit_spike_times_P(unit_spike_times_P >= during_decision_period(1) & unit_spike_times_P < during_decision_period(2));
spikes_post_decision_P = unit_spike_times_P(unit_spike_times_P >= post_decision_period(1) & unit_spike_times_P < post_decision_period(2));

% Calculate firing rates for each period
firing_rate_pre_P = length(spikes_pre_decision_P) / diff(pre_decision_period); % Firing rate before decision
firing_rate_during_P = length(spikes_during_decision_P) / diff(during_decision_period); % Firing rate during decision
firing_rate_post_P = length(spikes_post_decision_P) / diff(post_decision_period); % Firing rate after decision

% Display the results
fprintf('Firing Rate for Subject P (Pre-decision): %.2f spikes/second\n', firing_rate_pre_P);
fprintf('Firing Rate for Subject P (During-decision): %.2f spikes/second\n', firing_rate_during_P);
fprintf('Firing Rate for Subject P (Post-decision): %.2f spikes/second\n', firing_rate_post_P);

% Statistical comparison using t-test
[h_pre_during_P, p_pre_during_P] = ttest2(spikes_pre_decision_P, spikes_during_decision_P);
[h_during_post_P, p_during_post_P] = ttest2(spikes_during_decision_P, spikes_post_decision_P);

% Display the statistical results
fprintf('T-test (Pre-decision vs. During-decision) for Subject P: h = %d, p = %.3f\n', h_pre_during_P, p_pre_during_P);
fprintf('T-test (During-decision vs. Post-decision) for Subject P: h = %d, p = %.3f\n', h_during_post_P, p_during_post_P);

% Plot the firing rates
figure;
bar([firing_rate_pre_P, firing_rate_during_P, firing_rate_post_P]);
set(gca, 'XTickLabel', {'Pre-decision', 'During-decision', 'Post-decision'});
title('Firing Rates in Different Periods - Subject P');
xlabel('Periods');
ylabel('Firing Rate (spikes/second)');

% Load spike data
load('Wi170428_spikes.mat'); % Load data for subject W

% Extract spike times for specific unit
unit_spike_times_W = ex.EVENTS{2, 1, 1}; % Extract spike times for the second unit, first stimulus, first trial for subject W

% Define time periods (in seconds)
pre_decision_period = [0, 1]; % Example: 0 to 1 second before decision
during_decision_period = [1, 2]; % Example: 1 to 2 seconds during decision
post_decision_period = [2, 3]; % Example: 2 to 3 seconds after decision

% Extract spike times for each period
spikes_pre_decision_W = unit_spike_times_W(unit_spike_times_W >= pre_decision_period(1) & unit_spike_times_W < pre_decision_period(2));
spikes_during_decision_W = unit_spike_times_W(unit_spike_times_W >= during_decision_period(1) & unit_spike_times_W < during_decision_period(2));
spikes_post_decision_W = unit_spike_times_W(unit_spike_times_W >= post_decision_period(1) & unit_spike_times_W < post_decision_period(2));

% Calculate firing rates for each period
firing_rate_pre_W = length(spikes_pre_decision_W) / diff(pre_decision_period); % Firing rate before decision
firing_rate_during_W = length(spikes_during_decision_W) / diff(during_decision_period); % Firing rate during decision
firing_rate_post_W = length(spikes_post_decision_W) / diff(post_decision_period); % Firing rate after decision

% Display the results
fprintf('Firing Rate for Subject W (Pre-decision): %.2f spikes/second\n', firing_rate_pre_W);
fprintf('Firing Rate for Subject W (During-decision): %.2f spikes/second\n', firing_rate_during_W);
fprintf('Firing Rate for Subject W (Post-decision): %.2f spikes/second\n', firing_rate_post_W);

% Statistical comparison using t-test
[h_pre_during_W, p_pre_during_W] = ttest2(spikes_pre_decision_W, spikes_during_decision_W);
[h_during_post_W, p_during_post_W] = ttest2(spikes_during_decision_W, spikes_post_decision_W);

% Display the statistical results
fprintf('T-test (Pre-decision vs. During-decision) for Subject W: h = %d, p = %.3f\n', h_pre_during_W, p_pre_during_W);
fprintf('T-test (During-decision vs. Post-decision) for Subject W: h = %d, p = %.3f\n', h_during_post_W, p_during_post_W);

% Plot the firing rates
figure;
bar([firing_rate_pre_W, firing_rate_during_W, firing_rate_post_W]);
set(gca, 'XTickLabel', {'Pre-decision', 'During-decision', 'Post-decision'});
title('Firing Rates in Different Periods - Subject W');
xlabel('Periods');
ylabel('Firing Rate (spikes/second)');
