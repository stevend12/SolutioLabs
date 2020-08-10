% This script plots data from the gamma_1d.cpp test program
clear;

% Always assume the data is in the "Build" folder on the same level as the
% "Plot" directory.
file_name = '../Build/gamma_1d.txt';

% Read in data
data = dlmread(file_name, ' ');

% Experiment with different threshold settings
thresholds = [0 0.02 0.04 0.06 0.08 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.6 0.7 0.8 0.9];
max_dose = max(data(:,3));

for n = 1:length(thresholds)
  threshold_ind = find(data(:,3) > thresholds(n)*max_dose);
  pr_global(n) = 100.0 * sum(data(threshold_ind,6) <= 1)/ length(threshold_ind);
  pr_local(n) = 100.0 * sum(data(threshold_ind,7) <= 1)/ length(threshold_ind);
end

% Plot results
line_width = 3;
back_color = [0.15 0.15 0.15];

% Gamma vs. Dose-Difference vs. DTA
figure('position', [100 100 1000 800])
hold on;
plot(data(:,1), data(:,2), "r", "LineWidth", line_width)
plot(data(:,1), data(:,3), "b", "LineWidth", line_width)
plot(data(:,1), data(:,4), "w", "LineWidth", line_width)
plot(data(:,1), data(:,5), "m", "LineWidth", line_width)
plot(data(:,1), data(:,6), "g", "LineWidth", line_width)
hold off;
axis([0 10 0 1.25]) % Display only half-profiles
xlabel('Position (cm)', 'color', 'w')
ylabel('Relative Dose / DTA (cm) / Gamma', 'color', 'w');
set(gca, "fontsize", 20);
h = legend("Test","Reference","Dose Diff.","DTA","Gamma");
set(h, 'fontsize', 24);
% Invert remaining colors
set(h, 'color', back_color)
set(h, 'textcolor', 'w')
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)

% Gamma parameter analysis
figure('position', [120 120 1000 800])
hold on;
plot(data(:,1), data(:,2), "--m", "LineWidth", line_width)
plot(data(:,1), data(:,6), "g", "LineWidth", line_width)
plot(data(:,1), data(:,7), "r", "LineWidth", line_width)
plot(data(:,1), data(:,8), "b", "LineWidth", line_width)
plot(data(:,1), data(:,9), "c", "LineWidth", line_width)
hold off;
axis([0 10 0 2]) % Display only half-profiles
xlabel('Position (cm)', 'color', 'w');
ylabel('Gamma', 'color', 'w');
set(gca, "fontsize", 20);
h = legend({"Test Profile","Initial", "Local Max", "2 %/2 mm", ...
  "Resample (1x)"});
set(h, 'fontsize', 24);
% Invert remaining colors
set(h, 'color', back_color)
set(h, 'textcolor', 'w')
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)

% Gamma threshold analysis
figure('position', [140 140 1000 800])
hold on;
plot(thresholds, pr_global, "g", "LineWidth", line_width)
plot(thresholds, pr_local, "b", "LineWidth", line_width)
hold off;
%axis([0 1 0 1])
xlabel('Threshold', 'color', 'w')
ylabel('Gamma Pass Rate', 'color', 'w');
set(gca, "fontsize", 20);
h = legend({"Global Max", "Local Max"});
set(h, 'fontsize', 24);
% Invert remaining colors
set(h, 'color', back_color)
set(h, 'textcolor', 'w')
set(gca, 'color', back_color)
set(gca, 'xcolor', 'w')
set(gca, 'ycolor', 'w')
set(gcf, 'color', back_color)