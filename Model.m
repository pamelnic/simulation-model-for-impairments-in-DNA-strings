
% Monte Carlo simulation settings

numTrials = 10000;   %  trials

% Reference values
ref_del = 0.0295;
ref_ins = 0.0216;
ref_sub = 0.0175;

% Preallocate arrays
all_del = zeros(numTrials,1);
all_ins = zeros(numTrials,1);
all_sub = zeros(numTrials,1);


% Monte Carlo loop

for t = 1:numTrials
    [summary_del, summary_ins, summary_sub] = iid_channel();

    all_del(t) = summary_del;
    all_ins(t) = summary_ins;
    all_sub(t) = summary_sub;
end

% Store raw results

results.all_del = all_del;
results.all_ins = all_ins;
results.all_sub = all_sub;


% Statistical evaluation

% Means
results.mean_del = mean(all_del);
results.mean_ins = mean(all_ins);
results.mean_sub = mean(all_sub);

% Variances
results.var_del = var(all_del);
results.var_ins = var(all_ins);
results.var_sub = var(all_sub);

% Standard deviations
results.std_del = std(all_del);
results.std_ins = std(all_ins);
results.std_sub = std(all_sub);

% Absolute error
results.ae_del = abs(results.mean_del - ref_del);
results.ae_ins = abs(results.mean_ins - ref_ins);
results.ae_sub = abs(results.mean_sub - ref_sub);

% Relative error in %
results.re_del = ((results.mean_del - ref_del) / ref_del) * 100;
results.re_ins = ((results.mean_ins - ref_ins) / ref_ins) * 100;
results.re_sub = ((results.mean_sub - ref_sub) / ref_sub) * 100;

% 95% confidence intervals
alpha = 0.05;
tcrit = tinv(1 - alpha/2, numTrials - 1);

results.ci_half_del = tcrit * results.std_del / sqrt(numTrials);
results.ci_half_ins = tcrit * results.std_ins / sqrt(numTrials);
results.ci_half_sub = tcrit * results.std_sub / sqrt(numTrials);

results.ci_del = [results.mean_del - results.ci_half_del, results.mean_del + results.ci_half_del];
results.ci_ins = [results.mean_ins - results.ci_half_ins, results.mean_ins + results.ci_half_ins];
results.ci_sub = [results.mean_sub - results.ci_half_sub, results.mean_sub + results.ci_half_sub];

% One-sample t-test statistic
results.tstat_del = (results.mean_del - ref_del) / (results.std_del / sqrt(numTrials));
results.tstat_ins = (results.mean_ins - ref_ins) / (results.std_ins / sqrt(numTrials));
results.tstat_sub = (results.mean_sub - ref_sub) / (results.std_sub / sqrt(numTrials));

% Command window output

fprintf('\n');
fprintf('===============================================================\n');
fprintf('        STATISTICAL EVALUATION OF MONTE CARLO SIMULATION       \n');
fprintf('===============================================================\n');
fprintf('Number of trials: %d\n\n', numTrials);

fprintf('Reference values:\n');
fprintf('  Deletions      : %.6f\n', ref_del);
fprintf('  Insertions     : %.6f\n', ref_ins);
fprintf('  Substitutions  : %.6f\n\n', ref_sub);

fprintf('---------------------------------------------------------------\n');
fprintf('%-15s %-12s %-12s %-12s\n', 'Metric', 'Deletions', 'Insertions', 'Substitutions');
fprintf('---------------------------------------------------------------\n');
fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 'Mean', ...
    results.mean_del, results.mean_ins, results.mean_sub);
fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 'Variance', ...
    results.var_del, results.var_ins, results.var_sub);
fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 'Std. dev.', ...
    results.std_del, results.std_ins, results.std_sub);

fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 'Abs. error', ...
    results.ae_del, results.ae_ins, results.ae_sub);
fprintf('%-15s %-11.2f%% %-11.2f%% %-11.2f%%\n', 'Rel. error', ...
    results.re_del, results.re_ins, results.re_sub);
fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 't-statistic', ...
    results.tstat_del, results.tstat_ins, results.tstat_sub);
fprintf('%-15s %-12.6f %-12.6f %-12.6f\n', 'p-value', ...
    results.pval_del, results.pval_ins, results.pval_sub);
fprintf('---------------------------------------------------------------\n\n');

fprintf('95%% Confidence intervals for the mean:\n');
fprintf('  Deletions      : [%.6f, %.6f]\n', results.ci_del(1), results.ci_del(2));
fprintf('  Insertions     : [%.6f, %.6f]\n', results.ci_ins(1), results.ci_ins(2));
fprintf('  Substitutions  : [%.6f, %.6f]\n\n', results.ci_sub(1), results.ci_sub(2));


fprintf('===============================================================\n\n');

% Plot 1: Histogram of deletions

figure;
h = histogram(all_del);
xlabel('Deletion value');
ylabel('Frequency');
title('Distribution of deletions per trial');
grid on;

binCounts = h.BinCounts;
binEdges = h.BinEdges;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
ylim([0, max(binCounts)*1.25]);


hold on;
for i = 1:length(binCounts)
    if binCounts(i) > 0
        text(binCenters(i), binCounts(i), num2str(binCounts(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9, ...
            'FontWeight', 'bold');
    end
end
hold off;


% Plot 2: Histogram of insertions

figure;
h = histogram(all_ins);
xlabel('Insertion value');
ylabel('Frequency');
title('Distribution of insertions per trial');
grid on;

binCounts = h.BinCounts;
binEdges = h.BinEdges;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
ylim([0, max(binCounts)*1.25]);

hold on;
for i = 1:length(binCounts)
    if binCounts(i) > 0
        text(binCenters(i), binCounts(i), num2str(binCounts(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9, ...
            'FontWeight', 'bold');
    end
end
hold off;


% Plot 3: Histogram of substitutions

figure;
h = histogram(all_sub);
xlabel('Substitution value');
ylabel('Frequency');
title('Distribution of substitutions per trial');
grid on;

binCounts = h.BinCounts;
binEdges = h.BinEdges;
binCenters = binEdges(1:end-1) + diff(binEdges)/2;
ylim([0, max(binCounts)*1.25]);

hold on;
for i = 1:length(binCounts)
    if binCounts(i) > 0
        text(binCenters(i), binCounts(i), num2str(binCounts(i)), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', ...
            'FontSize', 9, ...
            'FontWeight', 'bold');
    end
end
hold off;
