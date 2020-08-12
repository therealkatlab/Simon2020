% -------------------------------------------------------------------------
%   This script accompanies the manuscript                                 
%   Simon et al., (2020) Developmental Cell                                
%   Repository available on https://github.com/therealkatlab               
%   Please consult READ_ME for more information                            
% -------------------------------------------------------------------------
% 
% 
% -------------------------------------------------------------------------
%   Optimal threshold finder to classify MEKi (low ERK) vs FGF (high ERK)               
% -------------------------------------------------------------------------

%% Load data that contains ERK-KTR C:N values for treated PD03 and FGF4 cells
% PD03 = MEKi (low ERK) treated ICM cells
% FGF4 = FGF (high ERK) treated ICM cells
% Do not include control (untreated) ICM cells

Xm = transpose(table2array(PD03(:, 2)));
Xf = transpose(table2array(FGF4(:, 2)));
x = [Xm Xf];
m = length(Xm) ;
f = length(Xf) ;
y = [zeros(1,m) ones(1,f)];

dt = 0.01; % the step among thresholds to be tested
t = min(x)+dt:dt:max(x)-dt; % the thresholds that will be tested
N = length(x); % the total number of cells

tpr = zeros(length(t),1); % to collect the true positive rate for each threshold 
fpr = zeros(length(t),1); % to collect the false positive rate for each threshold 
acc = zeros(length(t),1); % to collect the accuracy of each threshold as a classifier

num_pos = sum(y); % number of FGF-treated cells
num_neg = N - num_pos; % number of MEKi-treated cells

for i = 1:length(t) % looping through all thresholds
    class_pos = x >= t(i); % which cells are classified as high ERK with this threshold?
    
    num_true_pos  = sum(y .* class_pos); % number of cells correctly classified as high ERK
    num_false_pos = sum(class_pos) - num_true_pos; % number of cells incorrectly classified as high ERK
    
    tpr(i) = num_true_pos  / num_pos; % true positive rate (number of cells correctly classified as high ERK)/(number of FGF-treated cells)
    fpr(i) = num_false_pos / num_neg; % false positive rate (number of cells incorrectly classified as high ERK)/(number of MEKi-treated cells)
    
    acc(i) = sum(y == class_pos) / N; % accuracy (number of cells correctly classified)/(total number of cells)
end

best_thresholds = t(acc == max(acc)); % find the most accurate threshold(s)
if mod(length(best_thresholds),2) == 1
    best_thresh = median(best_thresholds); % if there are an odd number of maximally accurate thresholds, choose the median
else 
    best_thresh = median([best_thresholds(1) best_thresholds]); % if there are an even number of maximally accurate thresholds, choose the one just lower than the median
end

best_thresh

figure(1)
hold on
plot(x, y, '+')
plot([best_thresh best_thresh], [0 1], 'r--')
axis([min(x) max(x) 0 1])
xlabel 'ERK activity'
ylabel 'treatment'

figure(2)
hold on
plot(fpr, tpr)
plot(fpr(t == best_thresh), tpr(t == best_thresh), 'or')
axis([0 1 0 1])
axis square
title 'roc curve'
xlabel 'false positive rate'
ylabel 'true positive rate'

figure(3)
hold on
plot(t, acc)
plot(t(t == best_thresh), acc(t == best_thresh), 'or')
axis([min(x) max(x) 0 1])
title 'accuracy of different thresholds'
xlabel 'ERK activity threshold'
ylabel 'accuracy'

