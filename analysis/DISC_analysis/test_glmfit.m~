% this data is pulled from S3_1_revampHORRM, condition 1

g = [-1.2000
    -1.0588
    -0.9176
    -0.7765
    -0.6353
    -0.4941
    -0.3529
    -0.2118
    -0.0706
    0.0706
    0.2118
    0.3529
    0.4941
    0.6353
    0.7765
    0.9176
    1.0588
    1.2000]; % growth list

r = [0
    0
    0
    0
    0
    0
    0
    2
    3
    8
    16
    14
    7
    6
    3
    3
    3
    3]; % number of times subject responded "growing"

n = [     3
    3
    3
    3
    3
    3
    5
    8
    14
    26
    31
    19
    10
    6
    3
    3
    3
    4]; % total number of trials run at each growth step


% ryan's way (raw counts)
beta_from_counts = glmfit(g,[r n],'binomial', 'logit')

% chris' way (proportions)
beta_from_prop = glmfit(g,[r./n n./n],'binomial', 'logit')

% maybe it is because 'probit' and 'logit' account for the amount of data in each condition?
beta_from_eq_counts = glmfit(g,[(r./n)*max(n) repmat(max(n),size(n))],'binomial', 'logit')



% and again with data from matlab documentation:
x = [2100 2300 2500 2700 2900 3100 3300 3500 3700 3900 4100 4300]';
n = [48 42 31 34 31 21 23 23 21 16 17 21]';
y = [1 2 0 3 8 8 14 17 19 15 17 21]';
b_count = glmfit(x,[y n],'binomial','link','logit')
b_prop  = glmfit(x,[y./n n./n],'binomial','link','logit')






