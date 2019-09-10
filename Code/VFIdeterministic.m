close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A_h = 1; A_l = 678;
pi_hh = 0.977; pi_hl = 1-pi_hh;
pi_ll = 0.926; pi_lh = 1-pi_ll;
A=[A_h;A_l];

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons_h = A_h * k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; %consumption at A_h
cons_l = A_l * k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; %consumption at A_l

ret_h = cons_h .^ (1 - sigma) / (1 - sigma); % return function at A_h
ret_l = cons_l .^ (1 - sigma) / (1 - sigma); % return function at A_l


% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret_h(cons_h < 0) = -Inf;
ret_l(cons_l < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat_h = ret_h + beta * repmat(v_guess, [num_k 1]);%%%%%need to modify
    value_mat_l = ret_l + beta * repmat(v_guess, [num_k 1]);%%%%%need to modify

    % find the optimal k' for every k:
    [vfn_h, pol_indx] = max(value_mat_h, [], 2);
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

g = k(pol_indx); % policy function

plot(k,vfn)
figure
plot(k,g)


