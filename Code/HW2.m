%$$$ HW2 Q2 Q3 %%%%
close all

%%%% set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A_h = 1.1; A_l = 0.678;
pi_hh = 0.977; pi_hl = 1-pi_hh;
pi_ll = 0.926; pi_lh = 1-pi_ll;

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
v_guess_h = zeros(1, num_k);
v_guess_l = zeros(1, num_k);

while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat_h = ret_h + beta * (pi_hh * repmat(v_guess_h, [num_k 1]) + pi_hl * repmat(v_guess_l, [num_k 1]));
    value_mat_l = ret_l + beta * (pi_lh * repmat(v_guess_h, [num_k 1]) + pi_ll * repmat(v_guess_l, [num_k 1]));

    % find the optimal k' for every k:
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';

    % what is the distance between current guess and value function
    dis = [max(abs(vfn_h - v_guess_h)); max(abs(vfn_l - v_guess_l))];
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess_h = vfn_h;
    v_guess_l = vfn_l;
end

g_h = k(pol_indx_h); % policy function at A_h
g_l = k(pol_indx_l); % policy function at A_l

%%%% Plot VF over K at each stage A %%%%

figure (1);
plot(k,vfn_h,'DisplayName','A_h', 'LineWidth',2);
hold on;
plot(k,vfn_l,'DisplayName','A_l','LineWidth',2);
hold off;
xlabel('K');
ylabel('V(K)');
title('Value Function over K at stage A','FontSize',18);
lgd = legend;
lgd.Location ='Southeast';
lgd.FontSize = 14;

%%%% Plot PF over K at each stage A %%%%

figure (2);
plot(k,g_h,'DisplayName','A_h', 'LineWidth',2);
hold on;
plot(k,g_l,'DisplayName','A_l','LineWidth',2);
hold off;
xlabel('K');
ylabel('g(K)');
title('Policy Function over K at stage A','FontSize',18);
lgd = legend;
lgd.Location ='Southeast';
lgd.FontSize = 14;

%%%% Plot savings over K at each stage A %%%%

S_h = g_h - (1-delta) * k;
S_l = g_l - (1-delta) * k;

figure (3);
plot(k,S_h,'DisplayName','A_h', 'LineWidth',2);
hold on;
plot(k,S_l,'DisplayName','A_l','LineWidth',2);
hold off;
xlabel('K');
ylabel('Saving');
title('Savings over K at stage A','FontSize',18);
lgd = legend;
lgd.Location ='Southeast';
lgd.FontSize = 14;

%%%% The end of Q2 Q3 %%%%
