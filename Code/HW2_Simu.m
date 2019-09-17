close all;
%%%% set up parameter
pi_hh = 0.977; pi_hl = 1 - pi_hh;
pi_ll = 0.926; pi_lh = 1- pi_ll;
pi = [pi_hh,pi_hl;pi_lh,pi_ll];
pi_lr = pi^10000; %long run transition probability (Markov Chain)

%%%% Requirement 1, longrun mean of A = 1 %%%%

A_h = 1.1; %choose 1.1 as the initial value of A_h, adjust in the model
A_l = (1 - pi_lr(1,1)*A_h)/pi_lr(2,2);

%%%% Requirment 2, sd of Y_t in the model to be similar in the data%%%%
%%%% If sd(Y) is higher than it in data, lower the A_h

%%%% VFI Problem to find policy function first %%%%
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

%%%% The end of VFI problem %%%%

%%%% Simulation Problem %%%%

% generate the sequence of capital 
K_sim = zeros(1, 1000+1);
K_sim(1) = 5; % note ks are in (0,45]

for i = 2:1000+1
    [dif, lcn] = min(abs(K_sim(i-1)-k)); 
    K_sim(i) = g_h(lcn);
end

% generate the sequence of A
rng(9876);
A_sim = zeros(1,1000+1);
prob = rand(1,1000);
A_sim(1)= A_h;

for j = 2:1000+1
   if A_sim(j-1) == A_h
       if prob(j-1) < pi_hh
           A_sim(j) = A_l;
       else
           A_sim(j) = A_h;
       end
   else
       if prob(j-1) < pi_ll
           A_sim(j) = A_l;
       else
           A_sim(j) = A_h;
       end
   end
end


