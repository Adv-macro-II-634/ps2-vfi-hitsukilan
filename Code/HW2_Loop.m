%%%% Solve H2 with Loop instead %%%%
close all

%%%% set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;
A_h = 1.1; A_l = 0.678;
pi_hh = 0.977; pi_hl = 1-pi_hh;
pi_ll = 0.926; pi_lh = 1-pi_ll;
v_guess_h = 0;
v_guess_l = 0;
%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);
%%%% Set up consumption and return function
cons_h = zeros(num_k);
cons_l = zeros(num_k);
ret_h = zeros(num_k);
ret_l = zeros(num_k);

for i = 1: num_k
    for j = 1: num_k
        cons_h(i,j) = A_h * k(:,i)^ alpha + (1 - delta) * k(:,i) - k(:,j);
        cons_l(i,j) = A_l * k(:,i)^ alpha + (1 - delta) * k(:,i) - k(:,j);

        ret_h(i,j) = cons_h(i,j) ^ (1 - sigma) / (1 - sigma);
        ret_l(i,j) = cons_l(i,j)^ (1 - sigma) / (1 - sigma);

            if cons_h(i,j) < 0 
               ret_h(i,j) = -Inf;
            end
            
            if cons_l(i,j) < 0 
               ret_l(i,j) = -Inf;
            end
    end
end
%%%%Iteration
dis = 1; tol = 1e-06;
v_guess_h = zeros(1,num_k);
v_guess_l = zeros(1,num_k);
value_mat_h = zeros(num_k);
value_mat_l = zeros(num_k);

while dis>tol
    for i = 1:num_k
        for j = 1:num_k
            value_mat_h(i,j) = ret_h(i,j) + beta * (pi_hh * v_guess_h(:,j) + pi_hl * v_guess_l(:,j));
            value_mat_l(i,j) = ret_l(i,j) + beta * (pi_lh * v_guess_h(:,j) + pi_ll *  v_guess_l(:,j));

        end
    end
    
    [vfn_h, pol_indx_h] = max(value_mat_h, [], 2);
    vfn_h = vfn_h';
    
    [vfn_l, pol_indx_l] = max(value_mat_l, [], 2);
    vfn_l = vfn_l';
     dis = [max(abs(vfn_h - v_guess_h)); max(abs(vfn_l - v_guess_l))];
    
    v_guess_h = vfn_h;
    v_guess_l = vfn_l;
end

g_h = k(pol_indx_h); % policy function at A_h
g_l = k(pol_indx_l); % policy function at A_l



