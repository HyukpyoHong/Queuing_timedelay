function [results, acceptance] = MCMC_function(data, timespan, var_list, prior_mean, prior_var)
%prod(normpdf(data, mean_list, var_list));
% timespan must be positive. i.e., timespan(1) > 0;
theta_init = [5, 0.1, 3.6, 2];
burn = 0;
thin_rate = 1;
eff_num = 10;
total_num = burn + thin_rate * eff_num;

MCMC_sample = zeros(total_num, length(theta_init));
update_matrix = zeros(total_num, length(theta_init));

MCMC_sample(1,:) = theta_init;
update_matrix(2:end,1) = 1;

for rep = 2:total_num
    
    disp(111);
    % Direct sampling for theta(1);
    mean_trj1 = mean_trajectory(timespan, [1 MCMC_sample(rep-1, 2) MCMC_sample(rep-1, 3) MCMC_sample(rep-1, 4)]);
    
    post_var = 1 / (sum(mean_trj1.^2 ./ var_list) + prior_var(1)^(-1));
    post_mean = post_var * (sum(data .* mean_trj1 ./ var_list) + prior_mean(1)/prior_var(1));
    
    MCMC_sample(rep, 1) = normrnd(post_mean, post_var);
    
    % MCMC for theta(2): decay rate of X(t);
    % Random walk Metropolis algorithm (use normal proposal, symmetric)
    prop_mean2 = MCMC_sample(rep-1,2);
    prop_var2 = 0.01; % tune this parameter
    
    % a candidate for the next theta(2)
    
    theta2_cand = normrnd(prop_mean2, sqrt(prop_var2));
    if theta2_cand > 0
        gamma_pri_alpha_2 = prior_mean(2)^2 / prior_var(2);
        gamma_pri_beta_2 = prior_var(2) / prior_mean(2);
        
        log_pri = log(gampdf(MCMC_sample(rep-1, 2), gamma_pri_alpha_2, gamma_pri_beta_2));
        log_pri_st = log(gampdf(theta2_cand, gamma_pri_alpha_2, gamma_pri_beta_2));
        
        mean_trj2 = mean_trajectory(timespan, [MCMC_sample(rep,1), MCMC_sample(rep-1, 2), MCMC_sample(rep-1, 3), MCMC_sample(rep-1, 4)]);
        mean_trj2_st = mean_trajectory(timespan, [MCMC_sample(rep, 1), theta2_cand, MCMC_sample(rep-1, 3), MCMC_sample(rep-1, 4)]);
        
        log_lik = sum(log(normpdf(data, mean_trj2, sqrt(var_list))));
        log_lik_st = sum(log(normpdf(data, mean_trj2_st, sqrt(var_list))));
        
        acceptance_ratio2 = exp(log_pri_st - log_pri + log_lik_st - log_lik);
        
        if rand(1) < acceptance_ratio2
            MCMC_sample(rep, 2) = theta2_cand;
            update_matrix(rep,2) = 1;
        else
            MCMC_sample(rep, 2) = MCMC_sample(rep-1,2);
        end
    else
        MCMC_sample(rep, 2) = MCMC_sample(rep-1, 2);
    end
    
%     MCMC for theta(3) and theta(4): shape and rate parameters of delay
%     distribution;
    
    MCMC_sample(rep,3) = MCMC_sample(rep-1,3);
    MCMC_sample(rep,4) = MCMC_sample(rep-1,4);
    if rem(rep,10) == 0
       disp(rep) 
    end
end

end
