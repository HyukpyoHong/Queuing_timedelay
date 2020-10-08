function [results, acceptance] = MCMC_function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, thin_rate, effnum)
%prod(normpdf(data, mean_list, var_list));
% timespan must be positive. i.e., timespan(1) > 0;
total_num = burn + thin_rate * effnum;

MCMC_sample = zeros(total_num, length(theta_init));
update_matrix = zeros(total_num, length(theta_init));

MCMC_sample(1,:) = theta_init;
update_matrix(2:end,1) = 1;

for rep = 2:total_num
    
    % Direct sampling for theta(1);
    mean_trj1 = mean_trajectory(timespan, [1 MCMC_sample(rep-1, 2) MCMC_sample(rep-1, 3) MCMC_sample(rep-1, 4)]);
    
    post_var = 1 / (sum(mean_trj1.^2 ./ var_list) + prior_var(1)^(-1));
    post_mean = post_var * (sum(data .* mean_trj1 ./ var_list) + prior_mean(1)/prior_var(1));
    
    MCMC_sample(rep, 1) = normrnd(post_mean, sqrt(post_var));
    
    % MCMC for theta(2): decay rate of X(t);
    % Random walk Metropolis algorithm (use normal proposal, symmetric)
    prop_mean2 = MCMC_sample(rep-1,2);
    prop_var2 = 0.000001; % tune this parameter
    
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
    
%     MCMC for theta(3): the shape parameter of delay distribution

%     MCMC_sample(rep,3) = MCMC_sample(rep-1,3);

    prop_mean3 = MCMC_sample(rep-1,3);
    prop_var3 = 0.1; % tune this parameter
    
    theta3_cand = normrnd(prop_mean3, sqrt(prop_var3));
    if theta3_cand > 0
        gamma_pri_alpha_3 = prior_mean(3)^2 / prior_var(3);
        gamma_pri_beta_3 = prior_var(3) / prior_mean(3);
        
        log_pri = log(gampdf(MCMC_sample(rep-1, 3), gamma_pri_alpha_3, gamma_pri_beta_3));
        log_pri_st = log(gampdf(theta3_cand, gamma_pri_alpha_3, gamma_pri_beta_3));
        
        mean_trj3 = mean_trajectory(timespan, [MCMC_sample(rep,1), MCMC_sample(rep, 2), MCMC_sample(rep-1, 3), MCMC_sample(rep-1, 4)]);
        mean_trj3_st = mean_trajectory(timespan, [MCMC_sample(rep, 1), MCMC_sample(rep, 2), theta3_cand, MCMC_sample(rep-1, 4)]);
        
        log_lik = sum(log(normpdf(data, mean_trj3, sqrt(var_list))));
        log_lik_st = sum(log(normpdf(data, mean_trj3_st, sqrt(var_list))));
        
        acceptance_ratio3 = exp(log_pri_st - log_pri + log_lik_st - log_lik);
        
        if rand(1) < acceptance_ratio3
            MCMC_sample(rep, 3) = theta3_cand;
            update_matrix(rep,3) = 1;
        else
            MCMC_sample(rep, 3) = MCMC_sample(rep-1,3);
        end
    else
        MCMC_sample(rep, 3) = MCMC_sample(rep-1, 3);
    end
   
    
    % === MCMC for theta(4): the rate parameter of delay distribution ===

    %     MCMC_sample(rep,4) = MCMC_sample(rep-1,4);

    prop_mean4 = MCMC_sample(rep-1,4);
    prop_var4 = 0.003; % tune this parameter
    
    theta4_cand = normrnd(prop_mean4, sqrt(prop_var4));
    if theta3_cand > 0
        gamma_pri_alpha_4 = prior_mean(4)^2 / prior_var(4);
        gamma_pri_beta_4 = prior_var(4) / prior_mean(4);
        
        log_pri = log(gampdf(MCMC_sample(rep-1, 4), gamma_pri_alpha_4, gamma_pri_beta_4));
        log_pri_st = log(gampdf(theta4_cand, gamma_pri_alpha_4, gamma_pri_beta_4));
        
        mean_trj4 = mean_trajectory(timespan, [MCMC_sample(rep,1), MCMC_sample(rep, 2), MCMC_sample(rep, 3), MCMC_sample(rep-1, 4)]);
        mean_trj4_st = mean_trajectory(timespan, [MCMC_sample(rep, 1), MCMC_sample(rep, 2), MCMC_sample(rep, 3), theta4_cand]);
        
        log_lik = sum(log(normpdf(data, mean_trj4, sqrt(var_list))));
        log_lik_st = sum(log(normpdf(data, mean_trj4_st, sqrt(var_list))));
        
        acceptance_ratio4 = exp(log_pri_st - log_pri + log_lik_st - log_lik);
        
        if rand(1) < acceptance_ratio4
            MCMC_sample(rep, 4) = theta4_cand;
            update_matrix(rep,4) = 1;
        else
            MCMC_sample(rep, 4) = MCMC_sample(rep-1,4);
        end
    else
        MCMC_sample(rep, 4) = MCMC_sample(rep-1, 4);
    end
   
    if rem(rep,50) == 0
       disp(rep) 
    end
    
   results = MCMC_sample;
   acceptance = update_matrix;
end

end
