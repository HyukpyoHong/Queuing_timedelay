%
T_max = 400;
Iter_max = 1000000;

theta = [35.4 0.015 5.89 1/0.89];

theta_init = 1/2 * [35.4 0.015 5.89 1/0.89];

[Xt, tspan, Xbirth, Xdeath] = Gillespie_delayX(theta, T_max, Iter_max);
data = cumsum(Xbirth - Xdeath);
timespan = (1:T_max)';
var_list = data + 1;
prior_mean = 1/2 * theta;
prior_var = [1000, 1, 10, 1];

burn = 0;
thin_rate = 1;
effnum = 20;

[rr, aa] = MCMC_function(data, timespan, var_list, prior_mean, prior_var, theta_init, burn, thin_rate, effnum);

cc = clock;
timestamp = [num2str(cc(1)) num2str(cc(2),'%02d') num2str(cc(3),'%02d') num2str(cc(4),'%02d') num2str(cc(5),'%02d') num2str(floor(cc(6)),'%02d')];
save(['Queueing_MCMC_' timestamp]);

