function mean_trj = mean_trajectory(timespan, theta)
length_of_T = length(timespan);
mean_trj = zeros(size(timespan));
intfun1 = @(tau, theta) gampdf(tau, theta(3), theta(4)) .* exp(theta(2)* tau);
for ii = 1:length_of_T
    t = timespan(ii);
    mean_trj(ii) = theta(1)/theta(2) * (gamcdf(t, theta(3), theta(4)) - exp(-theta(2) * t) .* integral(@(tau) intfun1(tau,theta), 0, t));
end
end
