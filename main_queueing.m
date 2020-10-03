%%
clear; clc;

%%
maxT = 150;
data.time= (0:maxT)';   % x (mg / L COD) 
data.xdata = [0   0   1   2   4  11  16  22  30  36  42  46  57  64  74  78  84  89  93  94 104 109 118 122 125 128 138 138 145 148 153 155 150 153 153 157 158 155 160 171 183 188 191 190 191 189 192 191 194 202 191 187 191 193 190 190 187 178 183 179 170 171 177 177 179 184 186 180 178 182 183 184 185 181 192 187 185 181 185 190 190 187 190 181 182 186 185 182 187 188 196 196 190 193 199 200 205 203 209 209 204 205 206 215 216 224 230 226 231 233 232 235 241 237 237 228 225 215 221 219 219 214 214 214 218 214 212 214 214 213 217 213 213 211 221 220 222 216 216 222 217 223 215 212 203 204 203 207 207 204 199]'; % y (1 / h) 

%%
% Here is a plot of the data set.
figure(1); clf
plot(data.time,data.xdata,'s');
% xlim([0 maxT]); 
xlabel('time [min]'); ylabel('X(t) [number]');

%%

theta = [10, 0.05, 3.6, 0.6];


% theta(1): lambda_p, the birth rate of X
% theta(2): lambda_d, the death rate of X
% theta(3): alpha_X, the shape parameter of time delay gamma distribution of X
% theta(4): beta_X, the rate parameter of time delay gamma distribution of X

intfun1 = @(tau, theta) gampdf(tau, theta(3), 1/theta(4)) .* exp(theta(2)* tau);

modelfun = @(t,theta) theta(1)/theta(2) * (gamcdf(t, theta(3), 1/theta(4)) - exp(-theta(2) * t) .* integral(@(tau) intfun1(tau,theta), 0, t));

% modelfun = @(x,theta) theta(1)*x./(theta(2)+x);

ssfun    = @(theta,data) sum((data.xdata - modelfun(data.time,theta)).^2);

% plot(modelfun(data.xdata,[200 10]))
% hold on
% plot(data.xdata, data.ydata, 's');
%% 

fitted_value = zeros(size(data.time));
for tt = 1:length(data.time)
     fitted_value(tt) = modelfun(data.time(tt), theta);
end

figure(2);
plot(data.time,data.xdata,'s'); hold on;
plot(data.time, fitted_value); hold off;
% xlim([0 maxT]); 
xlabel('time [min]'); ylabel('X(t) [number]');

%% Initialization
% In this case the initial values for the parameters are easy to guess
% by looking at the plotted data. As we alredy have the sum-of-squares
% function, we might as well try to minimize it using |fminsearch|.
% [tmin,ssmin] = fminsearch(ssfun, [10;0.1;10;1], [] ,data);
% n = length(data.xdata);
% p = 2;
% mse = ssmin/(n-p) % estimate for the error variance


%%
% The Jacobian matrix of the model function is easy to calculate so we use
% it to produce estimate of the covariance of theta. This can be
% used as the initial proposal covariance for the MCMC samples by
% option |options.qcov| below.
J = [data.xdata./(tmin(2)+data.xdata), ...
     -tmin(1).*data.xdata./(tmin(2)+data.xdata).^2];
tcov = inv(J'*J)*mse


%%
% We have to define three structures for inputs of the |mcmcrun|
% function: parameter, model, and options.  Parameter structure has a
% special form and it is constructed as Matlab cell array with curly
% brackets. At least the structure has, for each parameter, the name
% of the parameter and the initial value of it. Third optional
% parameter given below is the minimal accepted value. With it we set
% a positivity constraits for both of the parameters.

params = {
    {'theta1', tmin(1), 0}
    {'theta2', tmin(2), 0}
    {'theta3', tmin(3), 0}
    {'theta4', tmin(4), 0}    
    };

%%
% The |model| structure holds information about the model. Minimally
% we need to set |ssfun| for the sum of squares function and the
% initial estimate of the error variance |sigma2|.

model.ssfun  = ssfun;
model.sigma2 = mse; % (initial) error variance from residuals of the lsq fit

%%
function [results, acceptance] = MCMC_function(mean_list, var_list, prior_mean, prior_var)
    theta = [10, 0.05, 3.6, 0.6];
    
    intfun1 = @(tau, theta) gampdf(tau, theta(3), 1/theta(4)) .* exp(theta(2)* tau);
    modelfun = @(t,theta) theta(1)/theta(2) * (gamcdf(t, theta(3), 1/theta(4)) - exp(-theta(2) * t) .* integral(@(tau) intfun1(tau,theta), 0, t));
    likelihood = prod(normpdf(data, mean_list, var_list));
    
    
end



x = 0:0.0001:0.9999;
sum((1-x).^15 + 15*x.*(1-x).^14 + 96*x.^2.*(1-x).^13 + 169*x.^3.*(1-x).^12 + 189*x.^4.*(1-x).^12+151*x.^5.*(1-x).^10 + 90*x.^6.*(1-x).^9 + 37*x.^7.*(1-x).^8+9*x.^8.*(1-x).^7 + x.^9.*(1-x).^6)
sum((1-x).^15 + 15*x.*(1-x).^14 + 87*x.^2.*(1-x).^13 + 17*x.^3.*(1-x).^12)

sum(15*x.*(1-x).^14 + 96*x.^2.*(1-x).^13 + 169*x.^3.*(1-x).^12 + 189*x.^4.*(1-x).^12+151*x.^5.*(1-x).^10 + 90*x.^6.*(1-x).^9 + 37*x.^7.*(1-x).^8+9*x.^8.*(1-x).^7 + x.^9.*(1-x).^6)
sum(15*x.*(1-x).^14 + 87*x.^2.*(1-x).^13 + 17*x.^3.*(1-x).^12)





