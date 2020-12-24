function logval = log_of_normpdf(x,mu,sigma)
    logval = -1/2 * log(2*pi) - log(sigma) - (x-mu).^2 ./ (2*sigma.^2);
end

