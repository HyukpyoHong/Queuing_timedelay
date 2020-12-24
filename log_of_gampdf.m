function logval = log_of_gampdf(x,shape,scale)
    logval = -log(gamma(shape-ceil(shape)+1)) - shape .* log(scale) + (shape-1) .* log(x) - 1./scale .* x;
    for d = 1:length(logval)
        for jj = 1:(ceil(shape(d)) - 1)
            logval(d) = logval(d) - log(shape(d) - jj);
        end 
    end
end

