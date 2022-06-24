function [ y ] = alpha_logistic( err_sqr,beta )

y=1/(1+exp(-beta*err_sqr));


end

