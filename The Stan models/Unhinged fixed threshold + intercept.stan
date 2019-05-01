data {
    int<lower=0> N_EBSm;
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0> bp; //Do I want this to be a vector? or a real number across the groups (like intercept)

    real<lower=0> mu_bp; //allow to find the mean value for the breakpoint
    real<lower=0, upper=100> sigma_bp; //allow to find the value of sd for the breakpoint
    
    real<lower=0, upper=100> error;    ///Leave a constant variance across fish individuals for now
    
    real<lower=0, upper=100> sigma_int1;
    real<lower=0, upper=150> sigma_int2;
    real intercept1; // two intercepts per fishID
    real intercept2;
    vector[2] beta;///two different beta values for each fish ID
}

transformed parameters {
    vector[N_EBSm] x2; //indicator variable
    vector[N_EBSm] yhat;
    
    for (i in 1:N_EBSm) {
        if (prev[i] < bp) {
            x2[i] = 0;
        } else {
            x2[i] = 1;
        }
    }
    
    for (i in 1:N_EBSm){
         yhat[i] = ((intercept1 + beta[1] * prev[i]) * (1-x2[i])) + ((intercept2 + beta[2]) * x2[i]);
        }
}

model {///make sure you have a distribution for each parameter defined 
    mu_bp ~ normal(0, 10);
    sigma_bp ~ cauchy(0, 5);
    bp ~ normal(mu_bp, sigma_bp);
    
    sigma_int1 ~ cauchy(0, 5);
    sigma_int2 ~ cauchy(0, 5);
    intercept1 ~ normal(0, sigma_int1);
    intercept2 ~ normal(0, sigma_int2);

    beta ~ normal(0, 10);
    
    error ~ cauchy(0, 10);

    oto_size ~ normal(yhat, error);
}
