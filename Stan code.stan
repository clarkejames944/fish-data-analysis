data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector[Ngroups] bp; //Do I want this to be a vector? or a real number across the groups (like intercept)

    real mu_bp; //allow to find the mean value for the breakpoint
    real<lower=0> sigma_bp; //allow to find the value of sd for the breakpoint
    
    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    vector[Ngroups] intercept; // one intercept per group
    vector[2] beta;///obviously then we have just two slope values for all the groups
}

transformed parameters {
    vector[N_EBSm] x2; //indicator variable
    
    for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]) {
            x2[i] = 0;
        } else {
            x2[i] = 1;
        }

}
}

model {///make sure you have a distribution for each parameter defined
    vector[N_EBSm] yhat;
    bp ~ normal(mu_bp, 100);
    mu_bp ~ normal (0, 10);

    intercept ~ normal(0, 1);
    beta ~ normal(0, 1);
    error ~ normal(0, 1);

    for (i in 1:N_EBSm){
         yhat[i] = intercept[fishID[i]] + beta[1] * prev[i] + beta[2] * (prev[i] - bp[fishID[i]]) * x2[i];
}

    oto_size ~ normal(yhat, error);
}