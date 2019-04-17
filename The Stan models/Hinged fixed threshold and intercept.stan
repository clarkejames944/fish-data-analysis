data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0> bp; //Do I want this to be a vector? or a real number across the groups (like intercept)
    
    real<lower=0> mu_bp;
    real<lower=0> sigma_bp;

    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    real<lower=0> sigma_int;
    real intercept; // one intercept per group
    real beta1;
    real beta2;
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
    yhat[i] = intercept + beta1 * prev[i] + beta2 * (prev[i] - bp) * x2[i];
    }
}

model {///make sure you have a distribution for each parameter defined
    mu_bp ~ normal(0, 5);
    sigma_bp ~ cauchy(0, 10);
    bp ~ normal(mu_bp, sigma_bp);

    sigma_int ~ cauchy(0, 10);
    intercept ~ normal(0, sigma_int);

    beta1 ~ normal(0, 5);
    beta2 ~ normal(0, 5);

    error ~ cauchy(0, 10);

    oto_size ~ normal(yhat, error);
}