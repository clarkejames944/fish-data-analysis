data {
    int<lower=0> N_EBSm;
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0> error;
    
    real<lower=0> sigma_int;
    real intercept;
    real beta;
}

transformed parameters {
    vector[N_EBSm] yhat;
    
    for (i in 1:N_EBSm){
         yhat[i] = intercept + beta * prev[i];
        }
}

model {
    sigma_int ~ cauchy(0,10);
    intercept ~ normal(0, sigma_int);
    beta ~ normal(0, 10);
    
    error ~ cauchy(0, 10);

    oto_size ~ normal(yhat, error);
}

generated quantities {
    vector[N_EBSm] sim_oto_size;

    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], error);
    }
}