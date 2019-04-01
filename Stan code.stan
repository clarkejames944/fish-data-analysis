data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0, upper=2> bp[Ngroups];
    
    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    real intercept[Ngroups]; //intercept (etaint)
    vector[2] beta;
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

model {
    vector[N_EBSm] yhat;   //The conditional mean
    bp ~ normal(0.5, 0.5);
    intercept ~ normal(0, 0.5);
    beta ~ normal(0, 1);
    error ~ normal(0, 1);

    for (i in 1:N_EBSm) {
        yhat[i] = intercept[fishID[i]] + beta[1] * prev[i] + beta[2] * (prev[i] - bp[fishID[i]]) * x2[i];
         }
    
    oto_size ~ normal(yhat, error);
}