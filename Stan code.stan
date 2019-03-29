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
    real slope_before[Ngroups]; //slope before breakpoint (etaslope1)
    real slope_after[Ngroups]; // part of slope after breakpoint (etaslope2)
}

transformed parameters {
    vector[N_EBSm] yhat;   //The conditional mean

    for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]) {
            yhat[i] = intercept[fishID[i]] + slope_before[fishID[i]] * (prev[i] - bp[fishID[i]]);
        } else {
            yhat[i] = intercept[fishID[i]] + slope_after[fishID[i]] * (prev[i] - bp[fishID[i]]);
        }
    }
}

model {
    bp ~ normal(0.5, 0.5);
    intercept ~ normal(0, 0.5);
    slope_before ~ normal(0.5, 1);
    slope_after ~ normal(0.5, 1);
    error ~ normal(0, 1);

    for (i in 1:N_EBSm) {
        oto_size[i] ~ normal(yhat[i], error);
         }
}