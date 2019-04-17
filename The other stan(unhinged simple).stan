data {
    int<lower=1> N_EBSm;
    int<lower=1> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector<lower=0>[Ngroups] bp;
    
    real<lower=0> mu_bp;
    real<lower=0> sigma_bp;
    
    real<lower=0> error;    ///Leave a constant variance across fish individuals for now
    
    real<lower=0> sigma_int_bef;
    real<lower=0> sigma_int_after;
    real interceptbefore; //This should be different now as I have 2 different values within each group
    real interceptafter;
    real beta; //slope will be common for each group
}

transformed parameters {
    vector[N_EBSm] yhat;

        for (i in 1:N_EBSm) {
        if (prev[i] < bp[fishID[i]]){
            yhat[i] = interceptbefore+ beta * (prev[i] - bp[fishID[i]]); //got the same slopes =beta
         } else {
            yhat[i] = interceptafter + beta * (prev[i] - bp[fishID[i]]); //intercepts should be different
         }
        }
}
model {
    mu_bp ~ normal(0, 5);
    sigma_bp ~ cauchy(0, 10);
    bp ~ normal(mu_bp, sigma_bp);
    
    sigma_int_bef ~ cauchy(0, 10);
    sigma_int_after ~ cauchy(0, 10);
    interceptbefore ~ normal(0, sigma_int_bef); //intercept varies but is kept constant for the mean
    interceptafter ~ normal(0, sigma_int_after);
    beta ~ normal(0, 10); ///beta is constant 

    error ~ cauchy(0, 10);

    oto_size ~ normal(yhat, error);
}
