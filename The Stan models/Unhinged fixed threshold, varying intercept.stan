data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0.25, upper=4> eta; 
    
    real<lower=0, upper=100> epsilon;    
    
    real<lower=0, upper=100> sigma_alpha1;
    real<lower=0, upper=150> sigma_alpha2;
    vector[Ngroups] alpha1;
    vector[Ngroups] alpha2; 
    real beta;
}

transformed parameters {
    vector[N_EBSm] yhat;
    
    for (i in 1:N_EBSm) {
        if (prev[i] < eta) {
            yhat[i] = alpha1[fishID[i]] + beta * prev[i];
        } else {
            yhat[i] = alpha2[fishID[i]] + beta * prev[i];
        }
    }
}

model {///make sure you have a distribution for each parameter defined 
    eta ~ normal(0, 10);
    
    sigma_alpha1 ~ cauchy(0, 5);
    sigma_alpha2 ~ cauchy(0, 5);
    alpha1 ~ normal(0, sigma_alpha1);
    alpha2 ~ normal(0, sigma_alpha2);
    
    beta ~ normal(0, 10);
    
    epsilon ~ cauchy(0, 10);

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    vector[N_EBSm] sim_oto_size;

    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], epsilon);
    }
}