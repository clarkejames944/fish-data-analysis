data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    vector<lower=0>[Ngroups] eta; 
    
    real<lower=0, upper=4> mu_eta;
    real<lower=0, upper=100> sigma_eta;

    real<lower=0, upper=100> epsilon;    
    
    real alpha; 
    real beta1;
    real beta2;
}

transformed parameters {
    vector[N_EBSm] tau; //indicator variable
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        if (prev[i] < eta[fishID[i]]) {
            tau[i] = 0;
        } else {
            tau[i] = 1;
        }
    }
    for (i in 1:N_EBSm){
    yhat[i] = alpha + beta1 * prev[i] + beta2 * (prev[i] - eta[fishID[i]]) * tau[i];
    }
}

model {///make sure you have a distribution for each parameter defined
    mu_eta ~ normal(0, 10);
    sigma_eta ~ cauchy(0, 5);
    eta ~ normal(mu_eta, sigma_eta);

    alpha ~ normal(0, 10);

    beta1 ~ normal(0, 10);
    beta2 ~ normal(0, 10);

    epsilon ~ cauchy(0, 10);

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    vector[Ngroups] intercept_after;
    real slope_after;
    
        for (i in 1:Ngroups){
            intercept_after[i] = alpha - eta[i]*beta2;
        }
        

        slope_after = beta1 + beta2;
}

generated quantities {
    vector[N_EBSm] sim_oto_size;

    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], epsilon);
    }
}