data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0, upper=4> eta; 
    
    real<lower=0, upper=100> epsilon;    
    
    real<lower=0, upper=100> sigma_alpha;
    real mu_alpha;
    real delta;
    
    vector[Ngroups] alpha;
    real beta1;
    real beta2;
}

transformed parameters {
    vector[N_EBSm] tau;
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        tau[i] = 1/(1 + exp(-100 * (prev[i] - eta)));
    }
    
    for (i in 1:N_EBSm){
         yhat[i] = (alpha[fishID[i]] + beta1 * prev[i]) * (1-tau[i]) + ((alpha[fishID[i]] + delta) + beta2 * prev[i]) * tau[i];
        }
}

model {///make sure you have a distribution for each parameter defined 
    eta ~ normal(1, 0.2);
    
    delta ~ normal(0, 0.5);
    sigma_alpha ~ cauchy(0, 0.5);
    mu_alpha ~ normal(0, 0.5);
    alpha ~ normal(mu_alpha, sigma_alpha);
    
    
    beta1 ~ normal(1, 0.5);
    beta2 ~ normal(1, 0.5);
    
    epsilon ~ cauchy(0, 1);

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    vector[N_EBSm] sim_oto_size;
    vector[N_EBSm] log_lik;


    for (i in 1:N_EBSm){
    sim_oto_size[i] = normal_rng(yhat[i], epsilon);
    }

    for (i in 1:N_EBSm){
        log_lik[i] = normal_lpdf(oto_size[i] | yhat[i], epsilon); 
    }
}