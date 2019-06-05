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
    
    real<lower=0, upper=100> sigma_alpha;
    
    vector[Ngroups] alpha1;
    vector[Ngroups] alpha2; 
    real beta;
}

transformed parameters {
    vector[N_EBSm] tau;
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        tau[i] = 1/(1 + exp(-100 * (prev[i] - eta)));
    }
    
    for (i in 1:N_EBSm){
         yhat[i] = (alpha1[fishID[i]] + beta * prev[i]) * (1-tau[i]) + (alpha2[fishID[i]] + beta*prev[i]) * tau[i];
        }
}

model {///make sure you have a distribution for each parameter defined 
    eta ~ normal(1, 0.025);
    
    sigma_alpha ~ cauchy(0.05, 0.1);
    alpha1 ~ normal(0.275, sigma_alpha);
    alpha2 ~ normal(0.29, sigma_alpha);
    
    beta ~ normal(1, 0.1);
    
    epsilon ~ cauchy(0, 0.1);

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