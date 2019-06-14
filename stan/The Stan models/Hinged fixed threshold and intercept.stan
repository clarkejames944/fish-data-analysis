data {
    int<lower=0> N_EBSm;
    int<lower=0> Ngroups;
    int<lower=1, upper=Ngroups> fishID[N_EBSm];
    vector<lower=0>[N_EBSm] prev;
    vector<lower=0>[N_EBSm] oto_size;
}

parameters {
    real<lower=0, upper=4> eta;   //breakpoint

    real<lower=0, upper=10> epsilon;   //error
    
    real alpha;   //intercept
    real beta1;  //gradient
    real beta2;
}

transformed parameters {
    vector[N_EBSm] tau; //indicator variable
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        tau[i] = 1/(1 + exp(-100 * (prev[i] - eta)));
    }

    for (i in 1:N_EBSm){
        yhat[i] = alpha + beta1 * prev[i] + beta2 * (prev[i] - eta) * tau[i];
    }
}

model {///make sure you have a distribution for each parameter defined

    eta ~ normal(1, 0.2);

    alpha ~ normal(0, 0.2);

    beta1 ~ normal(1, 0.5);
    beta2 ~ normal(0, 0.5);

    epsilon ~ cauchy(0, 0.1);

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    real intercept_after;
    real slope_after;
    vector[N_EBSm] sim_oto_size;
    vector[N_EBSm] log_lik;

        

        intercept_after = alpha - (eta)*beta2;

        slope_after = beta1 + beta2;
          

        for (i in 1:N_EBSm){
            sim_oto_size[i] = normal_rng(yhat[i], epsilon);
        }

        for (i in 1:N_EBSm){
        log_lik[i] = normal_lpdf(oto_size[i] | yhat[i], epsilon); 
        }
}