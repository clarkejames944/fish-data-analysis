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
    
    real alpha; 
    real beta1;
    real beta2;
    real beta3;
}

transformed parameters {
    vector[N_EBSm] tau; //indicator variable
    vector[N_EBSm] yhat;

    for (i in 1:N_EBSm) {
        if (prev[i] < eta) {
            tau[i] = 0;
        } else {
            tau[i] = 1;
        }
    }
    for (i in 1:N_EBSm){
    yhat[i] = alpha + beta1 * prev[i] + beta2 * (prev[i] - eta) * tau[i] + beta3*tau[i];
    }
}

model {///make sure you have a distribution for each parameter defined

    eta ~ normal(0, 10);

    alpha ~ normal(0, 10);

    beta1 ~ normal(0, 10);
    beta2 ~ normal(0, 10);
    beta3 ~ normal(0, 10);

    epsilon ~ cauchy(0, 10);

    oto_size ~ normal(yhat, epsilon);
}

generated quantities {
    real intercept_after;
    real slope_after;

        intercept_after = alpha - (eta)*beta2 + beta3;

        slope_after = beta1 + beta2;
}