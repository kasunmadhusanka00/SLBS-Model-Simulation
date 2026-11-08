function dydt = slbs_model(t, y, beta, alpha, gamma1, gamma2, gamma3, delta1, delta2, theta, mu1, mu2)
    % Extract variables
    S = y(1);  % Susceptible
    L = y(2);  % Latent
    B = y(3);  % Breaking-out
    
    % Infection probability function
    f_infection = 1 - (1 - beta)^(L + B);
    
    % SLBS ODE system 
    dSdt = mu1 + gamma1 * B + gamma2 * L - f_infection * S - (delta1 + theta) * S;
    dLdt = mu2 + f_infection * S + theta * S + gamma3 * B - (alpha + gamma2 + delta1) * L;
    dBdt = alpha * L - (gamma1 + gamma3 + delta1 + delta2) * B;
    
    dydt = [dSdt; dLdt; dBdt];
end