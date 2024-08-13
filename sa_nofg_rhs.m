function f = sa_nofg_rhs(t, state)
    
    sa = state(1);
    m = state(2);
    c = state(3);
    n = state(4);    

    %% Parameter values
    
    % Immune parameter values
    h = 3;
    gl = 0.12;
    r = 0.04; % intrinsic bacterial growth rate
    ksa = 10^9.5; % SA carrying capacity
    gm = 3e-7; % SA clearance by monocytes
    gn = 3e-7; % SA clearance by neutrophils
    nu = 0.3; % macrophage recruitment rate
    mmax = 1.25e6; % max macrophage count
    dm = 0.01; % natural macrophage death rate
    dms = 0.000000599/24; % death rate of macrophages due to SA
    eta = 9000/(24*1e9); % cytokine production
    kn = 1/(10^7.15); % cytokine inhibition by neutrophils
    km = 1/(10^7.15); % cytokine inhibition by macrophages
    theta = 1.66/10; % neutrophil recruitment rate
    nmax = 1.65e7; % max neutrophil count
    dn = 0.01; % neutrophil death rate
    dns = 0.000000000599/24; % neutrophil death rate due to SA
    dc = 10/24; % cytokine clearance rate
    
    
    f = zeros(4,1); % need to return column vector
    
    f(1) = r*sa*(1 - sa/ksa) - gl*(((m/sa)^h)/((0.01^h) + (m/sa)^h))*sa - gl*(((n/sa)^h)/((0.01^h) + (n/sa)^h))*sa;
    f(2) = nu*c*(1 - (m)/mmax) - dms*sa*m - dm*m;
    f(3) = eta*sa*((m)/(1 + kn*(n) + km*(m))) - dc*c;
    f(4) = theta*c*(1 - (n)/nmax) - dns*sa*n - dn*n;
    
end