function f = sa_perit_rhs(t, state)
    
    sa = state(1);
    coa = state(2);
    vw = state(3);
    st = state(4);
    m = state(5);
    mf = state(6);
    c = state(7);
    n = state(8);
    nf = state(9);
    fa = state(10);
    f1 = state(11);
    f2 = state(12);
    f3 = state(13);
    f4 = state(14);
    f5 = state(15);
    f6 = state(16);
    f7 = state(17);
    f8 = state(18);
    f9 = state(19);
    f10 = state(20);
    fn = state(21);
    fr = state(22);
    cfn = state(23);
    cfr = state(24);
    saf = state(25);
    fogen = fa + f1 + 2*f2 + 3*f3 + 4*f4 + 5*f5 + 6*f6 + 7*f7 + 8*f8 + 9*f9 + 10*f10 + cfn + cfr;
    fibrin = fogen - fa;
    large = cfn + cfr;
    

    %% Parameter values
    
    % Immune parameter values
    r = 0.3; % intrinsic bacterial growth rate
    ksa = 10^9.5; % SA carrying capacity
    gmf = 8e-6; % SA clearance by bound macrophages
    gnf = 8e-6; % SA clearance by bound neutrophils
    gm = 3e-7; % SA clearance by monocytes
    gn = 3e-7; % SA clearance by neutrophils
    nu = 0.3; % macrophage recruitment rate
    mmax = 1.25e6; % max macrophage count
    dm = 0.01; % natural macrophage death rate
    dms = 0.000000599/24; % death rate of macrophages due to SA
    eta = 9000/(24*1e9); % cytokine production
    kn = 1/(10^7.15); % cytokine inhibition by neutrophils
    km = 1/(10^7.15); % cytokine inhibition by macrophages
    theta = 1.65/10; % neutrophil recruitment rate
    nmax = 1.65e7; % max neutrophil count
    dn = 0.01; % neutrophil death rate
    dns = 0.000000000599/24; % neutrophil death rate due to SA
    dc = 10/24; % cytokine clearance rate
    
    % fibrin(ogen) parameter values
    kst = 3600;
    kpi = 4e-18*(1000*3600)*(1/340000)*(1/1.66e-18)*10000;
    kpg = 1e-16*(1000*3600)*(1/340000)*(1/1.66e-18)*10000;
    kfi = 1e-21*(1000*3600)*(1/340000)*(1/1.66e-18)*10000;
    kfg = 2e-17*(1000*3600)*(1/340000)*(1/1.66e-18)*10000;
    an = 1e1/5;
    bn = 1e2;
    am = 1e1/5;
    bm = 1e2;
    kcoa = 3e-9;
    dcoa = 3e-6;
    kvw = 3e-4;
    dvw = 3e-6;
    scoa = 1e-3;
    svw = 5e3;
    th = 1e5;
    mu = 5e1;
    rf = 350;
    asa = 8e1;
    bsa = 8e0;
    gfi = 1;
    rfi = r/1.2; % decreased growth rate from fibrin bound sa
    
    
    f = zeros(25,1); % need to return column vector
    
    f(1) = r*(sa + saf)*(1 - (sa + saf)/ksa) - gm*m*sa - gn*n*sa - gmf*mf*sa - gnf*nf*sa - asa*sa*large + bsa*saf;
    f(2) = kcoa*sa - dcoa*coa;
    f(3) = kvw*sa - dvw*vw;
    f(4) = scoa*coa*th + svw*vw*th - mu*st;
    f(5) = nu*c*(1 - (m+mf)/mmax) - am*m*fibrin + bm*mf - dms*sa*m - dm*m;
    f(6) = am*fibrin*m - bm*mf - dms*mf*sa - dm*mf;
    f(7) = eta*sa*((m+mf)/(1 + kn*(n+nf) + km*(m+mf))) - dc*c;
    f(8) = theta*c*(1 - (n+nf)/nmax) - dns*sa*n - dn*n - an*fibrin*n + bn*nf;
    f(9) = an*fibrin*n - bn*nf - dn*nf - dns*sa*nf;
    f(10) = rf*(1 - fogen/20) - kst*st*fa;
    f(11) = kst*st*fa - kpi*f1*(2*f1 + f2 + f3 + f4 + f5 + f6 + f7 + f8 + f9 + f10) - kpg*f1*fn;
    f(12) = kpi*f1*(f1 - f2);
    f(13) = kpi*f1*(f2 - f3);
    f(14) = kpi*f1*(f3 - f4);
    f(15) = kpi*f1*(f4 - f5);
    f(16) = kpi*f1*(f5 - f6);
    f(17) = kpi*f1*(f6 - f7);
    f(18) = kpi*f1*(f7 - f8);
    f(19) = kpi*f1*(f8 - f9);
    f(20) = kpi*f1*(f9 - f10);
    f(21) = kpi*f1*f10 - 2*kfi*fn*fn - kfg*fr*fn;
    f(22) = kfi*fn*fn;
    f(23) = 11*kpi*f1*f10 + kpg*fn*f1 - 2*kfi*fn*cfn - kfg*fr*cfn;
    f(24) = 2*kfi*fn*cfn + kfg*fr*cfn;
    f(25) = (rfi)*saf*(1 - (sa + saf)/ksa) + asa*sa*large - bsa*saf + gfi*(-gm*m*saf - gn*n*saf - gmf*mf*saf - gnf*nf*saf);
    
    
end