function dxdt = inf_model_lumped(t,x,mu,k_bind,k_deathT,k_deathI,k_rel,...
    K_s, Ys_T, Ys_I, k_degrV, Cin, r_bleed, Sin, D, scaling)
    
    T=x(1);
    V=x(2);
    I=x(3);
    S=x(4);

    b=r_bleed*D;

    growth_lim=max(S/(K_s+S),0);
    growthT=mu*T*growth_lim;
    bind1=k_bind*T*V;

    dxdt(1)=growthT-k_deathT*T-bind1+Cin*D/scaling-b*T;
    dxdt(2)=k_rel*I-k_bind*V*(T+I)-k_degrV*V-D*V;
    dxdt(3)=k_bind*T*V-k_deathI*I-b*I;
    dxdt(4)=D*(Sin-S)-(Ys_T*T+Ys_I *I)*S/(1e-2/scaling+S);

    dxdt=dxdt(:);
end