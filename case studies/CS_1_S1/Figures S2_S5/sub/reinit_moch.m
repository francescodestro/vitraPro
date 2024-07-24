function p = reinit_moch(p,age_bins)
   
    %% Mesh and state vectors indexing
    % Mesh
    n_bins=length(age_bins);
    
    age_bins_inf=age_bins(age_bins<=p.age_end_inf);
    [~,n_bins_inf]=max(age_bins_inf);
    
    endI1=2+n_bins;

    endB1=n_bins;
    
    %% Parameters infected cells
    f_rel_I1=zeros(1,length(age_bins));
    f_rel_I1(age_bins>=p.tau_rel_on_I1)=1;
    f_rel_I1(age_bins>=p.tau_rel_off_I1)=0;
    f_repl_I1=(0.5+0.5*tanh((age_bins-p.replStart_I1)/0.3)).*(0.5+0.5*tanh((p.replEnd_I1-age_bins)/1)); % smooth transition
    n1=0.5+0.5*(tanh((age_bins(1:n_bins_inf)-p.tau0_bind_I1)/.01));
    k_bind_I1=p.k_bind0_I1*((1-n1)+n1.*exp(-p.beta_bind_I1*(age_bins(1:n_bins_inf)-p.tau0_bind_I1)));
   
    %% Prepare updated parameters object to pass to solver
    p.endI1=endI1;
    p.endB1=endB1;
    p.n_bins=n_bins;
    p.n_bins_inf=n_bins_inf;
    p.k_bind_I1=k_bind_I1;
    p.f_repl_I1=f_repl_I1;
    p.f_rel_I1=f_rel_I1;
