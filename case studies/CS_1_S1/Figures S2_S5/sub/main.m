function [t,x,x_bind,x_nucl,sum_h] = main(t_vect, x0, ...
    x0_bind, x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling)

    t0=t_vect(1);
    t_final=t_vect(2);
    t=t0;

    x=x0;
    x_bind=x0_bind;
    x_nucl=x0_nucl;

    % get parameters from parameters object
    startI1=p.startI1;
    endI1=p.endI1;
    startI2=0;
    endI2=0;
    startCo=0;
    endCo=0;
    startB1=p.startB1;
    startB2=0;
    endB1=p.endB1;
    endB2=0;
    startB1Co=0;
    startB2Co=0;
    ind_scCo=0;
    age_end_inf=p.age_end_inf;
    a=p.a;
    
    ind=p.ind;
    ind_scB1Co=0;
    ind_scB2Co=0;
    
    n_bins=p.n_bins;
    
    n_bins_inf=p.n_bins_inf;
    age_bins=p.age_bins;
    age_co=0;
    
    tau_death_I1=p.tau_death_I1;
    k_death_I1=p.k_death_I1;
    tau_death_I2=0;
    k_death_I2=0;
    tau_death_Co=0;
    k_death_Co=0;
    f_death_Co1=0;
    f_death_Co2=0;
    
    k_bind_I1=p.k_bind_I1;
    k_bind_I2=0;
    k_bind_Co1=0;
    k_bind_Co2=0;
    
    mu_T=p.mu_T;
    mu_I1=p.mu_I1;
    mu_I2=0;
    mu_Co=0;

    k_deathT=p.k_deathT;
    
    k_int_I1=p.k_int_I1;
    k_int_I2=0;
    k_int_Co1=0;
    k_int_Co2=0;
    
    eta_I1=p.eta_I1;
    eta_I2=0;
    eta_Co1=0;
    eta_Co2=0;
    
    f_repl_I1=p.f_repl_I1;
    f_repl_I2=0;
    f_repl_Co1=0;
    f_repl_Co2=0;
    
    k_repl_I1=p.k_repl_I1;
    k_repl_I2=0;
    k_repl_Co1=0;
    k_repl_Co2=0;

    f_rel_I1=p.f_rel_I1;
    f_rel_I2=0;
    f_rel_Co=0;
    k_rel_I1=p.k_rel_I1;
    k_rel_I2=0;
    k_rel_Co=0;
    
    k_degrBV_1=p.k_degrBV_1;
    k_degrBV_2=0;
    
    randV2_I1=0;
    randV1_I2=0;

    k_d_N1=p.k_d_N1;
    k_d_N2=0; 

    K_s=p.K_s;
    Ys_T=p.Ys_T;
    Ys_I1=p.Ys_I1; 
    Ys_I2=0; 
    Ys_C=0; 
    
    threshold=p.threshold;
    rtol=p.rtol;
    Dtau=p.Dtau;
    hmax0=p.hmax0;

    % first step selection
    [k1,k1_bind,k1_nucl]=inf_model(t0,x0,x0_bind,x0_nucl,...
        r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
        mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
        f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
        tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
        eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
        k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
        f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
        f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
        k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
        endCo, startB1,startB2, n_bins_inf,age_bins,...
        ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling);
    r=max(norm(k1./max(abs(x0),threshold), inf),...
        norm(k1_bind./max(abs(x0_bind),threshold), inf));
    h = 0.8*rtol^(1/3)/r;
    
    nofailed=true;

    %%% start
    while t<t_final
    
        % step control
        hmin = eps; %*abs(t);
        hmax=min(min(hmax0,Dtau-sum_h),t_final-t); % CFL condition
        h=max(hmin,min(h,hmax));
    
        [k2, k2_bind, k2_nucl]=inf_model(t+h/2,x+k1*h/2,x_bind+k1_bind*h/2,x_nucl+k1_nucl*h/2,...
            r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
            mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
            f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
            tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
            eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
            k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
            f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
            f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
            k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
            endCo, startB1,startB2, n_bins_inf,age_bins,...
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling);
         
        [k3, k3_bind, k3_nucl]=inf_model(t+3*h/4,x+3*k2*h/4,x_bind+3*k2_bind*h/4,x_nucl+3*k2_nucl*h/4,...
            r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
            mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
            f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
            tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
            eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
            k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
            f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
            f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
            k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
            endCo, startB1,startB2, n_bins_inf,age_bins,...
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling);
    
        t_new=t+h;
        ages_new=age_bins+sum_h+h;
        x_new=real(x+h*(2*k1+3*k2+4*k3)/9);
        x_bind_new=real(x_bind+h*(2*k1_bind+3*k2_bind+4*k3_bind)/9);
        x_nucl_new=real(x_nucl+h*(2*k1_nucl+3*k2_nucl+4*k3_nucl)/9);
        
        [k4, k4_bind, k4_nucl]=inf_model(t_new,x_new,x_bind_new,x_nucl_new,...
            r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
            mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
            f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
            tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
            eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
            k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
            f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
            f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
            k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
            endCo, startB1,startB2, n_bins_inf,age_bins,...
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling);
    
        % estimate error
        e = h*(-5*k1 + 6*k2 + 8*k3 - 9*k4)/72;
        e_bind=h*(-5*k1_bind + 6*k2_bind + 8*k3_bind - 9*k4_bind)/72;
        e_nucl=h*(-5*k1_nucl + 6*k2_nucl + 8*k3_nucl - 9*k4_nucl)/72;
    
        err_est=max(max(norm(e./max(max(abs(x),abs(x_new)), threshold), inf),...
            norm(e_bind./max(max(abs(x_bind),abs(x_bind_new)), threshold), inf)),...
            norm(e_nucl./max(max(abs(x_nucl),abs(x_nucl_new)), threshold), inf));
    
        if err_est <= rtol
            t=t_new;
            x=x_new;
            x_bind=x_bind_new;
            x_nucl=x_nucl_new;

            x(x<0)=0;
            x_bind(x_bind<0)=0;
            x_nucl(x_nucl<0)=0;
    
            k1=k4;
            k1_bind=k4_bind;
            k1_nucl=k4_nucl;
    
            sum_h=sum_h+h;
    
            temp = 1.25*(err_est/rtol)^(1/3);
            h=h*min(5,1/temp);    
            nofailed=true;
    
            %% age switch
            if sum_h >= Dtau || t>=t_final
                t=round(t,2);
    
                sum_h=0;
    
                % reallocate cells and virus
                % I1
                x(endI1)=x(endI1)+x(endI1-1);
                x(startI1+1:endI1-1)=x(startI1:endI1-2);
                x(startI1)=0;
    
                % reallocate virus bound to I1
                x_bind(endB1)=x_bind(endB1)+x_bind(endB1-1);
                x_bind(startB1+1:endB1-1)=x_bind(startB1:endB1-2);
                x_bind(startB1)=0;
    
                % reallocate nuclear virus in I1
                x_nucl(endB1)=x_nucl(endB1)+x_nucl(endB1-1);
                x_nucl(startB1+1:endB1-1)=x_nucl(startB1:endB1-2);
                x_nucl(startB1)=0;
    
                % update first slope after re-allocation
                [k1,k1_bind,k1_nucl]=inf_model(t,x,x_bind,x_nucl,...
                    r_bleed, Cin, Sin, D, k_bind_I1, k_bind_I2, k_bind_Co1, k_bind_Co2, mu_T, ...
                    mu_I1, mu_I2, mu_Co, k_deathT, k_death_I1, k_death_I2, k_death_Co, ...
                    f_death_Co1, f_death_Co2, tau_death_I1, tau_death_I2, ...
                    tau_death_Co, k_int_I1, k_int_I2, k_int_Co1, ...
                    eta_I1, eta_I2, eta_Co1, eta_Co2, k_int_Co2, k_repl_I1, ...
                    k_repl_I2, k_repl_Co1, k_repl_Co2, f_repl_I1, ...
                    f_repl_I2, f_repl_Co1, f_repl_Co2, f_rel_I1, f_rel_I2, ...
                    f_rel_Co, k_rel_I1, k_rel_I2, k_rel_Co, k_degrBV_1, k_degrBV_2, k_d_N1,...
                    k_d_N2, K_s,Ys_T,Ys_I1,Ys_I2,Ys_C, randV2_I1, randV1_I2, startI1,endI1,startI2,endI2,startCo, ...
                    endCo, startB1,startB2, n_bins_inf,age_bins,...
                    ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling);
            end
    
        else
            if nofailed
                nofailed = false;
                h = max(hmin, h * max(0.5, 0.8*(rtol/err_est)^(1/3)));
            else
                h = max(hmin, 0.5 * h);
            end
        end
        if h <= hmin
            warning('Step size %e too small at t = %e.\n',h,t);
            break
        end
    end
    
end