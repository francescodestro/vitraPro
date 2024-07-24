function [t,x,x_bind,x_nucl,sum_h] = main_RK45(t_vect, x0, ...
    x0_bind, x0_nucl, D, r_bleed, Cin, Sin, p, sum_h,scaling)

    t0=t_vect(1);
    t_final=t_vect(2);
    t=t0;

    x=x0;
    x_bind=x0_bind;
    x_nucl=x0_nucl;

    %% RK45 parameters
    a2 = 1/5;
    a3 = 3/10;
    a4 = 4/5;
    a5 = 8/9;
    
    b11 = 1/5; 
    b21 = 3/40; 
    b31 = 44/45;
    b41 = 19372/6561;
    b51 = 9017/3168;
    b61 = 35/384;
    b22 = 9/40;
    b32 = -56/15;
    b42 = -25360/2187;
    b52 = -355/33;
    b33 = 32/9;
    b43 = 64448/6561;
    b53 = 46732/5247;
    b63 = 500/1113;
    b44 = -212/729;
    b54 = 49/176;
    b64 = 125/192;
    b55 = -5103/18656;
    b65 = -2187/6784;
    b66 = 11/84;
    
    e1 = 71/57600;
    e3 = -71/16695;
    e4 = 71/1920;
    e5 = -17253/339200;
    e6 = 22/525;
    e7 = -1/40;

    %% Get parameters from parameters object
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
    
    ind=0;
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

    k_lys=p.k_lys;

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
        ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);
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
        
        % step 2
        x2 = x + h .* (b11.*k1 );
        x_bind2 = x_bind + h .* (b11.*k1_bind );
        x_nucl2 = x_nucl + h .* (b11.*k1_nucl );
        t2 = t + h .* a2;

        [k2, k2_bind, k2_nucl]=inf_model(t2,x2,x_bind2,x_nucl2,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);

        % step3
        x3 = x + h .* (b21.*k1 + b22.*k2 );
        x_bind3 = x_bind + h .* (b21.*k1_bind + b22.*k2_bind);
        x_nucl3 = x_nucl + h .* (b21.*k1_nucl + b22.*k2_nucl);
        t3 = t + h .* a3;
        
        [k3, k3_bind, k3_nucl]=inf_model(t3,x3,x_bind3,x_nucl3,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);

        % step4
        x4 = x + h .* (b31.*k1 + b32.*k2 + b33.*k3);
        x_bind4 = x_bind + h .* (b31.*k1_bind + b32.*k2_bind + b33.*k3_bind);
        x_nucl4 = x_nucl + h .* (b31.*k1_nucl + b32.*k2_nucl + b33.*k3_nucl);
        t4 = t + h .* a4;
        
        [k4, k4_bind, k4_nucl]=inf_model(t4,x4,x_bind4,x_nucl4,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);

        % step5
        x5 = x + h .* (b41.*k1 + b42.*k2 + b43.*k3 + b44.*k4 );
        x_bind5 = x_bind + h .* (b41.*k1_bind + b42.*k2_bind + ...
            b43.*k3_bind + b44.*k4_bind );
        x_nucl5 = x_nucl + h .* (b41.*k1_nucl + b42.*k2_nucl + ...
            b43.*k3_nucl + b44.*k4_nucl );
        t5 = t + h .* a5;

        [k5, k5_bind, k5_nucl]=inf_model(t5,x5,x_bind5,x_nucl5,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);

        % step6
        x6 = x + h .* (b51.*k1 + b52.*k2 + b53.*k3 + b54.*k4 + b55.*k5 );
        x_bind6 = x_bind + h .* (b51.*k1_bind + b52.*k2_bind + b53.*k3_bind + ...
            b54.*k4_bind + b55.*k5_bind );
        x_nucl6 = x_nucl + h .* (b51.*k1_nucl + b52.*k2_nucl + b53.*k3_nucl + ...
            b54.*k4_nucl + b55.*k5_nucl );
        t6 = t + h;

        [k6, k6_bind, k6_nucl]=inf_model(t6,x6,x_bind6,x_nucl6,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);

        % Final step
        t_new=t+h;

        x_new=real(x + h.* (b61.*k1 + b63.*k3 + b64.*k4 + b65.*k5 + b66.*k6));
        x_bind_new=real(x_bind + h.* (b61.*k1_bind + b63.*k3_bind + ...
            b64.*k4_bind + b65.*k5_bind + b66.*k6_bind));
        x_nucl_new=real(x_nucl + h.* (b61.*k1_nucl + b63.*k3_nucl + ...
            b64.*k4_nucl + b65.*k5_nucl + b66.*k6_nucl));

        [k7, k7_bind, k7_nucl]=inf_model(t_new,x_new,x_bind_new,x_nucl_new,...
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
            ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);
    
        % estimate error
        e = h* (k1*e1 + k3*e3 + k4*e4 + k5*e5 + k6*e6 + k7*e7);
        e_bind = h* (k1_bind*e1 + k3_bind*e3 + k4_bind*e4 + k5_bind*e5 + ...
            k6_bind*e6 + k7_bind*e7);
        e_nucl = h* (k1_nucl*e1 + k3_nucl*e3 + k4_nucl*e4 + k5_nucl*e5 + ...
            k6_nucl*e6 + k7_nucl*e7);
    
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
    
            k1=k7;
            k1_bind=k7_bind;
            k1_nucl=k7_nucl;
    
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
                    ind_scB1Co,ind_scB2Co,n_bins,age_co,ind,startB1Co,startB2Co,scaling,k_lys);
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