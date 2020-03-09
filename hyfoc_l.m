function y = hyfoc_l(l)
global sigma phi 
global beta g_A tr r pen delta_h Psi Xi w tau tre tw
global it ih ia is
global a_grid h_grid s_grid e_grid
global V_agrid V_hgrid na nh ns
global max_a max_h min_a min_h ex_n
    
    at1=s_grid(is).*(1+g_A);
    ht1=(h_grid(ih).*e_grid(it,is+ex_n,ih)).^(Psi).*Xi+h_grid(ih).*(1-delta_h);
            %Define A
            V=squeeze(V_agrid(it+1,:,:));
            
            A=beta.*h_grid(ih).*w.*(1-tau)/(1+g_A).*griddata(a_grid,h_grid,V,at1,ht1,'v4')/(1-phi);
            %Solve c-1(1-l-e) define B
            c=(1-e_grid(it,is+ex_n,ih)-l)/A;
            B=beta/(1+g_A).*griddata(a_grid,h_grid,V,at1,ht1,'v4')/A/phi;
            %Solve l
            y=-A+(c.^phi.*(1-l-e_grid(it,is+ex_n,ih)).^(1-phi)).^(-sigma).*c.^phi.*(1-l-e_grid(it,is+ex_n,ih)).^(-phi);
            %y=B-c.^(-1).*(1-l-e_grid(it,is+ex_n,ih));
end