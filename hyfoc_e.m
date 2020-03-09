function y = hyfoc_e( e )
global sigma phi 
global beta g_A tr r pen delta_h Psi Xi w tau tre tw
global it ih ia is
global a_grid h_grid s_grid e_grid
global V_agrid V_hgrid na nh ns hgr
global max_a max_h min_a min_h ex_n 
    a1=(1+g_A).*s_grid(is);
    h1=(h_grid(ia,ih).*e).^(Psi).*Xi+h_grid(ia,ih).*(1-delta_h);
    %g_e=Xi.*Psi.*(e.*h_grid(ih)).^(Psi-1).*h_grid(ih);
    Va=squeeze(V_agrid(it+1,:,:));
    Vh=squeeze(V_hgrid(it+1,:,:));

    y=w.*(1-tau).*griddata(a_grid,h_grid,Va,a1,h1,'v4')/(1+g_A)-griddata(a_grid,h_grid,Vh,a1,h1,'v4').*Xi.*Psi.*(e.*h_grid(ia,ih)).^(Psi-1);
end