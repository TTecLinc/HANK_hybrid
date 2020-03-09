function y=inter_va(a,h,grid,a_grid,h_grid)
global sigma phi 
global beta g_A tr r pen delta_h Psi Xi w tau tre tw
global it ih ia is
global   s_grid e_grid
global V_agrid V_hgrid na nh ns
global max_a max_h min_a min_h ex_n

%Define the current grid
    %current_grid=V_agrid(it+1,:,:);
    current_grid=grid;
%Find the position of a h
    positionh=sum(h_grid<h);
    if positionh<1;
        positiona=sum(a_grid(:,1)<a);
    else
        positiona=sum(a_grid(:,positionh)<a);
    end
%Judge the position 
    if positiona<1;
        if positionh<1;
            y=current_grid(1,1,1);
        elseif positionh==nh;
            y=current_grid(1,1,nh);
        else
            q=(h_grid(positionh+1)-h)/(h_grid(positionh+1)-h_grid(positionh));
            y=q.*current_grid(1,1,positionh)+(1-q).*current_grid(1,1,positionh+1);
        end
    elseif positiona<na&&positiona>1;
        if positionh<1;  
            p=(a_grid(positiona+1,1)-a)/(a_grid(positiona+1,1)-a_grid(positiona,1));
            y=p.*current_grid(1,positiona,1)+(1-p).*current_grid(1,positiona+1,1);
        elseif positionh==nh;
            p=(a_grid(positiona+1,nh)-a)/(a_grid(positiona+1,nh)-a_grid(positiona,nh));
            y=p.*current_grid(1,positiona,nh)+(1-p).*current_grid(1,positiona+1,nh);
        else
            p=(a_grid(positiona+1,positionh)-a)/(a_grid(positiona+1,positionh)-a_grid(positiona,positionh));
            m=(a_grid(positiona+1,positionh+1)-a)/(a_grid(positiona+1,positionh+1)-a_grid(positiona,positionh+1));
            q=(h_grid(positionh+1)-h)/(h_grid(positionh+1)-h_grid(positionh));
            y=p.*q.*current_grid(1,positiona,positionh)+(1-p).*q.*current_grid(1,positiona+1,positionh)+...
                m.*(1-q).*current_grid(1,positiona,positionh+1)+(1-m).*(1-q).*current_grid(1,positiona+1,positionh+1);
        end
    else
        y=current_grid(1,na,nh);
    end   
end
