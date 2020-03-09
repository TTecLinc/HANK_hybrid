function y=hyinter_avalue(a,h)
global V_grid V_agird V_hgird na nh 
global theta rho r delta alpha gamma w phi beta 
global a_grid h_grid
global ia ih it
%Define the current grid
    current_grid=V_agird(it+1,:,:);
%Find the position of a h
    positionh=sum(h_grid<h);
    positiona=sum(a_grid(:,positionh)<a);
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