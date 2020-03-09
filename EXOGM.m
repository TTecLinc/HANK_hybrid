%Parameter
clear;

global V_grid V_agird V_hgird na nh 
global theta rho r delta alpha gamma w phi beta 
global a_grid h_grid
global ia ih it
global c_grid e_grid

theta=0.5;
rho=0.04;
r=0.05;
delta=0.05;
alpha=0.35;
gamma=1;
w=0.1;
phi=0.5;
beta=0.95;
nh=25;
na=25;
nT=6; 
%Step1: Discrete Space
h_grid=linspace(1,3,nh);
a_grid=linspace(0.01,10,na);
c_grid=zeros(nT,na,nh);
e_grid=zeros(nT,na,nh);
h1_grid=zeros(nT,na,nh);
a1_grid=zeros(nT,na,nh);
V_grid=zeros(nT,na,nh);
V_agird=zeros(nT,na,nh);
V_hgird=zeros(nT,na,nh);


c_profile=zeros(nT,1);
e_profile=zeros(nT,1);
a_profile=zeros(nT,1);
h_profile=zeros(nT,1);
%Step2: The last Period
for ia=1:na;
    for ih=1:nh;
        c_grid(nT,ia,ih)=a_grid(ia)+w.*h_grid(ih);
        if c_grid(nT,ia,ih)<=0;
            c_grid(nT,ia,ih)=1e-8;
        end
        V_grid(nT,ia,ih)=c_grid(nT,ia,ih).^(1-theta)/(1-theta);
        V_agird(nT,ia,ih)=c_grid(nT,ia,ih).^(-theta);
        V_hgird(nT,ia,ih)=w.*c_grid(nT,ia,ih).^(-theta);
    end
end 
%Step3: Iterate Backwards
for it=nT-1:-1:1;
    for ia=1:na;
        for ih=1:nh;               
            options = optimoptions(@fsolve,'Algorithm','levenberg-marquardt','MaxFunEvals',1500);
            solution=fsolve(@(x)opt_equation(x(1),x(2)),[0.5;0.6],options);
            c_grid(it,ia,ih)=solution(1);
            e_grid(it,ia,ih)=solution(2);
            at1=(1+r).*(a_grid(ia)+w.*h_grid(ih)-c_grid(it,ia,ih)-e_grid(it,ia,ih));
            ht1=(1-delta).*(h_grid(ih)+gamma/alpha.*e_grid(it,ia,ih).^alpha);
            h1_grid(it+1,ia,ih)=ht1;
            a1_grid(it+1,ia,ih)=at1;
            V_grid(it,ia,ih)=1/(1-theta).*c_grid(it,ia,ih).^(1-theta)+beta.*(1-phi/(1+ht1)).*inter_value(at1,ht1);
            V_agird(it,ia,ih)=c_grid(it,ia,ih).^(-theta);
            V_hgird(it,ia,ih)=(w+1/gamma.*e_grid(it,ia,ih).^(1-alpha)).*c_grid(it,ia,ih).^(-theta);            
        end
    end
end

for i=1:nT;
    c_profile(i)=mean(mean(c_grid(i,:,:)));
    e_profile(i)=mean(mean(e_grid(i,:,:)));
    a_profile(i)=mean(mean(a1_grid(i,:,:)));
    h_profile(i)=mean(mean(h1_grid(i,:,:)));
end

plot(c_profile);
title('Consumption Profile');

figure();
plot(e_profile);
title('Education Profile');

figure();
plot(a_profile);
title('Asset Profile');

figure();
plot(h_profile);
title('Human Capital Profile');


