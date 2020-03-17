clear all;
close all;
clc;
tic;
survivalprobs=xlsread('survival_probs_China.xlsx','A23:P42');
popgrowth=xlsread('survival_probs_China.xlsx','A47:P47');
timespan=linspace(2020,2100,16);
nrate=[timespan' popgrowth'];
save('survivalprobs','popgrowth','nrate');
%% Step1 : 
year1=1995;
movavperiods=4;			

case_productivity=0;	
case_growth=1;			
periodlength=5;			
%periodlength=1;
case_pen=1;				
						

case_UNscen=1;          
case_beta_calib=1;		

benchmark=1;			

save_results=1;         
maxit=50;               
phi=0.9;                
tol = 1e-10;            

if periodlength==1
	nage=75;
	Rage=46;            % first period of retirement 
	nr=nage-Rage+1;     % number of retirement years 
	nw=Rage-1;          % number of working years 	
elseif periodlength==5
	nage=15;
	Rage=10;            % first period of retirement 
	nr=nage-Rage+1;     % number of retirement years 
	nw=Rage-1;          % number of working years 	
else
	disp('wrong parameter period length');
	pause; 
end

%% Step2 :
nage1=75;
age=linspace(20,94,nage1);
efage=xlsread('efficiency_profile.xls',1,'A1:A45');
efage=efage/mean(efage);
year0=(year1-1950)/5+1;

if year0==round(year0)
	popgrowth=nrate(year0-movavperiods+1:year0,case_UNscen+1);
	popgrowth=mean(popgrowth);
	popgrowth=popgrowth/100;        % data is expressed in percentage numbers
	sp=survivalprobs(1:nage,year0-movavperiods+1:year0);
	sp=mean(sp,2);
else
	disp('year1 must be a multiple of 5');
end
if periodlength==5;		% transformation of annual data to 5-year data: survival probs and efficiency age
	if case_productivity==0
		efage1=zeros(nw,1);
		for i=1:1:nw
			efage1(i)=mean(efage((i-1)*5+1:i*5));		% average productivity
        end
		efage=efage1;
	elseif case_productivity==1
		efage=ones(nw,1);
    else
		disp('wrong parameter case_productivity');
		pause;
    end
end

ef=efage;
%% Step 3: 
alpha1=0.35;
delta=0.083;
rbbar=1.04;		% annual real interest rate on bonds
taun=0.28;
taunbar=taun;
taulbar=0.28;		% both taul+taup=0.28!!, see Mendoza, Razin (1994), p. 311
tauk=0.36;
taukbar=tauk;
tauc=0.05;
taucbar=tauc;
taup=0.124;
taupbar=0.124;
%replacement_ratio=0.50;		
% gross replacement ratio
replacement_ratio=0.352;
bybar=0.63;
%bybar=1.049;
%bybar=1;
if year1==2015
	gybar=0.18;
else
	gybar=0.239;
	gybar=0.18;
end

varphi=1;		% Frisch labor supply elasticity
varphi=0.3;
lbar=0.25;		% steady-state labor supply
eta1=2.0;		% 1/IES
kappa=3.619;
kappa=5.0;
kappa=10.0;
kappa=21.5;

if case_growth==1
	ygrowth=1.02;		% annual growth factor
else
	ygrowth=1.00;
end

if periodlength==5      % transformation of annual paramaters to 5-year values
	delta=1-(1-delta)^5;
	rbbar=rbbar^5;
	ygrowth=ygrowth^5;
	bybar=bybar/5;
	popgrowth=(1+popgrowth)^5-1;
end

if case_beta_calib==1
	beta1=ygrowth^(eta1)/rbbar;
	betaold=beta1;
elseif case_beta_calib==0
    load('output.mat','betacalib');
    beta1=betacalib;
else
	disp('wrong parameter case_beta_calib'); 
end
%% Step 4: Computation of initial value for capital-age profile
%			assuming constant labor n^s=0.3

nbar=0.2;
kbar=(alpha1/(rbbar+delta))^(1/(1-alpha1))*nbar;
kold=100;
nold=2;
wbar=wagerate(kbar,nbar);
dbar=interest(kbar,nbar);
ybar=production(kbar,nbar);
bigb=bybar*ybar;
kshare=kbar/(kbar+bigb);
% initialization of trbar
trbar=0.01*ybar;
% initial guess for pensions
pen=taup;		
%
% 	iteration over nr=1,..
%	to find an initial value for {k^1,..,k^nr}

nagegrid=linspace(nw+1,nw+nr,nr);
for iage=1:1:nr
	nage=nagegrid(iage);
	
	% computation of cohort mass
	mass=ones(nage,1);
	for i=2:1:nage
		mass(i)=mass(i-1)*sp(i-1)/(1+popgrowth);
    end
	mass=mass/sum(mass);

	% agents' policy function  
	aopt=zeros(nw+iage,1);      % optimal asset 
	copt=zeros(nw+iage,1);      % optimal consumption 
	nopt=lbar*ones(nw,1);       % optimal labor supply 

	% initial capital age-profile
	% assumption linear increase until retirement, 
	% thereafter linear decrease
	if iage==1
		kmax=kbar*5;
		kinitial = linspace(kmax/nw,(nw-1)*kmax/nw,nw-1);	% capital stock of the 2-year - nw-year old (capital stock of the one-year old is zero)
        kinitial = kinitial';       % column vector
        kinitial1=kmax/2;
		kinitial=[kinitial; kinitial1];
    else
		kinitial=[kinitial; kinitial(nage-2)/2];	
    end
    
	disp('initial values ss eoxgenous labor: ');
	y=ss_exog_labor(kinitial);
    
	kold=2*kbar;

	for q=1:1:maxit
		krit=abs((kbar-kold)/kbar);
		krit0=abs((nbar-nold)/nbar);
		wbar=wagerate(kbar,nbar); 
		dbar=interest(kbar,nbar);   
		kold=kbar;
		nold=nbar;

        xf = fsolve(@ss_exog_labor,kinitial);
        clc;

		if q>3 
            kinitial=xf; 
        end
		aopt=[0; xf; 0];
		k=kshare*aopt;
		b=(1-kshare)*aopt;
		
	
		for i=1:1:nw	
			copt(i)=(1-taun-taup)*wbar*ef(i)*nopt(i)+(1-tauk)*(dbar-delta)*k(i)+k(i)-ygrowth*k(i+1)+rbbar*b(i)-ygrowth*b(i+1)+trbar;
			copt(i)=copt(i)/(1+tauc);
        end

		if nr>0;
			for i=1:1:iage-1
				copt(i+nw) = pen+(1-tauk)*(dbar-delta)*k(nw+i)+k(nw+i)-ygrowth*k(nw+i+1)+rbbar*b(nw+i)-ygrowth*b(nw+i+1)+trbar;
				copt(i+nw) = copt(i+nw)/(1+tauc);
            end
			copt(iage+nw) = pen+(1-tauk)*(dbar-delta)*k(nw+iage)+k(nw+iage)+rbbar*b(nw+iage)+trbar;
			copt(iage+nw) = copt(iage+nw)/(1+tauc);
        end
	
% computation of the aggregate capital stock and employment nbar 
		anew=mass'*aopt(1:nage);	% aggregate wealth
		knew=kshare*anew;
		kbar=phi*kold+(1-phi)*knew;
		bigc=mass'*copt;            % total consumption
		taxes=taun*wbar*nbar+tauk*(dbar-delta)*kbar+tauc*bigc;
		temp=(1-sp(1:nage-1)).*mass(1:nage-1);
		bequests=temp'*aopt(2:nage);	
		bequests=bequests*(1+(1-tauk)*(dbar-delta)); 
		
		ybar=production(kbar,nbar);
		bigg=gybar*ybar;
		bigb=bybar*ybar;
		transfernew=taxes+bequests+bigb*((1+popgrowth)*ygrowth-rbbar)-bigg;
		ksharenew=kbar/(kbar+bigb);
% update of N, K, transfer, kshare, and pen
		trbar=phi*trbar+(1-phi)*transfernew;
		nopt1=ef.*nopt;
		nnew=nopt1'*mass(1:nw);
		nbar=phi*nold+(1-phi)*nnew;
		kshare=phi*kshare+(1-phi)*ksharenew;
		
		if case_beta_calib==1
			betanew=beta1*(1+(1-tauk)*(dbar-delta))/rbbar;
			beta1=phi*beta1+(1-phi)*betanew;
        end
			
		mean_labor=nopt1'*mass(1:nw);
		mean_labor=mean_labor/sum(mass(1:nw));

		if case_pen==1
			pennew=replacement_ratio*wbar*mean_labor;
			taupnew=pennew*sum(mass(nw+1:nage))/(wbar*nbar);
			taup=phi*taup+(1-phi)*taupnew;
		elseif case_pen==0
			pennew=taup*nbar*wbar/sum(mass(nw+1:nage));
        else			
			disp('wrong parameter for case_pen');
			pause;
        end	

		pen=phi*pen+(1-phi)*pennew;	
    end   % q */
    
end    % iage
    y=pen/(wbar*mean_labor);
    
% Step 5: Initial value for capital-age profile
%			assuming endogenous labor supply
ninitial=nopt;
xinitial = [kinitial; ninitial];
disp('initial values ss endogenous labor: ');
y=ss_endog_labor(xinitial);
y
%pause;

q=0;
while (q<maxit) || ((q>10) && abs((kbar-kold)/kbar)<tol )
    krit=abs((kbar-kold)/kbar);
    krit0=abs((nbar-nold)/nbar);
    q=q+1;
    wbar=wagerate(kbar,nbar); 
    dbar=interest(kbar,nbar);   
    kold=kbar;
    nold=nbar;

    xinitial = [kinitial; ninitial];
    xf = fsolve(@ss_endog_labor,xinitial);
    clc;
	kinitial=xf(1:nage-1);
	ninitial=xf(nage:nage+nw-1);

    aopt = [0; kinitial; 0];
	nopt = ninitial;
	
    % computation of copt	
	k=kshare*aopt;
	b=(1-kshare)*aopt;
	
	for i=1:1:nw	
		copt(i) = (1-taun-taup)*wbar*ef(i)*nopt(i)+(1-tauk)*(dbar-delta)*k(i)+k(i)-ygrowth*k(i+1)+rbbar*b(i)-ygrowth*b(i+1)+trbar;
		copt(i) = copt(i)/(1+tauc);
    end

	if nr>0
		for i=1:1:nr-1
			copt(i+nw) = pen+(1-tauk)*(dbar-delta)*k(nw+i)+k(nw+i)-ygrowth*k(nw+i+1)+rbbar*b(nw+i)-ygrowth*b(nw+i+1)+trbar;
			copt(i+nw) = copt(i+nw)/(1+tauc);
        end
		copt(nr+nw) = pen+(1-tauk)*(dbar-delta)*k(nw+nr)+k(nw+nr)+rbbar*b(nw+nr)+trbar;
		copt(nr+nw) = copt(nr+nw)/(1+tauc);
    end

    % computation of the aggregate capital stock and employment nbar 
    anew=mass'*aopt(1:nage);	% aggregate wealth
	knew=kshare*anew;
    kbar=phi*kold+(1-phi)*knew;
	bigc=mass'*copt;            % total consumption
	taxes=taun*wbar*nbar+tauk*(dbar-delta)*kbar+tauc*bigc;
	temp=(1-sp(1:nage-1)).*mass(1:nage-1);
	bequests=temp'*aopt(2:nage);
	bequests=bequests*(1+(1-tauk)*(dbar-delta)); 
	ybar=production(kbar,nbar);
	bigg=gybar*ybar;
	bigb=bybar*ybar;
	transfernew=taxes+bequests+bigb*((1+popgrowth)*ygrowth-rbbar)-bigg;
	ksharenew=kbar/(kbar+bigb);
    % update of N, K, transfer, kshare, and pen
	trbar=phi*trbar+(1-phi)*transfernew;
	nopt1=ef.*nopt;
	nnew=mass(1:nw)'*nopt1;
    nbar=phi*nold+(1-phi)*nnew;
	kshare=phi*kshare+(1-phi)*ksharenew;
	if case_beta_calib==1
		betanew=beta1*(1+(1-tauk)*(dbar-delta))/rbbar;
		beta1=phi*beta1+(1-phi)*betanew;
    end
	
	mean_labor=nopt1'*mass(1:nw);
	mean_labor=mean_labor/sum(mass(1:nw));
	if case_pen==1
		pennew=replacement_ratio*wbar*mean_labor;
		taupnew=pennew*sum(mass(nw+1:nage))/(wbar*nbar);
		taup=phi*taup+(1-phi)*taupnew;
	elseif case_pen==0
		pennew=taup*nbar*wbar/sum(mass(nw+1:nage));
    else			
		disp('wrong parameter for case_pen');
		pause;
    end	

	pen=phi*pen+(1-phi)*pennew;	
	
    
	mean_labor
	taup
    disp('nbar kbar transfers: ');
    [nbar, kbar, trbar]
	disp('kshare~q: ');
    [kshare, q]

end     % q 

disp('1. endogenous labor: ');
disp('=======================');
disp('aggregate consistency condition?');		
disp('y ');
ybar
disp('C+[(1+n)(1+g)-1+delta] K+g: ');
bigc+((1+popgrowth)*ygrowth-1+delta)*kbar+bigg;
disp('bequests: ');
bequests 
trbar
disp('transfers/GDP: ');
trbar/ybar
taup
disp('dependency ratio: ');
sum(mass(Rage:nage))/sum(mass(1:nw))
kshare
kbar
nbar
bigb
rbbar
disp('1+(1-tauk)*(d-delta): '); 
1+(1-tauk)*(dbar-delta)
disp('tax revenues: '); 
taun*wbar*nbar
beta1;
if periodlength==5;
	disp('annualized beta: ');
    beta1^(1/5)
end
mean_labor

xss=xf;


%
% the complete non-linear eqs system with aggregate variables
%
% aggregate non-linear eqs system: 
%
x0=[kinitial; ninitial; kbar; nbar; kshare; taup; trbar];
debtoutputratio=bybar;

y=ss_aggregate(x0);
y


%	Calibration
%
% the complete non-linear eqs system with aggregate variables
%
% 
disp('aggregate non-linear eqs system for calibration: ');
x0=[kinitial; ninitial; kbar; nbar; bigb; taup; trbar; beta1];
y=ss_aggregate_calib(x0);
y

% pause;
xagg = fsolve(@ss_aggregate_calib,x0);
y = ss_aggregate_calib(xagg);
	disp('mean(abs(y)), ss_aggregate_calib: ');
    mean(abs(y))
	beta1=xagg(nage+nw+5);
	kbar=xagg(nage+nw);
	nbar=xagg(nage+nw+1);
	bigb=xagg(nage+nw+2);
	trbar=xagg(nage+nw+4);
	taupbar=xagg(nage+nw+3);
	kshare=kbar/(kbar+bigb);
	bigg=gybar*production(kbar,nbar);
    ybar = production(kbar,nbar);
	disp('calibration of beta, debt, and G');
    beta1
    bigb
    bigg
 %   pause;
    
    
	y=ss_aggregate(xagg(1:nage+nw+4));
    disp('mean(abs(y)), ss_aggregate: ');
    y=mean(abs(y));
    y
	xstartss=xagg;
    
	debtb=bigb;
    x0=xagg(1:nage+nw+4);
	[cbench, savbench, penbench, taxesbench]=ss_computec(x0);
	
	lbench=xagg(nage:nage+nw-1);
	abench=[0; xagg(1:nage-1)];

% pause;

	trbench=trbar;
	ybench=ybar; 
	bbench=bigb;
	gbench=bigg;
    % save benchmark values for use in the transition dynamics
    save('Ch7_US_debt','trbench','ybench','cbench','penbench','taxesbench','xstartss','bbench','gbench');
    save('Ch7_US_debt','lbench','beta1','-append');
    
if case_beta_calib==1
     betacalib=beta1;
    save('Ch7_output.mat','betacalib');
end
    
    periods=linspace(20,90,15);
    periods1=linspace(20,60,9);
    figure
    subplot(2,2,1);
    plot(periods,abench);
    xlabel('Age');
    title('Wealth');
    subplot(2,2,2);
    plot(periods1,lbench);
    xlabel('Age');
    title('Working hours');
    subplot(2,2,3);
    plot(periods,cbench);
    xlabel('Age');
    title('Consumption');
    subplot(2,2,4);
    plot(periods,savbench*100);
    xlabel('Age');
    title('Savings rate');
    pause;
close all;

