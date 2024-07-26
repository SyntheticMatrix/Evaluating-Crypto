var 
 pi     ${\hat{\pi_m}}$      (long_name='price level in traditional money sector')
 z      ${z}$                (long_name='monetary demand shock')
 x      ${\hat{y}}$          (long_name='output gap')
 ej     ${\hat{\delta^m}}$   (long_name='traditional money demand')
 ejd    ${\hat{\delta^c}}$   (long_name='cryptocurrency demand') 
 CC     ${\hat{cc}}$         (long_name='real balance of cryptocurrency')
 M      ${\hat{m}}}$         (long_name='real balance of traditional money')
 rc     ${i}$                (long_name='official interest rate')
 hi     ${\hat{\chi}}$       (long_name='exchange rate')
 fi     ${\hat{\nu}_t}$       (long_name='cryptocurrency production cost')
 ksit   ${\hat{\xi}^e}$     (long_name='cryptocurrency production cost AR(1)')
 muu    ${e}$               (long_name='monetary policy shock') 
  x_obs                     (long_name='data of output gap')
 pi_obs                     (long_name='data of price level in traditional money')
 r_obs                      (long_name='data of federal fund rate') 
 cc_obs                     (long_name='data of real cryptocurrency coin accumulation') 
 m_obs                      (long_name='data of M2') 
 pji    ${\hat{\pi_c}}$      (long_name='price level in cryptocurrency setor')       
 pc_obs                      (long_name='data of cryptocurrency sectoral price')
 a       ${a}$                (long_name='TFP')
 xx      ${y}$                (long_name='output')
 y_nat    ${y^nat}$         (long_name='natural rate of output') 
  r_nat   ${r^nat}$         (long_name='natural rate of interest')
 c     ${\hat{c}}$       (long_name='consumption')
 n    ${\hat{n}}$          (long_name='labor supply')
 s    ${\hat{s}}$          (long_name='relative price')
 ymt_gap ${\hat{y}_m}$     (long_name='traditional money sectoral output gap')
 mcm      ${\hat{mc}_m}$    (long_name='traditional money sectoral marginal cost')
 mcc      ${\hat{mc}_c}$    (long_name='cryptocurrency sectoral marginal cost')
 ymt      ${y_m}$          (long_name='traditional money sectoral output')
 p_m_hat  ${\hat{p}_m}$    (long_name='unilateral relative price of traditional money sector')
 p_c_hat  ${\hat{p}_c}$    (long_name='unilateral relative price of cryptocurrency sector')
 r_d      ${r^d}$          (long_name='deposit rate') 
 w_hat    ${\hat{w}}$      (long_name='real wage')
pi_agg    ${\pi^{agg}}$    (long_name='Aggregated price level')
d         ${d}$            (long_name='deposit')
er        ${er}$           (long_name='excess return of holding cryptocurrency')
r_b       ${r^b}$          (long_name='government bond rate')
er1       ${er_1}$         (long_name='1st argument of excess return')
er2       ${er_2}$         (long_name='2nd argument of excess return')
CH                          (long_name='cryptocurrency sectoral output gap for welfare analysis')
div       ${div}$           (long_name='dividend of investor')
;

varexo epsx   ${\epsilon^{\delta^m}}$  (long_name='traditional money demand shock')
       epsz   ${\epsilon^{z}}$         (long_name='monetary demand shock')
       epsdn   ${\epsilon^{}}$          (long_name='cryptocurrency production shock - exchange rate depreciation')
       epsdj   ${\epsilon^{\delta^c}}$  (long_name='cryptocurrency demand shock') 
       epsa    ${\epsilon^a}$           (long_name='TFP shock')
       epsm    ${\epsilon^m}$           (long_name='monetary policy shock')
;


parameters ksi nu_d nu_ci zeta_cb_para alpha sigma phi niz betta theta1 theta2 eta  gamma9  ni rhoa  rhoz rhoe varphi rho ksiphi rhoksi rhonit omega2 omega3 gamma4 gamma8 gamma2 gamma3 gamma6 gamma7 rhor rhoy rhopi rhomu gamma1 gamma5 zeta rhom rhoed zetaa ;


zetaa=0.05;    % basic externality cost
alpha = 0.25; %share of labor input in production
sigma = 1;    %coefficient of risk aversion
phi = 1;     %inverse frisch elasticity of labor supply
betta = 0.99;  %discount factor
theta1= 0.75; %probability of not adjusting prices in dollar sector
theta2= 0.25; %probability of not adjusting prices in non-dollar sector
eta = 9;      %elsticity of substitution among consumption goods
ni = 0.4;     %size of non-dollar sector
rhoa =  0.9;  % techonology shock AR(1) parameter
rhoz = 0.5;  % monetary demand shock AR(1) parameter
rhoe = 0.5;  % traditional money demand shock AR(1) parameter
varphi  = 5;     % natural rate of output definition following Gali (2015)
ksiphi=0.5;     % cryptocurrency production cost deviation
rho=0.9;        % cryptocurrency production cost parameter
rhoksi=0.6;  % cryptocurrency production cost AR(1) parameter
rhonit=0.6;  % cryptocurrency production cost AR(1) parameter
omega2=0.2; % output elasticity to real balance of traditional money
omega3=0.05; % output elasticity to real balance of cryptocurrency
zeta=0.2;   % share of crypto-medium of exchange
gamma1=0.015;   % income elasticity of traditional money demand
gamma2= 0.15;  % interest semi-elasticity of traditional money demand
gamma3=0.9;  % elasticity of real balance of traditional money w.r.t traditional money demand shock
gamma4=0.5;  % cross elasticity of traditional money demand and cryptocurrency demand
gamma5=0.015;   % income elasticity of cryptocurrency demand
gamma6=0.15;  % Interest semi-elasticity of cryptocurrency demand
gamma7=0.8;  % elasticity of real balance of cryptocurrency w.r.t cryptocurrency demand shock
gamma8=0.8 ;  % cross elasticity of cryptocurrency demand and traditional money demand
gamma9=0.6;  % elasticity of real balance of cryptocurrency w.r.t cryptocurrency exchange rate periodic difference

rhor=0.8;    % interest rate smoothing
rhoy=0.2;     % Taylor rule coefficient on output
rhopi=1.8;    % Taylor rule coefficient on inflation
rhomu=0.2;    % Taylor rule coefficient on change in traditional money supply

rhom=0.5;    % monetary policy shock parameter
rhoed=0.7;   % cryptocurrency demand shock AR(1) parameter

niz=0.1;        % cryptocurrency externality cost adjustment
zeta_cb_para=0.1; % transaction cost
nu_ci=0.4;        % share or probability of investor on government bond
nu_d=0.7;       % share of deposit over M2 (based on U.S. June 2024 )
ksi=0.1;      % cryptocurrency demand Euler equation SoC related product



model(linear);

//(1)Cryptocurrency demand Euler equation
CC=gamma5*(x-zeta*(CC)-zeta*hi)+gamma7*ejd+gamma8*ej-gamma9*(hi-hi(+1))-gamma8*M-gamma6*rc+gamma9*(zetaa/niz-ni)*hi(+1);

//(2)Traditional money demand Euler equation
M=gamma1*(x-zeta*(CC)-zeta*hi)+gamma3*ej+gamma4*ejd-gamma4*(hi)-gamma4*CC-gamma2*rc;

//(3)IS curve w.r.t aggregated output and price level
x = x(+1)-(sigma)*(rc-pi_agg-r_nat)+omega2*(M-M(+1)-ej+ej(+1))+omega3*(hi-hi(+1)+CC-CC(+1)-ejd+ejd(+1))-(zeta)*(CC(+1)-CC-hi+hi(+1));

//(4)&(5)NKPC w.r.t traditional money sector
pi = betta*pi(+1)+mcm+((1-theta1)*(1-betta*theta1)/theta1)*(1-ni)*(s);
mcm=((1-theta1)*(1-theta1*betta)/theta1)*((1-alpha)/(1-alpha+alpha*eta))*((sigma))*(x-zeta*CC-zeta*hi-omega2*(M-ej)-omega3*(hi+CC-ejd))+((1-theta1)*(1-theta1*betta)/theta1)*((1-alpha)/(1-alpha+alpha*eta))*((phi+alpha)/(1-alpha))*ymt;

//(6)&(7)NKPC w.r.t cryptocurrency sector full anonymity
pji = betta*pji(+1)-((betta-ksi)/(1+betta-ksi))*((1-theta2)*(1-betta*theta2)/theta2)*ni*s+mcc;
mcc=((1-theta2)*(1-theta2*betta)/theta2)*((1-alpha)/(1-alpha+alpha*eta))*((sigma))*(x-zeta*CC-zeta*hi-(M-ej+CC+hi-ejd))*(ksi)+((1-theta2)*(1-theta2*betta)/theta2)*((1-alpha)/(1-alpha+alpha*eta))*((phi+alpha)/(1-alpha))*(zeta*CC+zeta*hi);

//(8)Aggregate inflation
pi_agg=(1-ni)*pji+(ni)*pi;

//(9)&(10)&(11)Relative price equation
s=s(-1)+pji-pi+hi-hi(-1);

p_c_hat=(1-ni)*s;

p_m_hat=-ni*s;

//(12)Cryptocurrency production FoC
hi=-rho*fi;

//(13)&(14)Cryptocurrency production cost AR(1) system
fi=ksiphi*ksit+(1-ksiphi)*er;

ksit=rhoksi*ksit(-1)+epsdn;

//(15)&(16)&(17)excess return
er1=nu_ci*(r_b-hi(-1));

er2=(1-nu_ci)*(r_d-hi(-1));

er=er1+er2;

//(18)Taylor rule
rc=rhor*rc(-1)+(1-rhor)*rhoy*x+(1-rhor)*rhopi*(pi)+muu+(1-rhor)*rhomu*(M-M(-1));

//(19)Monetary policy shock
muu=rhom*muu(-1)+epsm;

//(20)Traditional money demand shock
ej = rhoe*ej(-1)+epsx;

//(21)Cryptocurrency demand shock
ejd = rhoed*ejd(-1)+epsdj;

//(22)monetary demand shock
z = rhoz*z(-1) + epsz;

//(23)Production function
xx = a+(1-alpha)*n;

//(24)market clearing
xx = c+zeta*(hi+CC);

//(25)Real wage
w_hat=(1-alpha)*(a-alpha*n);


//(26)TFP shock
a=rhoa*a(-1)+epsa;

//(27)natural rate of interest
r_nat=-sigma*((1+varphi)/(sigma*(1-alpha)+varphi+alpha))*(1-rhoa)*a+(1-rhoz)*z;

//(28)natural rate of output
y_nat=((1+varphi)/(sigma*(1-alpha)+varphi+alpha))*a;

//(29)output gap
x=xx-y_nat;

//(30)&(31)sectoral market clearing

ymt=c-eta*p_m_hat;
ymt_gap=ymt-a;

//(32)&(33)Non-arbitrage
r_b=r_d+M*zeta_cb_para;
r_d=rc;

//(34)Monetary authority clearing:
d=M*nu_d;

//(35)Investor constraint
div=xx-w_hat-n;

//(36)cryptocurrency sectoral output gap
CH=zeta*(CC+hi)-a;


//estimation identification strategy
x=x_obs;
pi=pi_obs;
rc=r_obs;
M=m_obs;
CC=cc_obs;
pji=pc_obs;




end;

steady;


//Estimation activation

varobs x_obs, pi_obs, r_obs, cc_obs, m_obs pc_obs; 

estimated_params;

zeta, NORMAL_PDF, 0.2, 0.01;
zeta_cb_para, NORMAL_PDF, 0.1, 0.05;
omega2, NORMAL_PDF, 0.2, 0.1;
omega3, NORMAL_PDF, 0.05, 0.02;
theta2, GAMMA_PDF, 0.25, 0.1;
gamma1, GAMMA_PDF, 0.015, 0.005;
gamma2, GAMMA_PDF, 0.15, 0.01;
gamma3, GAMMA_PDF, 0.8, 0.1;
gamma4, GAMMA_PDF, 0.5, 0.01;
gamma5, GAMMA_PDF, 0.015, 0.005;
gamma6, GAMMA_PDF, 0.15, 0.01;
gamma7, GAMMA_PDF, 0.7, 0.1;
gamma8, GAMMA_PDF, 0.8, 0.1;
gamma9, GAMMA_PDF, 0.6, 0.1;
ksi, GAMMA_PDF, 0.1, 0.05;

rhoe,           BETA_PDF,0.5,0.01;
rhonit,           BETA_PDF,0.6,0.05;
rhoz,           BETA_PDF,0.5,0.05;
rhoed,           BETA_PDF,0.7,0.01;
rhom,          BETA_PDF,0.5,0.05;
rhoa,          BETA_PDF,0.9,0.05;

rhoy, BETA_PDF, 0.2, 0.005;
rhomu, NORMAL_PDF, 0.2, 0.05;
rhor, BETA_PDF, 0.8, 0.05;
rhopi, GAMMA_PDF, 1.8, 0.005;

stderr epsz,              INV_GAMMA2_PDF,0.5,INF;
stderr epsx,              INV_GAMMA2_PDF,0.5,INF;
stderr epsdn,              INV_GAMMA2_PDF,0.5,INF;
stderr epsdj,              INV_GAMMA2_PDF,0.5,INF;
stderr epsm,              INV_GAMMA2_PDF,0.5,INF;
stderr epsa,              INV_GAMMA2_PDF,0.5,INF;
end;

estimation(datafile=ccestrrnew,diffuse_filter, mode_check,mode_compute=6,mh_replic=20000,mh_nblocks=2, mh_jscale=0.8, mh_drop=0.2) x pi rc M CC pji; 
//stoch_simul(order=1, irf=40)x pi rc M CC pji s ymt  c mcm mcc er hi d;
stoch_simul(order=1, irf=40)ymt ymt_gap pi pji CH;
//model_diagnostics(M_,options_,oo_)
y_pos=strmatch('ymt',var_list_ ,'exact');
y_gap_pos=strmatch('ymt_gap',var_list_ ,'exact');
pi_pos=strmatch('pi',var_list_ ,'exact');
pji_pos=strmatch('pji',var_list_ ,'exact');
CH_pos=strmatch('CH',var_list_ ,'exact');



    %----------------------------------------------------------------
    % Traditional money demand shock
    %----------------------------------------------------------------

    %phi_pi_vec=[1.5 1.5 5 1.5];
    %phi_y_vec=[0.125 0 0 1];
     phi_pi_estimated = get_param_by_name('rhopi');
    phi_y_estimated = get_param_by_name('rhoy');
    phi_pi_vec = phi_pi_estimated + [-0.2, 0, 0.2];  % Adjust range as needed
    phi_y_vec = phi_y_estimated + [-0.1, 0, 0.1];    % Adjust range as needed

shocks;
        var epsx=epsx; 
        var epsa=0; %see description p. 113
    end;   

    variance.y_gap=NaN(1,length(phi_pi_vec));
    variance.y=NaN(1,length(phi_pi_vec));
    variance.pi=NaN(1,length(phi_pi_vec));
    variance.pji=NaN(1,length(phi_pi_vec));
    variance.CH=NaN(1,length(phi_pi_vec));
    
    L=NaN(1,length(phi_pi_vec));
    for ii=1:length(phi_pi_vec)
        set_param_value('rhopi',phi_pi_vec(ii));
        set_param_value('rhoy',phi_y_vec(ii));
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %loop over stoch_simul
        if ~info(1)
            %read out current parameter values
            par.theta1=M_.params(strmatch('theta1',M_.param_names,'exact'));
            par.alpha=M_.params(strmatch('alpha',M_.param_names,'exact'));
            par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
            par.eta=M_.params(strmatch('eta',M_.param_names,'exact'));
            par.sigma=M_.params(strmatch('sigma',M_.param_names,'exact'));
            par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
            par.theta2=M_.params(strmatch('theta1',M_.param_names,'exact'));
            par.zeta=M_.params(strmatch('zeta',M_.param_names,'exact'));
            par.lambda=(1-par.theta1)*(1-par.betta*par.theta1)/par.theta1*(1-par.alpha)/(1-par.alpha+par.alpha*par.eta);
            par.lambda2=(1-par.theta2)*(1-par.betta*par.theta2)/par.theta2*(1-par.alpha)/(1-par.alpha+par.alpha*par.eta);

            variance.y_gap(ii)=oo_.var(y_gap_pos,y_gap_pos);
            variance.y(ii)=oo_.var(y_pos,y_pos);
            variance.pi(ii)=oo_.var(pi_pos,pi_pos);
            variance.pji(ii)=oo_.var(pji_pos,pji_pos);
            variance.CH(ii)=oo_.var(CH_pos,CH_pos);
            
         % L(ii)=0.5*((par.sigma)*variance.y_gap(ii)+par.eta/par.lambda*variance.pi(ii))/100;
         L(ii)=0.5*((par.sigma)*(variance.y_gap(ii))+((par.varphi+par.alpha)/(1-par.alpha)*(1+ksi)*variance.CH(ii))+par.eta/par.lambda*variance.pi(ii)+par.eta/par.lambda2*(variance.pji(ii)))/100;
            end
    end
    %Print result

    labels={'rhopi';'rhoy';'sigma(y)';'sigma(tilde y)';'sigma(pi)';'L'};
    headers={' ';' ';' ';' '};
    values=[phi_pi_vec;phi_y_vec;sqrt(variance.y);sqrt(variance.y_gap);sqrt(variance.pi);L];
    options_.noprint=0;
    dyntable(options_,'Traditional Money',headers,labels,values,size(labels,2)+2,4,3)
    options_.noprint=1;

    %----------------------------------------------------------------
    % Cryptocurrency Demand shock
    %----------------------------------------------------------------
    shocks;
        var epsdj=epsdj; 
        var epsa=0; %see description p. 113
    end;   

    variance.y_gap=NaN(1,length(phi_pi_vec));
    variance.y=NaN(1,length(phi_pi_vec));
    variance.pi=NaN(1,length(phi_pi_vec));
    variance.pji=NaN(1,length(phi_pi_vec));
    variance.CH=NaN(1,length(phi_pi_vec));
  
    L=NaN(1,length(phi_pi_vec));
    for ii=1:length(phi_pi_vec)
        set_param_value('rhopi',phi_pi_vec(ii));
        set_param_value('rhoy',phi_y_vec(ii));
        [info, oo_, options_] = stoch_simul(M_, options_, oo_, var_list_); %loop over stoch_simul
        if ~info(1)
            %read out current parameter values
            par.theta1=M_.params(strmatch('theta1',M_.param_names,'exact'));
            par.alpha=M_.params(strmatch('alpha',M_.param_names,'exact'));
            par.betta=M_.params(strmatch('betta',M_.param_names,'exact'));
            par.eta=M_.params(strmatch('eta',M_.param_names,'exact'));
            par.sigma=M_.params(strmatch('sigma',M_.param_names,'exact'));
            par.varphi=M_.params(strmatch('varphi',M_.param_names,'exact'));
             par.theta2=M_.params(strmatch('theta1',M_.param_names,'exact'));
            par.zeta=M_.params(strmatch('zeta',M_.param_names,'exact'));
            par.lambda=(1-par.theta1)*(1-par.betta*par.theta1)/par.theta1*(1-par.alpha)/(1-par.alpha+par.alpha*par.eta);
            par.lambda2=(1-par.theta2)*(1-par.betta*par.theta2)/par.theta2*(1-par.alpha)/(1-par.alpha+par.alpha*par.eta);

            variance.y_gap(ii)=oo_.var(y_gap_pos,y_gap_pos);
            variance.y(ii)=oo_.var(y_pos,y_pos);
            variance.pi(ii)=oo_.var(pi_pos,pi_pos);
            variance.pji(ii)=oo_.var(pji_pos,pji_pos);
            variance.CH(ii)=oo_.var(CH_pos,CH_pos);
         
        %  L(ii)=0.5*((par.sigma)*variance.y_gap(ii)+par.eta/par.lambda*variance.pi(ii))/100;
         L(ii)=0.5*((par.sigma)*(variance.y_gap(ii))+((par.varphi+par.alpha)/(1-par.alpha)*(1+ksi)*variance.CH(ii))+par.eta/par.lambda*variance.pi(ii)+par.eta/par.lambda2*(variance.pji(ii)))/100;

            end
    end
    %Print result
    values=[phi_pi_vec;phi_y_vec;sqrt(variance.y);sqrt(variance.y_gap);sqrt(variance.pi);L];
    options_.noprint=0;
    dyntable(options_,'Cryptocurrency',headers,labels,values,size(labels,2)+2,4,3);