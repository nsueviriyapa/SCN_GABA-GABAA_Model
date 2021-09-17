clear all; close all; clc;  


%% Set up GABA/GABAA conditions

xaa_t          = [ 1.00  1.00   1.00   1.00    1.00];                %Changing only delta GABAA receptor expression levels
xaa_p          = [ 1.00  1.00   1.00   1.00    1.00];                %Changing only gamma GABAA receptor expression levels
xaa_receptor   = [ 1.00  0.75   0.50   0.25    0.00];                %Changing both detla & gammma GABAA receptors (KD_t & KD_p) simultaneously 

xaa_t2          = [ 1  1  1  1  1 ];                  % Controlling tonic GABA signaling pathway
xaa_p2          = [ 1  1  1  1  1 ];                  % Controlling phasic GABA signaling pathway
xaa_signaling   = [ 1  1  1  1  1 ];                  % Adjusting for GABA blockade vs. application

xaa_VIP         = [ 1  1  1  1  1];                   % Controlling overall VIP coupling activities 

xaa_scenario    = [ 0  25  50 75  100];               % Scenario: different varitions in GABAA receptor expression or GABA signaling levels

           
%% Run simulations
           

for w= 1:length(xaa_scenario)           %Number of scenarios
for ii= 1:5                             %Number of independent runs for network realization
    

tic



%%%%%%%% Basics %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


options = odeset('RelTol',1e-3,'AbsTol',1e-6);


ncell_core=   10*10; 
ncell_shell=  12*12;  
ncell = ncell_core+ncell_shell;


ns = 24;%# of states



%%%%%%%% Create Heterogeneity %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 sd = 0.01;
 
 vsP0=(0.94)*ones(ncell,1)+sqrt(((0.94)*sd*6)^2)*randn(ncell,1);
 vsB = ((1).*ones(ncell,1))+ sqrt(((1).*sd*(1))^2)*randn(ncell,1);
 vmB = ((0.8)*ones(ncell,1))+ sqrt(((0.8)*sd*(1))^2)*randn(ncell,1);
 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% GABA and VIP %%%%%%%%%%%%%%%


SF_VIP=1*xaa_VIP(w); 

SF_GABA_t= 1*(xaa_signaling(w).^xaa_t2(w));
SF_GABA_p= 1*(xaa_signaling(w).^xaa_p2(w));

KD_t=1.*(xaa_t(w)).*(xaa_receptor(w));
KD_p=1.*(xaa_p(w)).*(xaa_receptor(w));

vClo_Pm = 5.85.*ones(ncell,1); 
vClo_Tm = 13.5;

Vspill_m=0.05;
v_tGABA=0.4;
v_pGABA=0.684.*ones(ncell,1);
K_pGABA=22.*ones(ncell,1);
n_pGABA=(3.4).*ones(ncell,1);

Cl_core  =0.4;  
Cl_shell =0.4;   

vClo_core  = (Cl_core*ones(ncell_core,1)); 
vClo_shell = (Cl_shell*ones(ncell_shell,1)); 

vClo(1:ncell_core,1)=vClo_core;
vClo(ncell_core+1:ncell,1)=vClo_shell;

gPT=0.8; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  Connectivity Matrix %%%%%%%%%%%%%%%

% CORE %
Perc_VIP =0.5 ;  
bita=0.05;      
[A A1 VIP_prod]=adjacency_core(ncell_core,bita,Perc_VIP);

% SHELL %
A_shell=adjacency_shell(ncell_shell);

% ADD LINKS FROM CORE-SHELL%
bitacs=0.05; 
[A_vip A_gaba]=coretoshell(bitacs,ncell_core,ncell_shell,A, A1,A_shell,VIP_prod);

sumal_vip=1./sum(A_vip,2)';
sumal_vip(isinf(sumal_vip))=0; 

sumal_pgaba=1./sum(A_gaba,2)'; 
sumal_pgaba(isinf(sumal_pgaba))=0; 

sumal_tgaba=ones(1,ncell)./ncell;


clear A A1 A_shell bitacs VIP_prod 



%%%%%% Solve ODEs %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



tic

t = [0:0.1:600];
[t,y]=ode23(@ODEs,t,SCN_IC(ncell,ncell_core,ns),options,SCN_Param,A_vip,sumal_vip,A_gaba,sumal_pgaba,ncell,ns,vsP0,vsB,vmB,vClo,v_pGABA,K_pGABA,n_pGABA,sumal_tgaba,vClo_Pm,vClo_Tm,KD_t,KD_p,v_tGABA,Vspill_m,gPT,SF_VIP,SF_GABA_t,SF_GABA_p);

clear options sd  bita bita_core b r 

toc 





end
end


                 

