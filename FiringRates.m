function [v E_inhib I_inhib  I_gClo  I_g_pGABA  I_g_tGABA  IPSC  Vrest theta g_inhib  g_pGABA g_tGABA gbasal ]= FiringRates(Ca_in,F,CC,BC,MP,Cl_in,osc,KD_t,KD_p,gPT,SF_VIP,SF_GABA_t,SF_GABA_p,switz2,b_tGABA,b_pGABA)

% PARAMETERS

Ca_out = F(1);
Cl_ex=F(2);
% Cl_o=F(3);
Cp=F(4);
E_ex=F(5);
EK_o=F(6);
EL_o=F(7);
ENa_o=F(8);
% g_inhib= F(10);
gKo=F(11);
gNa =F(12);
k = F(13);
K_R=F(14);
KCa= F(15);
% KCl1=F(16);
% KCl2=F(17);
Kex1=F(18);
Kex2=F(19);
Kgk=F(20);
KKCa=F(21);
nca=F(24);
% nCl=F(25);
nex1=F(26);
nex2=F(27);
nKCa=F(28);
q = F(33);
T = F(35);
T_room= F(36);
V_R=F(37);
V_theta=F(38);
vCa=F(39);
% vCl1=F(40);
% vCl2=F(41);
Vex1= F(42);
Vex2=F(43);
vgk=F(44);
vKCa=F(45);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REVERSAL POTENTIALS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ENa = ENa_o * T/(T_room);%
EK = (EK_o)*T/(T_room); %
EL =EL_o* T/(T_room); 
ECa =k*T/(2*q)*log(Ca_out./Ca_in)*1000; %
E_inhib = -k*T/(q)*log(Cl_ex./(Cl_in+0))*1000';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEMBRANE PROPERTIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vrest =properties (Ca_in,ENa,EK,q,k,T,Ca_out,Cl_ex,Cl_in,F,BC); % 
theta =Vrest + V_theta ; 
Vreset= Vrest+4;
R=V_R*(Vrest)./(K_R + Vrest);  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONDUCTANCES & CURRENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

INa = gNa.*(Vrest-ENa);

gK=(gKo+MP./(Kgk+MP)*vgk);
g_ex=(Vex1*abs(INa).^nex1./(Kex1+abs(INa).^nex1)+ ((Ca_in).^nex2)./(Kex2+(Ca_in).^nex2).*Vex2);
gL = 1./R ;
gCa=vCa*(MP.^nca./(KCa+MP.^nca));
gKCa =vKCa*(CC.^nKCa./(KKCa+CC.^nKCa));

nClo_P=1.25;
nClo_T=1.25;  

expP=(MP.^nClo_P./(1+MP.^nClo_P)).^osc;
expT=(CC.^nClo_T./(20+CC.^nClo_T)).^osc;

vgPT=39.7782; 
vgP=gPT*vgPT.*(KD_p.^switz2);
vgP_r=vgP.*(expP.^1);
g_pGABA  = 0.1.*vgP_r   +  0.05.*vgP_r  + (0.1.*vgP_r  -  0.05.*vgP_r ).*(1-((b_pGABA.^(2))./(1.25^(2) +  (b_pGABA.^(2)))));

vgT=(1-gPT)*vgPT.*(KD_t.^switz2);
vgT_r=vgT.*(expT.^1);
g_tGABA  = 0.5.*vgT_r   +  0.1.*vgT_r   + (0.6*vgT_r  -  0.1.*vgT_r)   .*(1-((b_tGABA.^(3))./(0.01^(3) +  (b_tGABA.^(3)))));

gClo= 7.92;
g_inhib =(gClo + g_pGABA+g_tGABA);
gbasal=g_inhib - g_pGABA-g_tGABA; %This is gClo

I_gClo = gClo.*(Vrest- E_inhib);
I_g_pGABA = g_pGABA.*(Vrest- E_inhib);
I_g_tGABA = g_tGABA.*(Vrest- E_inhib);
IPSC = I_gClo+ I_g_pGABA+I_g_tGABA;
I_inhib = g_inhib.*(Vrest- E_inhib);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINAL CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I_star = (-g_inhib.*E_inhib -g_ex.*E_ex+gNa.*ENa + gCa.*ECa + gK.*EK + gL.*EL+ gKCa.*EK) ;
R_star = 1./(gNa+ gK + gL + gCa + gKCa- g_inhib - g_ex) ;
tau_m =Cp.*R_star;

v= -(tau_m.*log( (theta-(R_star.*I_star))./(Vreset-R_star.*I_star))).^-1 ;

