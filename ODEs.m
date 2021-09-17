function dydt= ODEs(t,y,p,A_vip,sumal_vip,A_gaba,sumal_pgaba,ncell,ns,vsP0,vsB,vmB,vClo,v_pGABA,K_pGABA,n_pGABA,sumal_tgaba,vClo_Pm,vClo_Tm,KD_t,KD_p,v_tGABA,Vspill_m,gPT,SF_VIP,SF_GABA_t,SF_GABA_p)

dydt=zeros(ns*ncell,1);

i=1:ncell;

Ca=y((i-1)*ns+1);
Ca_store =  y((i-1)*ns+2);
MP = y((i-1)*ns+3);
MC = y((i-1)*ns+4);
MB = y((i-1)*ns+5);
PC = y((i-1)*ns+6);
CC = y((i-1)*ns+7);
PCP = y((i-1)*ns+8);
CCP = y((i-1)*ns+9);
PCC = y((i-1)*ns+10);
PCN = y((i-1)*ns+11);
PCCP = y((i-1)*ns+12);
PCNP = y((i-1)*ns+13);
BC = y((i-1)*ns+14);
BCP = y((i-1)*ns+15);
BN = y((i-1)*ns+16);
BNP = y((i-1)*ns+17);
IN = y((i-1)*ns+18);
CB = y((i-1)*ns+19);
vVIP = y((i-1)*ns+20);
pGABA=y((i-1)*ns+21);
tGABA=y((i-1)*ns+22);
Cl_in=y((i-1)*ns+23);
osc=y((i-1)*ns+24);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%PARAMETERS%%%%%%%%%%%%%%%%%%%%%

k1=  p(1);
k2 = p(2);
k3 = p(3);
k4 = p(4);
k5 = p(5);
k6 = p(6);
k7 = p(7);
k8 = p(8);
KAP = p(9);
KAC = p(10);
KIB = p(11);
kdmb = p(12);
kdmc = p(13);
kdmp = p(14);
kdnc = p(15);
kdn = p(16);
Kd = p(17);
Kdp = p(18);
Kp = p(19);
KmB = p(20);
KmC = p(21);
KmP = p(22);
ksB = p(23);
ksC = p(24);
ksP = p(25);
n = p(26);
m = p(27);
Vphos = p(28);
V1P = Vphos;
V1PC = Vphos;
V3PC = Vphos;
V1B = p(29);
V1C = p(30);
V2B = p(33);
V2C = p(34);
V2P = p(35);
V2PC = p(36);
V3B = p(37);
V4B = p(39);
V4PC = p(40);
vdBC = p(41);
vdBN = p(42);
vdCC = p(43);
vdIN = p(44);
vdPC = p(45);
vdPCC = p(46);
vdPCN = p(47);
%vmB= p(48);%
vmC = p(49);
vmP = p(50);
%vsB = p(51);
vsC = p(52);
KD = p(55);%
b_IP3=p(58);
vP = p(59);
VMK=p(60);
K_1 = p(62);
K_2 = p(63);
WT = p(64);
CT = p(65);
KC = p(66);
%vsP0 = p(67);%
kf=p(69);
IP3=p(70);
VM3= p(71);
M3= p(72);
KR= p(73);
KA = p(74);
pA = p(75);
VM2=p(76);
K2=p(77);
M2= p(78);
kMK=p(79);
V_b=p(80);
k_b=p(81);
GABA_o=p(82);
v_kk=p(83);
K_kk=p(84);
n_kk=p(85);
v_VIP=p(86);
K_VIP=p(87);
n_VIP=p(88);
k_dVIP=p(89);
n_dVIP=p(90);
%v_GABA=p(91);
% K_GABA=p(92);
% n_GABA=p(93);
k_dpGABA=p(94);
n_dpGABA=p(95);
v_glu = p(96);
K_glu =p(97);
K_gluR =p(98);
v_gluR=p(99);
% vClo =p(100);
vCl1 =p(101);
nCl1 =p(102);
KCl1 =p(103);
%vCl2 =p(104);
%nCl2 =p(105);
%KCl2 =p(106);
vCl3 =p(107);
nCl3 =p(108);
KCl3 =p(109);
nClout =p(110);
%v_tGABA=p(111);
K_tGABA=p(112);
n_tGABA=p(113);
k_dtGABA=p(114);
n_dtGABA=p(115);
V_th=p(116);
K_th=p(117);



%%%%%%%%% Neurotransmitters %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if t<150
    
    VKO=0;  
    switz=0;
    switz2=0;
    beta1=zeros(ncell,1);
    S_pGABA=zeros(ncell,1); 
    S_tGABA=zeros(ncell,1); 
  


elseif t>150 && t<250 
    
    VKO=1.*SF_VIP;  
    switz=1;
    switz2=0;
    
    %%%%%%% VIP calcs %%%%%%%
    
    VIP=vVIP';
    VIP=VIP(ones(1,ncell),:);
    S_VIP=sum(VIP.*A_vip,2)'.*sumal_vip;
    clear VIP
    beta1 = (0+S_VIP)'./(KD+(0+S_VIP))';

    %%%%%%% GABA calcs %%%%%%%

    GABAp=pGABA';   
    GABAp=GABAp(ones(1,ncell),:);
    S_pGABA=sum(GABAp.*A_gaba,2)'.*sumal_pgaba; 

    GABAt=tGABA';    
    GABAt=GABAt(ones(1,ncell),:);
    A_tgaba=ones(ncell,ncell);
    S_tGABA= sum(GABAt.*A_tgaba,2)'.*sumal_tgaba  ; 


    
else
    
    VKO=1.*SF_VIP;  
    switz=1;
    switz2=1;
    
    %%%%%%% VIP calcs %%%%%%%
    
    VIP=vVIP';
    VIP=VIP(ones(1,ncell),:);
    S_VIP=sum(VIP.*A_vip,2)'.*sumal_vip;
    clear VIP
    beta1 = (0+S_VIP)'./(KD+(0+S_VIP))';

    %%%%%%% GABA calcs %%%%%%%

    GABAp=pGABA';    
    GABAp=GABAp(ones(1,ncell),:);
    S_pGABA=sum(GABAp.*A_gaba,2)'.*sumal_pgaba; 

    GABAt=tGABA';    
    GABAt=GABAt(ones(1,ncell),:);
    A_tgaba=ones(ncell,ncell);
    S_tGABA= sum(GABAt.*A_tgaba,2)'.*sumal_tgaba  ; 

    clear GABA
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                                                                       
glu=16*MP./(4+MP);
b_gluR=glu./(2+glu); 


b_tGABA=switz.*(S_tGABA).*(SF_GABA_t.^switz2);
b_pGABA=switz.*(S_pGABA).*(SF_GABA_p.^switz2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[vv E_inhib I_inhib  I_gClo  I_g_pGABA  I_g_tGABA  IPSC  Vrest theta g_inhib g_pGABA g_tGABA gbasal ]= FiringRates((Ca), Fir,CC,BC,MP,Cl_in,osc,KD_t,KD_p,gPT,SF_VIP,SF_GABA_t,SF_GABA_p,switz2,b_tGABA,b_pGABA);



VGCCm= (0.0525).*MB.^0.5./(1+MB.^0.5);
VGCC=  VGCCm.*vv.^1.9./(20+vv.^1.9) ;

vv3 =(VM3.*(Ca_store.^M3)./(KR^M3+Ca_store.^M3)).*(Ca.^pA)./(KA^pA+Ca.^pA);
vv2 =VM2 *(Ca.^M2)./(K2^M2+Ca.^M2);
kk1=v_kk*CC.^n_kk./(K_kk+CC.^n_kk);
vo= b_gluR.*0.04;



nCa  =(1.1.*osc) + (1.2625.*(1-osc));
nVMK =(1.1.*osc) + (1.2625.*(1-osc));

nbeta1=(1.1.*osc) + (1.2625.*(1-osc));
nV_b  =(1.1.*osc) + (1.2625.*(1-osc));

nCT   =((1.5).*osc) + ((1.5).*(1-osc));
nvsP0 =((0.985).*osc) + ((0.985).*(1-osc));


vK = switz*( ( (  (VMK.*nVMK)   ).*(Ca.^nCa)./((kMK.^nCa)+(Ca.^nCa)))+ (VKO.*( (V_b.*nV_b)  ).*(beta1.^nbeta1)./((k_b.^nbeta1)+(beta1.^nbeta1))) );%

vsPc =nvsP0.*vsP0  + nCT.*CT.*CB./(KC+CB);


nClo_P=1.25;
nClo_T=1.25;


vClo_P =((vClo_Pm.*(KD_p.^switz2)) ).*((MP.^nClo_P./(1+MP.^nClo_P)).^osc);
vClo_T = ((vClo_Tm.*(KD_t.^switz2))).*((CC.^nClo_T./(20+CC.^nClo_T)).^osc);



Vspill = Vspill_m*(MP.^2.5./(1+MP.^2.5));



Vup_Am=0.1;
n_up_A=1;
K_up_A=1;
Vup_A = Vup_Am*MP.^n_up_A/(K_up_A+MP.^n_up_A);




for i =1:ncell
    
    
    dydt((i-1)*ns+1,1) =  [VGCC(i)+vo(i)+b_IP3*IP3-vv2(i)+vv3(i)+kf*Ca_store(i) - kk1(i)*Ca(i).^(2)];
    dydt((i-1)*ns+2,1) =  [(vv2(i)-vv3(i)-kf*Ca_store(i))]; 
    dydt((i-1)*ns+3,1) = vsPc(i)*BN(i)^n/(KAP^n+BN(i)^n)-vmP*MP(i)/(KmP+MP(i))-kdmp*MP(i);
    dydt((i-1)*ns+4,1) = vsC*BN(i)^n/(KAC^n+BN(i)^n)-vmC*MC(i)/(KmC+MC(i))-kdmc*MC(i);
    dydt((i-1)*ns+5,1) = vsB(i)*KIB^m/(KIB^m+BN(i)^m)-vmB(i)*MB(i)/(KmB+MB(i))-kdmb*MB(i);
    dydt((i-1)*ns+6,1) = ksP*MP(i)-V1P*PC(i)/(Kp+PC(i))+V2P*PCP(i)/(Kdp+PCP(i))+k4*PCC(i)-k3*PC(i)*CC(i)-kdn*PC(i);
    dydt((i-1)*ns+7,1) = ksC*MC(i)-V1C*CC(i)/(Kp+CC(i))+V2C*CCP(i)/(Kdp+CCP(i))+k4*PCC(i)-k3*PC(i)*CC(i)-kdnc*CC(i);
    dydt((i-1)*ns+8,1) = V1P*PC(i)/(Kp+PC(i))-V2P*PCP(i)/(Kdp+PCP(i))-vdPC*PCP(i)/(Kd+PCP(i))-kdn*PCP(i);
    dydt((i-1)*ns+9,1) = V1C*CC(i)/(Kp+CC(i))-V2C*CCP(i)/(Kdp+CCP(i))-vdCC*CCP(i)/(Kd+CCP(i))-kdn*CCP(i);
    dydt((i-1)*ns+10,1) = -V1PC*PCC(i)/(Kp+PCC(i))+V2PC*PCCP(i)/(Kdp+PCCP(i))-k4*PCC(i)+k3*PC(i)*CC(i)+k2*PCN(i)-k1*PCC(i)-kdn*PCC(i);
    dydt((i-1)*ns+11,1) = -V3PC*PCN(i)/(Kp+PCN(i))+V4PC*PCNP(i)/(Kdp+PCNP(i))-k2*PCN(i)+k1*PCC(i)-k7*BN(i)*PCN(i)+k8*IN(i)-kdn*PCN(i);
    dydt((i-1)*ns+12,1) = V1PC*PCC(i)/(Kp+PCC(i))-V2PC*PCCP(i)/(Kdp+PCCP(i))-vdPCC*PCCP(i)/(Kd+PCCP(i))-kdn*PCCP(i);
    dydt((i-1)*ns+13,1) = V3PC*PCN(i)/(Kp+PCN(i))-V4PC*PCNP(i)/(Kdp+PCNP(i))-vdPCN*PCNP(i)/(Kd+PCNP(i))-kdn*PCNP(i);
    dydt((i-1)*ns+14,1) = ksB*MB(i)-V1B*BC(i)/(Kp+BC(i))+V2B*BCP(i)/(Kdp+BCP(i))-k5*BC(i)+k6*BN(i)-kdn*BC(i);
    dydt((i-1)*ns+15,1) = V1B*BC(i)/(Kp+BC(i))-V2B*BCP(i)/(Kdp+BCP(i))-vdBC*BCP(i)/(Kd+BCP(i))-kdn*BCP(i);
    dydt((i-1)*ns+16,1) = -V3B*BN(i)/(Kp+BN(i))+V4B*BNP(i)/(Kdp+BNP(i))+k5*BC(i)-k6*BN(i)-k7*BN(i)*PCN(i)+k8*IN(i)-kdn*BN(i);
    dydt((i-1)*ns+17,1) = V3B*BN(i)/(Kp+BN(i))-V4B*BNP(i)/(Kdp+BNP(i))-vdBN*BNP(i)/(Kd+BNP(i))-kdn*BNP(i);
    dydt((i-1)*ns+18,1) = -k8*IN(i)+k7*BN(i)*PCN(i)-vdIN*IN(i)/(Kd+IN(i))-kdn*IN(i);
    dydt((i-1)*ns+19,1) = (vP/WT)*(vK(i)/vP*(1-CB(i))/(K_1+1-CB(i))-CB(i)/(K_2+CB(i)));
    dydt((i-1)*ns+20,1) = v_VIP*vv(i).^n_VIP/(K_VIP+vv(i).^n_VIP)-k_dVIP*vVIP(i).^n_dVIP;
    dydt((i-1)*ns+21,1) = v_pGABA(i)*vv(i).^n_pGABA(i)/(K_pGABA(i)+vv(i).^n_pGABA(i))- Vspill(i).*pGABA(i)*[1./[1+ exp((theta(i)-V_th)/K_th)]] -k_dpGABA*pGABA(i).^n_dpGABA;
    dydt((i-1)*ns+22,1) = v_tGABA*1000*Ca(i).^n_tGABA/(K_tGABA+1000*Ca(i).^n_tGABA)+ 1*Vspill(i).*pGABA(i)*[1./[1+ exp((theta(i)-V_th)/K_th)]] -(k_dtGABA+ Vup_A(i))*tGABA(i).^n_dtGABA;
    dydt((i-1)*ns+23,1) = [vClo(i)]+ switz.*(  [ 0.1*0.5*vClo_T(i) +   ( 0.2*0.5*vClo_T(i)- (0.1*0.5)*vClo_T(i)).*(1-  ((b_tGABA(i).^3)./(0.01^3 +  (b_tGABA(i).^3))) )]       +    [  (0.125)*vClo_P(i)  + ( (0.4)*vClo_P(i) - (0.125)*vClo_P(i)).*( 1-  ((b_pGABA(i).^5)./(1.25^5 +  (b_pGABA(i).^5)) ) ) ]     ) + (vCl1 )*(MP(i).^nCl1./(KCl1+MP(i).^nCl1)) -((0.2*0.5)+(0.5*vCl3 )*(MP(i).^nCl3./(KCl3+MP(i).^nCl3)))*Cl_in(i).^nClout; %mM 
    dydt((i-1)*ns+24,1) = 0;
   
 
    
end




