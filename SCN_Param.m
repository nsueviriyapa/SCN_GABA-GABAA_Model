function params = SCN_Param

k1 =0.3875;
k2 = 0.1;
k3 = 0.4;
k4 = 0.2;
k5 = 0.4;
k6 = 0.2;
k7 = 0.5;
k8 = 0.1;
KAP =0.6;
KAC = 0.6;
KIB = 2.2;
kdmb = 0.01;
kdmc = 0.01;
kdmp = 0.01;
kdnc = 0.12;
kdn = 0.01;
Kd = 0.3;
Kdp = 0.1;
Kp = 0.1;
KmB = 0.4;
KmC = 0.4;
KmP = 0.31;
ksB = 0.12;
ksC = 1.6;
ksP = 0.6;
n = 4;
m = 2;
Vphos = 0.4;
V1B = 0.5;
V1C = 0.66; 
V1P = Vphos;        
V1PC = Vphos;
V2B = 0.1;
V2C = 0.1;
V2P = 0.3;
V2PC = 0.1;
V3B = 0.5;
V3PC = Vphos;
V4B = 0.2;
V4PC = 0.1;
vdBC = 0.5;
vdBN = 0.6;
vdCC = 0.7;
vdIN = 0.8;
vdPC = 0.7;
vdPCC = 0.7;
vdPCN = 0.7;
vmB = 0.8;
vmC = 1.0;
vmP = 1.1;
vsB = 1;
vsC = 1.1;
vsP = 1.5;
RT = 1;
KD = 0.06;
k = 10;
v0 = 0.5;
b_IP3=0.0003;
vP = 1;
VMK =5.175; 
Ka = 2.5;
K_1 = 0.01;
K_2 = 0.01;
WT = 1;
CT =1.8;
KC =0.15;
vsP0 = 1;
v2= 5;
kf = 0.001;
IP3=0.5;
VM3=400;
M3=6;
KR=3;
KA = 0.67;
pA = 4.2;
VM2=149.5;
K2=5;
M2= 2.2;
kMK=4.2;
V_b= 1.8 ;
K_b=1.8;


GABA_o=0.00001;
v_kk=3.3; 
K_kk=.02;
n_kk=0.1;
v_VIP=0.375; 
K_VIP=20;
n_VIP=1.9;
k_dVIP=0.5;
n_dVIP=0.2;
v_pGABA=0.5;
K_pGABA=20;
n_pGABA=1.9;

k_dpGABA=0.5;
n_dpGABA=0.7;

v_glu = 10;
K_glu =4;
K_gluR =2;
v_gluR=0.1;

vClo =0.1;
 
vCl1 =13.95; 
nCl1 =1;
KCl1 =1;
 
vCl2 =20; 
nCl2 =1; 
KCl2 =1;
 
vCl3 = 21.6;
nCl3 =1; 
KCl3 =20; 
nClout =1; 

v_tGABA=0.4;
K_tGABA=3.5;
n_tGABA=0.5;
k_dtGABA=1;
n_dtGABA=0.2;
V_th=1;
K_th=1;


params = [k1 k2 k3 k4 k5 k6 k7 k8 KAP KAC KIB kdmb kdmc kdmp kdnc kdn Kd Kdp Kp KmB KmC KmP ksB,...
ksC ksP n m Vphos V1B V1C V1P V1PC V2B V2C V2P V2PC V3B V3PC V4B V4PC vdBC vdBN vdCC vdIN vdPC vdPCC...
vdPCN vmB vmC vmP vsB vsC vsP,...
RT KD k v0 b_IP3 vP VMK Ka K_1 K_2 WT CT KC vsP0 v2 kf IP3 VM3 M3 KR KA pA VM2 K2 M2 kMK V_b K_b...
GABA_o v_kk K_kk n_kk v_VIP K_VIP n_VIP k_dVIP n_dVIP v_pGABA K_pGABA n_pGABA k_dpGABA n_dpGABA...
v_glu K_glu  K_gluR v_gluR,...
vClo,vCl1,nCl1,KCl1,vCl2,nCl2,KCl2,vCl3,nCl3,KCl3,nClout,v_tGABA,K_tGABA,n_tGABA,k_dtGABA,n_dtGABA,V_th,K_th];
