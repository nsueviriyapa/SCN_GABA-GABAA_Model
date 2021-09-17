function IC2 = SCN_IC(ncell,ncell_core,ns)

IC1(1) =0.1;
IC1(2)= 0.1;
IC1(3) = 2.79566;
IC1(4) = 1.964254;
IC1(5) = 7.946372;
IC1(6) = .4;
IC1(7) = 12;
IC1(8) = 0.126363;
IC1(9) = 8.986265;
IC1(10) = 1.25558;
IC1(11) = 0.165471;
IC1(12) = 0.197619;
IC1(13) = 0.090808;
IC1(14) = 2.408975;
IC1(15) = 0.479535;
IC1(16) = 1.941518;
IC1(17) = 0.32736;
IC1(18) = 0.049602;
IC1(19)= 0.12;
IC1(20)=0;
IC1(21)=1;
IC1(22)=0.001;
IC1(23)=0.1;


IC2=zeros(1,ncell*ns);


for i = 1:ncell_core
    for j = 1:23
        IC2((i-1)*ns+j) = IC1(j); 
    end
    
    for j = 24
        IC2((i-1)*ns+j) = 1;        
    end
end



for i = ncell_core+1:ncell
    for j = 1:23
        IC2((i-1)*ns+j) = IC1(j); 
    end
    
    for j = 24
        IC2((i-1)*ns+j) = 0;        
    end
end

 

clear y IC1 sol temp options step lags

