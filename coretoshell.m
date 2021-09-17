function [A_vip A_gaba]=coretoshell(bita,ncell_core,ncell_shell,A, A1,A_shell,VIP_prod)

ncell=ncell_core+ncell_shell;



% VIP matrix

A_vip=zeros(ncell,ncell);
A_vip(1:ncell_core,1:ncell_core)=A;



% GABA matrix

A_gaba=zeros(ncell,ncell);
A_gaba(1:ncell_core,1:ncell_core)=A1;
A_gaba(ncell_core+1:ncell,ncell_core+1:ncell)=A_shell;



% CREATE CONNECTIONS FROM CORE TO SHELL

for i=ncell_core+1:ncell
    r=1:length(VIP_prod);   
    b=rand(length(r),1);   
    b(b<bita)=0;
    r=r(b==0);
    
    if ~isempty(r)
        for k=1:length(r)
            A_vip(i,VIP_prod(r(k)))=1;  
            A_gaba(i,VIP_prod(r(k)))=1;
            A_gaba(VIP_prod(r(k)),i)=1; 
        end
    end
end
