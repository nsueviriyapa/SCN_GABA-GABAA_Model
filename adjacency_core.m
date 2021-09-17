function [A A1 VIP_prod]=adjacency_core(ncell,bita,Perc_VIP)

kk=1:ncell;
X=kk(ones(1,ncell),:);
s=sqrt(ncell);
kk=kk';
Var1=fix((kk-1)/s);
Var1=Var1(:,ones(1,ncell));
Var2=rem((kk-1),s);
Var2=Var2(:,ones(1,ncell));
A = (sqrt((Var1-fix((X-1)/s)).^2+(Var2-rem((X-1),s)).^2)).\1;
A(A~=1)=0;

clear Var1 Var2 

%%%% CREATE PERIODIC BOUNDARY CONDITIONS %%%%

% FOLD THE TOPS

for i=1:s   
A(i,s*(s-1)+i)=1;
A(s*(s-1)+i, i)=1;
end

% FOLD SIDES

w=1:s:ncell; 
for o=1:length(w)
A(w(o),(s-1)+w(o))=1;
A((s-1)+w(o),w(o))=1;
end
%%%%
 
% EXCLUDE SELF CONNECTIONS

for i=1:ncell
    A(i,i)=NaN;
end


% ADD CONNECTIONS RANDOMLY

for i=1:ncell
    r=X(i,A(i,:)==0);   
    b=rand(length(r),1);
    b(b<bita)=0;
    r=r(b==0);
    if ~isempty(r)
        for k=1:length(r)
            A(i,r(k))=1;  
            A(r(k),i)=1;
        end
    end
end

A(isnan(A))=0; 

%%%% GABAergic NETWORK %%%%

A1=A; 
clear b r k i w X 


%%%% VIP PERCENTAGE %%%%

b=rand(1,ncell);       
b(b <Perc_VIP)=0;     
b(b~=0)=1;
wa=kk(b==1);

VIP_prod=kk(b==0);    

    for o = 1:length(wa)
        A(:,wa(o))=0;    
    end
    
clear Perc_VIP b wa s kk
    
