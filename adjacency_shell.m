function [A ]=adjacency_shell(ncell)

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

%CREATE PERIODIC BOUNDARY CONDITIONS


%FOLD THE TOPS

for i=1:s   
A(i,s*(s-1)+i)=1;
A(s*(s-1)+i, i)=1;
end


%FOLD SIDES

w=1:s:ncell; 
for o=1:length(w)
A(w(o),(s-1)+w(o))=1;
A((s-1)+w(o),w(o))=1;
end

 
return



