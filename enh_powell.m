%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Engin Aybey
% Enhanced Powell's Method
% Requiring funcpwll.m,goldsec.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=2;
a=-10.0;
b=10.0;
delta=10.0^-6;%tolerance for golden search
maxiter=4;
X=[0.99;1.1];% our estimation value
U=zeros(n,n);
%creating U matrix
for k=1:1:n
  E=zeros(1,n);
  E(k)=E(k)+1;
  U(:,k)=U(:,k)+E';
end
pM=zeros(n,n);
i=0;
r=0;
while i<maxiter
  P=X;
  pM(:,1)=P; %P_0=X
  for k=1:1:n
    %f=x^2-4*x+y^2-y-x*y
    f='( (P(1,1)) + (gama) * (U(1,k)) )^2 - ( (P(1,1)) + (gama) * (U(1,k)) )*4 + ( (P(2,1)) + (gama) * (U(2,k)) )^2 - ( (P(2,1)) + (gama) * (U(2,k)) ) - ( (P(1,1)) + (gama) * (U(1,k)) )*( (P(2,1)) + (gama) * (U(2,k)) )';
    f = strrep(f,'P(1,1)',num2str(P(1)));
    f = strrep(f,'P(2,1)',num2str(P(2)));
    f = strrep(f,'U(1,k)',num2str(U(1,k)));
    f = strrep(f,'U(2,k)',num2str(U(2,k)));
    func1 = inline(f,'gama');
    gama=goldsec(func1,a,b,delta); % gama that minimizes f(P_{k-1}+gama*U_k)
    pM(:,k)=P; %pM holds the columns of P
    P=P+gama*U(:,k);
    pM(:,k+1)=P;
    del_f(k)=funcpwll(pM(:,k+1))-funcpwll(pM(:,k)); %delta f_k=f(P_k)-f(P_{k-1})
    clear gama;
  end
  
  del_fr=max(abs(del_f)); % magnitude of the max decrease f
  %finding subscript r for using in U_r
  for j=1:1:n
    if del_fr==abs(del_f(j))
       r=j;
    end
  end
  %increment the counter
  i=i+1;
  %f_k=f(P_k) for k=0,1,2,...,n
  for k=0:1:n
    fnc(k+1)=funcpwll(pM(:,k+1));
  end
  %the function vale in the extended direction 2(P_n-P_0) from P_0 
  fncE=funcpwll(2*pM(:,n+1)-pM(:,1));
  if fncE>fnc(1)
    X=pM(:,n+1); %X_i=P_n
    continue;
  elseif 2*(fnc(1)-2*fnc(n+1)+fncE)*(fnc(1)-fnc(n+1)-del_fr)^2>=del_f*(fnc(1)-fncE)^2
    X=pM(:,n+1); %X_i=P_n
    continue;
  else
    break;
  end
  %finding U_r
  U(:,r)=pM(:,n+1)-pM(:,1);
  %f2=x^2-4*x+y^2-y-x*y
  f2='( (P_0(1,1)) + (gama1) * (U(1,2)) )^2-( (P_0(1,1)) + (gama1) * (U(1,2)) )*4+ ( (P_0(2,1)) + (gama1) * (U(2,2)) )^2-( (P_0(2,1)) + (gama1) * (U(2,2)) )-( (P_0(1,1)) + (gama1) * (U(1,2)) )*( (P_0(2,1)) + (gama1) * (U(2,2)) )';
  f2 = strrep(f2,'P_0(1,1)',num2str(pM(:,1)));
  f2 = strrep(f2,'P_0(2,1)',num2str(PM(:,2)));
  f2 = strrep(f2,'U(1,2)',num2str(U(1,2)));
  f2 = strrep(f2,'U(2,2)',num2str(U(2,2)));
  func2 = inline(f2,'gama1');

  gama1=goldsec(func2,a,b,delta);% gama1 that minimizes f(P_0+gama1*U_r)
  
  X=pM(:,1)+gama1*U(:,r);%X_i=P_0+gama1*U_r
  clear gama1;
end
X
min_value_of_function=funcpwll(X)
