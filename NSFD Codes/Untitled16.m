T=1.0; % simulation time
format long
Fo=0.18;
f=Fo;
g=1+(4*Fo);
temp=100;
%THIS SECTION CALCULATES TEMPs FOR N=11
ms=11; %matrix size
dxa=1/(ms-1); %length of dx, same as dy
dta=Fo*dxa^2;
ta=(dta:dta:T); %sample times
xa=(0:dxa:1); %length
nta=length(ta); %time marching length
ma=ms^2; %sizing of the matrix
na=ms^2;
TPa(ma,1)=0;
TP1a(ma,nta)=0;
aa(1:ma,1:na)=0; %declaration of the sparse matrix
for e=(1:ma)
    aa(e,e)=1; %inserting 1's on the diagonal elements
end
for ee=(1:(ms-2))
    for e=(2:(ms-1))
        qa=(ms*ee)+e;
        TPa(qa)=temp;
        aa(qa,qa)=g;
        aa(qa,qa-1)=-f;
        aa(qa,qa+1)=-f;
        aa(qa,qa-ms)=-f;
        aa(qa,qa+ms)=-f;
    end
end
aa=inv(aa); % Inversion of matrix aa
TPa(((ma-ms)+1):ma)=temp;
%below is the time-marching calculation
TP1a(:,1)=TPa;
for d=(2:nta)
    TP1a(:,d)=aa*TPa;
    TPa=TP1a(:,d);
end
z1a=(ma-ms)/2+1; %declaring the range for the x-center line nodes
z2a=(ma+ms)/2;
z3a(1:ms)=(ms+1)/2; %declaring the range for the y-center line nodes
for k=(2:ms)
    z3a(k)=z3a(1)+((k-1)*ms);
end
ca(ms)=0;
for k=(1:ms)
    ca(k)=TP1a(z3a(k),ma);
end
t11a=ta; %n=11
TP111a=TP1a(((ma+1)/2),:); %n=11
x11b=xa; %n=11
TP111b=TP1a(z1a:z2a,ma); %n=11
c11c=ca; %n=11
figure(1),plot(t11a,TP111a,'-.')
grid on
hold on
figure(2),plot(x11b,TP111b,'-.')
grid on
hold on
figure(3),plot(x11b,c11c,'-.')
grid on
hold on
%THIS SECTION CALCULATES TEMPs FOR N=21
ms=21; %matrix size
dxb=1/(ms-1); %length of dx, same as dy
dtb=Fo*dxb^2;
tb=(dtb:dtb:T); %sample times
xb=(0:dxb:1); %length
ntb=length(tb); %time marching length
mb=ms^2; %sizing of the matrix
nb=ms^2;
TPb(mb,1)=0;
TP1b(mb,ntb)=0;
ab(1:mb,1:nb)=0; %declaration of the sparse matrix
for e=(1:mb)
    ab(e,e)=1; %inserting 1's on the diagonal elements
end
for ee=(1:(ms-2))
    for e=(2:(ms-1))
        qb=(ms*ee)+e;
        TPb(qb)=temp;
        ab(qb,qb)=g;
        ab(qb,qb-1)=-f;
        ab(qb,qb+1)=-f;
        ab(qb,qb-ms)=-f;
        ab(qb,qb+ms)=-f;
    end
end
ab=inv(ab); % Inversion of matrix ab
TPb(((mb-ms)+1):mb)=temp;
%below is the time-marching calculation
TP1b(:,1)=TPb;
for d=(2:ntb)
    TP1b(:,d)=ab*TPb;
    TPb=TP1b(:,d);
end
z1b=(mb-ms)/2+1; %declaring the range for the x-center line nodes
z2b=(mb+ms)/2;
z3b(1:ms)=(ms+1)/2; %declaring the range for the y-center line nodes
for k=(2:ms)
    z3b(k)=z3b(1)+((k-1)*ms);
end
cb(ms)=0;
for k=(1:ms)
    cb(k)=TP1b(z3b(k),mb);
end
t21a=tb; %n=21
TP121a=TP1b(((mb+1)/2),:); %n=21
x21b=xb; %n=21
TP121b=TP1b(z1b:z2b,mb); %n=21
c21c=cb; %n=21
figure(1),plot(t21a,TP121a,'r-.')
grid on
figure(2),plot(x21b,TP121b,'r-.')
grid on
figure(3),plot(x21b,c21c,'r-.')
grid on
%THIS SECTION CALCULATES TEMPs FOR N=41
ms=41; %matrix size
dxc=1/(ms-1); %length of dx, same as dy
dtc=Fo*dxc^2;
tc=(dtc:dtc:T); %sample times
xc=(0:dxc:1); %length
ntc=length(tc); %time marching length
mc=ms^2; %sizing of the matrix
nc=ms^2;
12
TPc(mc,1)=0;
TP1c(mc,ntc)=0;
ac(1:mc,1:nc)=0; %declaration of the sparse matrix
for e=(1:mc)
    ac(e,e)=1; %inserting 1's on the diagonal elements
end
for ee=(1:(ms-2))
    for e=(2:(ms-1))
        qc=(ms*ee)+e;
        TPc(qc)=temp;
        ac(qc,qc)=g;
        ac(qc,qc-1)=-f;
        ac(qc,qc+1)=-f;
        ac(qc,qc-ms)=-f;
        ac(qc,qc+ms)=-f;
    end
end
ac=inv(ac); % Inversion matrix of ac
TPc(((mc-ms)+1):mc)=temp;
%below is the time-marching calculation
TP1c(:,1)=TPc;
for d=(2:ntc)
    TP1c(:,d)=ac*TPc;
    TPc=TP1c(:,d);
end
z1c=(mc-ms)/2+1; %declaring the range for the x-center line nodes
z2c=(mc+ms)/2;
z3c(1:ms)=(ms+1)/2; %declaring the range for the y-center line nodes
for k=(2:ms)
    z3c(k)=z3c(1)+((k-1)*ms);
end
cc(ms)=0;
for k=(1:ms)
    cc(k)=TP1c(z3c(k),mc);
end
t41a=tc; %n=41
TP141a=TP1c(((mc+1)/2),:); %n=41
x41b=xc; %n=41
TP141b=TP1c(z1c:z2c,mc); %n=41
c41c=cc; %n=41
figure(1),plot(t41a,TP141a,'g-.')
grid on
xlabel('time (s) ')
ylabel('nodal temperature (C)')
title('Temp. Profile for center node T(5,5) (Implicit Scheme)')
legend('11x11 nodes','21x21 nodes','41x41 nodes')
figure(2),plot(x41b,TP141b,'g-.')
grid on
xlabel('length (x-axis) ')
ylabel('nodal temperature (C)')
title('Temp. Profile for the x-axis center nodes (x=0.5) (Implicit Scheme)')
legend('11x11 nodes','21x21 node','41x41 nodes')
figure(3),plot(x41b,c41c,'g-.')
grid on
xlabel('length (y-axis) ')
ylabel('nodal temperature (C)')
title('Temp. Profile for the y-axis center nodes (y=0.5) (Implicit Scheme)')
legend('11x11 nodes','21x21 nodes','41x41 nodes')