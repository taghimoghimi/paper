clear all
close all
clc
format longeng
nmax=70;
sr1=0.495; 
% st1=0.495; 
st1=0.2475;
sr2=0.33; st2=0.33;
sr3=1.8; st3=1.8;
sr4=0.00498;
% st4=0.00498;
st4=0.00996;
sr5=0.33; st5=0.33;
r1=0.052; r2=0.072; r3=0.0736; r4=0.0776; r5=0.08;
rr=0.09;
g=0.5*10^-2;
rQ=0.068;
Q=1e-10;
mu0=4*pi*(1*10^-7);
phi=pi/3;
phiQ=pi/3;

%**************************************
for i=1:nmax
    b15(i)=1;
    a15(i)=((i+1)/i)*(r5^(-2*i-1));
end
%**************************************

%**************************************
for i=1:nmax
    v4=sqrt((i*(i+1)*(st4/sr4))+0.25)-0.5;
    z1=r4^(v4); z2=r4^(-1*(v4+1)); z3=sr4*v4*r4^(v4-1); z4=(-1*sr4)*(v4+1)*r4^(-1*(v4+2));
    y1=a15(i)*(r4^(i))+r4^(-1*(i+1)); y2=sr5*i*a15(i)*(r4^(i-1))-sr5*(i+1)*r4^(-1*(i+2));
              
    a14(i)=((z4*y1)-(z2*y2))/((z1*z4)-(z2*z3));
    b14(i)=((z1*y2)-(z3*y1))/((z1*z4)-(z2*z3));
end
%**************************************

%**************************************
for i=nmax:-1:1
    v4=sqrt((i*(i+1)*(st4/sr4))+0.25)-0.5;
    x(i)=(((((i+1)*sr3+(i*sr4))*r3^(v4-1)*a14(i))+((sr3-sr4)*(i+1)*r3^(-1*(v4+2))*b14(i)))*((g^i)*factorial(i+1)))/(sr3*(2*i+1)*r3^(i-1));
    w=0;
    for j=(i+1):nmax
        w=w+(((-1)^(j-i))*(factorial(j+1))*(g^j)*a13(j))/(factorial(j-i));
    end
    a13(i)=(x(i)-w)/(g^i*factorial(i+1));
end

for i=1:nmax
    y(i)=((((i/(2*i+1))*(1-(sr4/sr3))*a14(i)*(r3^v4))+(1/(2*i+1))*(i+((i+1)*(sr4/sr3)))*b14(i)*(r3^(-1*(v4+1))))*(r3^(i+1)))/((g^(i+1))*(factorial(i-1)));
    z=0;
    for j=1:(i-1)
        z=z+(b13(j)/((g^(j+1))*(factorial(j-1))*factorial(i-j)));
    end
    b13(i)=(y(i)-z)*((g^(i+1))*factorial(i-1));
end
%**************************************

%**************************************
for i=1:nmax
    z1=r2^(i); z2=r2^(-1*(i+1)); z3=sr2*i*r2^(i-1); z4=(-1*sr2)*(i+1)*r2^(-1*(i+2));
    y1=a13(i)*(r2^(i))+b13(i)*(r2^(-1*(i+1))); y2=sr3*i*a13(i)*(r2^(i-1))-sr3*(i+1)*b13(i)*r2^(-1*(i+2));
        
        
    a12(i)=((z4*y1)-(z2*y2))/((z1*z4)-(z2*z3));
    b12(i)=((z1*y2)-(z3*y1))/((z1*z4)-(z2*z3));
end
%**************************************

%**************************************
for i=1:nmax
    v1=sqrt((i*(i+1)*(st1/sr1))+0.25)-0.5;
    z1=r1^(i); z2=r1^(-1*(i+1)); z3=sr2*i*r1^(i-1); z4=(-1*sr2)*(i+1)*r1^(-1*(i+2));
    y1=r1^(v1); y2=sr1*v1*(r1^(v1-1));
    
    
    a02(i)=((z4*y1)-(z2*y2))/((z1*z4)-(z2*z3));
    b02(i)=((z1*y2)-(z3*y1))/((z1*z4)-(z2*z3));
end
%**************************************

%**************************************
for i=1:nmax
    v1=sqrt((i*(i+1)*(st1/sr1))+0.25)-0.5;
    z1=r1^(v1); z2=r1^(-1*(v1+1)); z3=sr1*v1*r1^(v1-1); z4=(-1*sr1)*(v1+1)*r1^(-1*(v1+2));
    y1=r1^(i); y2=r1^(-1*(i+1)); y3=sr2*i*r1^(i-1); y4=(-1*sr2)*(i+1)*r1^(-1*(i+2));
        
    b11(i)=(((z1*y3-z3*y1)*a12(i))+((z1*y4-z3*y2)*b12(i)))/(z1*z4-z2*z3);
end
%**************************************

%**************************************
t=16.5*pi/180;
w=cos(t);
for i=1:nmax    
    pp1=legendre(i,w);

    p1(i)=pp1(2);
    
end


%**************************************

%**************************************

m=1;
for t=1:179
    teta=t*pi/180;
    vteta_a=0;
    vteta_b=0;
    vphi_a=0;
    vphi_b=0;
    for n=1:nmax
        pp1=legendre(n,(cos(teta)));
        pp2=legendre((n+1),(cos(teta)));
        pq1(n)=pp1(2);
        pq2(n)=pp2(2);
        
        Q3a(n)=Q*cos(phiQ)*(sqrt(2*n+1))*(sqrt(n*(n+1)/(8*pi)))*(rQ^(n-1)*(a02(n))+rQ^(-1*(n+2))*(b02(n)))/((2*n+1)*sr2*(a02(n)*b12(n)-a12(n)*b02(n)));
        Q3b(n)=(-1)*Q*sin(phiQ)*(sqrt(2*n+1))*(sqrt(1/(8*pi*n*(n+1))))*(rQ^(n-1)*(a02(n))+rQ^(-1*(n+2))*(b02(n)))/((2*n+1)*sr2*(a02(n)*b12(n)-a12(n)*b02(n)));
        

        
        
        zza(n)=(Q3a(n)*(b13(n)+((r2^(2*n+1))*a13(n)))*sqrt(((2*n+1)*factorial(n-1))/(2*pi*factorial(n+1)))*pq1(n))/((2*n+1)*rr^(n+1));
        zzb(n)=(Q3b(n)*(b13(n)+((r2^(2*n+1))*a13(n)))*sqrt(((2*n+1)*factorial(n-1))/(2*pi*factorial(n+1)))*pq1(n))/((2*n+1)*rr^(n+1));
        yya(n)=(Q3a(n)*(b13(n)+((r2^(2*n+1))*a13(n)))*sqrt(((2*n+1)*factorial(n-1))/(2*pi*factorial(n+1)))*((n+1)*cos(teta)*pq1(n)-n*pq2(n)))/((2*n+1)*rr^(n+1));
        yyb(n)=(Q3b(n)*(b13(n)+((r2^(2*n+1))*a13(n)))*sqrt(((2*n+1)*factorial(n-1))/(2*pi*factorial(n+1)))*((n+1)*cos(teta)*pq1(n)-n*pq2(n)))/((2*n+1)*rr^(n+1));
        vteta_a=vteta_a+zza(n);
        vteta_b=vteta_b+zzb(n);
        vphi_a=vphi_a+yya(n);
        vphi_b=vphi_b+yyb(n);
        
    end
    Bteta_a=(vteta_a*sin(phi))/sin(teta);
    Bteta_b=(vteta_b*cos(phi))/sin(teta);
    Bteta=Bteta_a+Bteta_b;
    Bphi_a=(vphi_a*cos(phi))/sin(teta);
    Bphi_b=(vphi_b*sin(phi))/sin(teta);
    Bphi=Bphi_a+Bphi_b;
    vvv=sqrt((Bphi^2)+(Bteta^2));
    Btot(m)=vvv;

    m=m+1;
end
rrr=(-1)*mu0*(sr2-sr3);
m=1:179;

plot(m,rrr*Btot(m))
hold on
load('datax.mat');
m=1:32400;
for s=1:180
    m=((s-1)*180)+1;
    x(s)=data(m);
end
s=1:179;
r=1.3*10^-14;
plot(s,r*x(s))
xlabel('\theta (deg)')
ylabel('Magnetic field (T)')
legend('Analytical model','Numerical model (FEM)')

