clear all
close all
clc
format longeng



smax=25;
phi=pi/3;
sig1=0.33;
sig2=8.3*(10^(-3));
sig3=0.33;
k1=sig1/sig2;
k2=sig2/sig3;
f=6.8*10^-2;
a=0.25*10^-2;
b=7.36*10^-2;    % b=R1      
c=7.76*10^-2;    % c=R2 
R=8*10^-2;       % R=R3
rr=9*10^-2;
Roh0=a/c;
Roh3=R/c;
for j=1:smax
    
    C1a(j)=(1-k1)*(b^(2*j+1));
    C1b(j)=((k1*j)+j+1).*(b^(2*j+1));
    
%     C2a(j)=((1-k1)*((k2*j)+j+k2).*(b^(4*j+2)))+((1-k2)*((k1*j)+j+1).*(b^(2*j+1)).*(c^(2*j+1)));
%     C2b(j)=((1-k1)*(1-k2)*(+j.^2+j).*(b^(4*j+2)))+(((k1*j)+j+1).*((k2*j)+j+1).*(b^(2*j+1)).*(c^(2*j+1)));
    dd1=[((k2*j+j+k2).*(b^(2*j+1))) ((1-k2)*(c^(2*j+1))); ((1-k2)*(j.^2+j).*(b^(2*j+1))) ((k2*j+j+1).*(c^(2*j+1)))];
    dd2=[0 ((1-k1)*(b^(2*j+1))); 0 ((k1*j+j+1).*(b^(2*j+1)))];
    dd3=[1;1];
    dd=dd1*dd2*dd3;
    C2a(j)=dd(1);
    C2b(j)=dd(2);

    m11(j)=(R^(2*j+1)).*((k2*j)+j+1);
    m12(j)=(k2-1)*(c^(2*j+1)).*(j+1);
    m21(j)=(k2-1)*(R^(2*j+1)).*j;
    m22(j)=(c^(2*j+1)).*((k2*j)+j+k2);
    
    mp11(j)=(Roh3.^(j+1)).*((k2*j)+j+1);
    mp12(j)=(k2-1)*(j+1)./(Roh3.^j);
    mp21(j)=(k2-1)*(Roh3.^(j+1)).*j;
    mp22(j)=(k2*j+j+k2)./(Roh3.^j);
    
end




%******************* Eccentric Four Layer ************************************
sig1=0.33;
sig2=1.8;
sig3=8.3*(10^(-3));
sig4=0.33;
k3=sig3/sig4;
R1=7.2*10^-2;        
R2=7.36*10^-2;     
R3=7.76*10^-2;
R4=8*10^-2;

for j=1:smax
    
    t3=[((k2*j+j+1).*(R3^(2*j+1))) ((k2-1)*(j+1)*(R2^(2*j+1))); ((k2-1)*(j).*(R3^(2*j+1))) ((k2*j+j+k2).*(R2^(2*j+1)))];
    t4=[((k3*j+j+1).*(R4^(2*j+1))) ((k3-1)*(j+1)*(R3^(2*j+1))); ((k3-1)*(j).*(R4^(2*j+1))) ((k3*j+j+k3).*(R3^(2*j+1)))];
    
    dd=t3*t4;
    
    m11(j)=dd(1);
    m12(j)=dd(3);
    m21(j)=dd(2);
    m22(j)=dd(4);
end

for s=1:smax
    for j=1:s
        G(s,j)=(j.*(R1^(2*j+1)).*C1a(j))./(C1b(j).*(a^(j+1)).*(factorial(j-1)).*(factorial(s-j)));
    end
    zy(s)=((-1)*(s.^2).*(m11(s)+m12(s)))./(((factorial(s+1)).^2).*(m21(s)+m22(s)).*(Roh0^(2*s+1)));
    for j=s:smax
        H(s,j)=zy(s).*(((-1)^(j-s)).*factorial(j+1).*(a^j))./(factorial(j-s));
    end
end            
P=G+H;

for s=1:smax
    Qp=0;
    for j=1:s
    Qp=Qp+((-1)*(f^(j-1).*(2*j+1).*(R1^(2*j+1)))/(C1b(j).*(a^(j+1)).*(factorial(j-1)).*factorial(s-j)));
    end
    Qpp(s)=Qp;
end

Q=transpose(Qpp);
Z=(P^-1)*Q;
ZT=transpose(Z);

m=1;
for t=1:180
    teta=t*pi/180;
    vr=0;
    vphi=0;
    vteta=0;
    rp=sqrt(rr^2+a^2-(2*a*rr*cos(teta)));
    for s=1:smax
        
        pp1=legendre(s,(cos(teta)));
        pp2=legendre((s+1),(cos(teta)));
        p1(s)=pp1(2);
        p2(s)=pp2(2);
        
        zzr(s)=(a*sin(phi))*(f^(s-1)+((R1^(2*s+1))*ZT(s)))/(((k1+1)*s+1)*rr^(s+1));
        vr=vr+(zzr(s)*p1(s));
        
        zzphi(s)=(cos(phi))*(f^(s-1)+((R1^(2*s+1))*ZT(s)))/(((k1+1)*s+1)*rr^(s+1)*sin(teta));
        q(s)=((s+1)*cos(teta)*p1(s))-(s*p2(s));
%         q(s)=((s+1)*((rr*cos(teta)-a)/rp(t))*p1(s))-(s*p2(s));
        vphi=vphi+(zzphi(s)*q(s));
        
        zzteta(s)=(sin(phi))*((rr-a*cos(teta))/(rr*sin(teta)))*(f^(s-1)+((R1^(2*s+1))*ZT(s)))/(((k1+1)*s+1)*rr^(s+1));
        vteta=vteta+(zzteta(s)*p1(s));
                      
    end
    vtot=sqrt((vr^2)+(vphi^2)+(vteta)^2);
    v3(m)=vtot;
    m=m+1;
end
m=1:179;
re=(k1-1)*((0.4*10^-16)/3);
%--------------------------------

load('datax.mat');
m=1:32400;
for s=1:180
    m=((s-1)*180)+1;
    x(s)=data(m);
end
s=1:179;
r=(8*10^-14)/3;
t=1:179;
err(t)=(abs(re*v3(t)-r*x(t))./abs(re*v3(t)))*10;
u=1:179;
errs(u)=err(180-u);

plot(u,errs(u))
xlabel('\theta (deg)')
ylabel('Relative error (%)')
