clear
clc
close all

global n delta2 c Pmax
n=7;
delta2=10^(-6); Pmax=3;
In=diag(ones(1,n));
A=[0 1 1 0 1 1 0;
   1 0 1 0 0 1 0;
   1 1 0 1 1 0 0;
   0 0 1 0 1 0 1;
   1 0 1 1 0 0 0;
   1 1 0 0 0 0 1;
   0 0 0 1 0 1 0];
D=sum(A,2);
L=diag(D)-A;
BarL=kron(L,In);
BarA=kron(A,In);
R=diag(ones(1,n));
BarR=blkdiag(R(1,:),R(2,:),R(3,:),R(4,:),R(5,:),R(6,:),R(7,:));


R1=[0  4  4  6  5  6  8;
    4  0  7  9  8  4  6;
    4  7  0  3  3  7  6;
    6  9  3  0  3  6  4;
    5  8  3  3  0  9  7;
    6  4  7  6  9  0  3;
    8  6  6  4  7  3  0];


R1=10^6*R1;
pi=3.14; D_T=0.15; D_R=0.3; theta_T=.2; theta_R=.2; lambda=1.55*10^(-6);
for i=1:n
    for j=1:n
        if j~=i
            g(i,j)=(pi^2*D_T^2*D_R^2)/(16*lambda^2*R1(i,j)^2)...
                *exp(D_T^2 *theta_T^2+D_R^2*theta_R);
        end
    end
end
c=0.1*ones(1,n);  b=1*[1 3 2 1 2 3 1 ]; %

Bard=A*D;

kappa1=3.7;kappa2=0.2;zeta=10*ones(n,2); beta=80; delta=0.4;

load ('My_AttackDefence.mat')

x=x0;
chex=x;
N0=4000;
xend=[2.2474;2.4878;2.5916;1.9764;1.5995;1.8826;0.9101]; %Optimal NE points
xo=kron(xend, ones(1,N0+1));
for k=1:N0
    mu(k)=0.05/((k+1)^0.01);
    muk=mu(k);
    if k==1
        chex=x;
        for i=1:n
            Trs(i,k)=1;
        end
    else
        for i=1:n
            
            temp=BarL((i-1)*n+1:i*n,:);
            temp2=temp*chex(:,k-1);
            omega(i,k)=temp2'*temp2;
            Chi_min=min(1+muk^2,2*beta*muk^2*(1/delta-1));
            xi(i,k)=(2*beta*muk^2*(1/delta-1)-Chi_min)*(D(i)^3+D(i)^2+D(i)+Bard(i))...
                +(1+muk^2-Chi_min)*Bard(i)+Chi_min*D(i)^2;
            epshatxi(i,k)=beta*(1-delta)*muk^2*omega(i,k)/xi(i,k);
            
        end
        for i=1:n
            hatxi=x((i-1)*n+1:i*n,k);
            chexik_1=chex((i-1)*n+1:i*n,k-1);
            exik=chexik_1-hatxi;
            
            PsiValue(i,k)=kappa1*(exik'*exik-epshatxi(i,k))-zeta(i,k);
            if PsiValue(i,k)>=0
                chex((i-1)*n+1:i*n,k)=hatxi;
                Trs(i,k)=1;
            else
                chex((i-1)*n+1:i*n,k)=chex((i-1)*n+1:i*n,k-1);
                Trs(i,k)=0;
            end
            eik=chex((i-1)*n+1:i*n,k)-hatxi;
            zeta(i,k+1)=(1-kappa2)*zeta(i,k)+epshatxi(i,k)-eik'*eik;
            eabs2(i,k)=eik'*eik;
            Threshold(i,k)=epshatxi(i,k)+zeta(i,k)/kappa1;
        end
        
    end
    e(:,k)=chex(:,k)-x(:,k);
    x(:,k+1)=x(:,k)-mu(k)*(BarR'*BarFFun(x(:,k),g,A,b)+BarL*x(:,k)+BarA*e(:,k));
 
end
xx=BarR*x;
barI=[In In In In In In In ];
Tildex=barI*x/n-xx; 
TrsRatio=(sum(Trs,2)/N0*100)'




figure1=figure;
axes1 = axes('Parent',figure1);
hold on
for i=1:1:n
    plot(xo(i,:),'--','linewidth',2)
    plot(xx(i,:),'linewidth',2)
end
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes1,'FontSize',12,'FontName','Times New Roman');
h3=xlabel('Iteration (k) ');
set(h3,'Interpreter','none','FontSize',12,'FontName','Times New Roman');
h4=ylabel('${x}_{i,k}$ ');
set(h4,'Interpreter','latex','FontSize',14)
xlim([0 N0])




figure1=figure;
axes1 = axes('Parent',figure1);
hold on
for i=1:1:n
    plot(Tildex(i,:),'linewidth',2)
end
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes1,'FontSize',12,'FontName','Times New Roman');
h3=xlabel('Iteration (k) ');
set(h3,'Interpreter','none','FontSize',12,'FontName','Times New Roman');
h4=ylabel('$\vec{x}_{k}$ ');
set(h4,'Interpreter','latex','FontSize',14)
xlim([0 N0])


figure1=figure;
axes1 = axes('Parent',figure1);
plot(eabs2(1,:),'--','linewidth',2)
hold on
plot(Threshold(1,:),'-','linewidth',2)
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes1,'FontSize',12,'FontName','Times New Roman');
h3=xlabel('Iteration (k) ');
set(h3,'Interpreter','none','FontSize',12,'FontName','Times New Roman');
h4=ylabel(' Amplitude ');
set(h4,'Interpreter','none','FontSize',12,'FontName','Times New Roman');
legend1=legend('$\|{e}_{1,k}\|^2$ ','$\epsilon_{\hat{x}_k^1}+\zeta_{1,k}/\kappa_1$ ');
set(legend1,'Interpreter','latex','FontSize',14)
xlim([0 100])

axes2=axes('Position',[0.42 0.28 0.28 0.25]);
plot(eabs2(1,80:100),'--','linewidth',2)
hold on
plot(Threshold(1,80:100),'-','linewidth',2)
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes2,'XTickLabel',{'80','85','90','95','100'},...
    'XTick',[1 6 11 16 21],'FontSize',12,'FontName','Times New Roman');
xlim([1 21])




figure1=figure;
axes1 = axes('Parent',figure1);
hold on
for i=1:1:n
    plot(zeta(i,:),'linewidth',2)
end
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes1,'FontSize',12,'FontName','Times New Roman');
h3=xlabel('Iteration (k) ');
set(h3,'Interpreter','none','FontSize',12,'FontName','Times New Roman');
h4=ylabel('$\zeta_{i,k}$');
set(h4,'Interpreter','latex','FontSize',14)
xlim([0 100])
legend1=legend('$\zeta_{1,k}$','$\zeta_{2,k}$','$\zeta_{3,k}$','$\zeta_{4,k}$','$\zeta_{5,k}$','$\zeta_{6,k}$','$\zeta_{7,k}$');
set(legend1,'Interpreter','latex','FontSize',12)
axes2=axes('Position',[0.42 0.28 0.28 0.25]);
hold on
for i=1:1:n
    plot(zeta(i,80:100),'linewidth',2)
end
grid on
set(gca,'GridLineStyle',':','GridColor','k','GridAlpha',1);
box on
set(axes2,'XTickLabel',{'80','85','90','95','100'},...
    'XTick',[1 6 11 16 21],'FontSize',12,'FontName','Times New Roman');
xlim([1 21])
ylim([-0.002 0.015])



%%  Functions
function output=Fenzi(p,g,i,j)
global n delta2
Sum=0;
for s=1:n
    if (s~=i) &&(s~=j)
        Sum=Sum+g(s,j)*p(s);
    end
end
Sum=Sum+delta2;
output=g(i,j)/Sum;
end

%---------------------------------
function Fi=FFun(p,g,A,i,b)
global c n Pmax delta2
Temp=A(i,1:4);
pi=p((i-1)*n+1:i*n,1);
Label=find(Temp==1);
Sum=0;
for s=1:size(Label,2)
    j=Label(s);
    temp1=g(i,j)/( g(:,j)'*pi+delta2);
    Sum=Sum+temp1;
end
Fi=c(i)+.1*(-pi(i)^(-2)+(Pmax-sum(pi(i)))^(-2))-Sum*b(i);%
end

function Fi=FFunAttack(p,g,A,i,b)
global c n Pmax 

Temp=A(i,1:4);
pi=p((i-1)*n+1:i*n,1);
Label=find(Temp==1);
Sum=0;
for s=1:size(Label,2)
    j=Label(s);
    Sum=Sum+FFunAttackSub(p,g,A,i,j);
end

Fi=c(i)+.1*(-pi(i)^(-2)+(Pmax-sum(pi(i)))^(-2))-Sum*b(i);%
end

function Fi=FFunAttackSub(p,g,A,i,j)
global  n  delta2

Temp=A(j,1:4);
pi=p((j-1)*n+1:j*n,1);
Label=find(Temp==1);
Sum=0;
for s=1:size(Label,2)
    m=Label(s);
    temp1=g(i,s)*g(j,s)*pi(j)/(( g(:,m)'*pi+delta2)*(g(:,m)'*pi+delta2-g(j,s)*pi(j)));
    Sum=Sum+temp1;
end

Fi=Sum;
end

function F=BarFFun(p,g,A,b)
global n

for i=1:n
    
    if i<=4
        F(i,1)=FFun(p,g,A,i,b);
    else
        F(i,1)=FFunAttack(p,g,A,i,b);
    end
end
end


%}

