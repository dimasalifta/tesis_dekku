Khw=5.81;
gr=4;
%AlphacKhwgr=10000.3720
%AlphacKhwgr=0.3720
%AlphacKhwgr=1.3778 
AlphacKhwgr=0; 
%AlphavKhwgr=0.114
%AlphavKhwgr=0.214
AlphavKhwgr=0.914;

JrFL=0.0675; %Full Load
JrHL=0.0404; %Half Load
JrNL=0.0133;  %No Load

%JrFL=0.0675 %Full Load
%JrHL=0.0675+0.0271 %Half Load
%JrNL=0.0675+0.0542  %No Load

JrC1=0.0675; %(mwd=0.8 , rwd=0.05, mwl=2.0, rwd=0.05)
JrC2=0.1033; %(mwd=2.0 , rwd=0.05, mwl=0.8, rwd=0.05)
JrC3=0.1114; %(mwd=2.0 , rwd=0.05, mwl=0,   rwd=0.05)

%JrC1=0.0675 %Full Load
%JrC2=0.0404 %Half Load
%JrC3=0.0133  %No Load

K=Khw*gr/AlphavKhwgr;
tau=JrC1/AlphavKhwgr;

%Ao1=[0 1  ; 0 -1/0.5919]
%B=[0 ; 203.8596/0.5919]
%C=[1 0]

Ao1=[0 1  ; 0 -1/tau];
B=[0 ; K/tau];
C=[1 0];

%BACKLASH=20000;
%BACKLASH=200;
%BACKLASH=50;
BACKLASH=0;
%%BACKLASH=50000;

c=10;

%syms s

%PID=tf([0.2 0.3 0.1],[1 0])
G2=tf((1/(Khw*gr))*(JrC1)*(c^2),1);
%G2=tf([-(1/(Khw*gr))*(JrC1)*(c) 0 0],[1 0])
G1=tf([(1/(Khw*gr))*(JrC1) 0 0],1);



[PLANTnum, PLANTden]=ss2tf(Ao1,B,C,0,1);
PLANT=tf(PLANTnum, PLANTden);
PLANT2=tf(PLANTnum*(0.0675), PLANTden*(0.0675));
 
wcf=2;
q=tf(2*pi*wcf,[1 (2*pi*wcf)]) ;

H=PLANT*(G1-G2)/(1-G2*PLANT);
H=minreal(H);
Gff=1.9;
HFLOWPASSFILTER=tf(2*pi*Gff,[1 (2*pi*Gff)]);
%HFLOWPASSFILTER=tf([(2*pi*Gff)^2],[1 2*(2*pi*Gff) (2*pi*Gff)^2])
%Gf=(1/H)*HFLOWPASSFILTER;
%Gf=minreal(Gf);
Gf=tf(1,1);
[numGf,denGf]=tfdata(Gf,'v');


s=tf('s');



HPOLE=pole(H);
HSTABILITY=isstable(H);

QTIMEONEMINUSGFH=q*(1-Gf*H);
NORM=norm(q*(1-Gf*H),Inf)  ;
bode(q,inv(1-Gf*H)) ;
%QTIMEINVONEPLUSG=q*inv(1+G)
%NORM=norm(q*inv(1+G),Inf)  
%bode(q,1+G) 

E=[-1/0.07 0 (0.6287/0.07)-(-8); 0 -6 0-(-15); -1/(0.0271)  -1/(0.0271)  0];
lamda1=eig(E);        %symetric
lamda2=eig((E+E')/2);  %non-symetric
 

%E=[-1/3 0 (0.206/3)-(-200); 0 -600 0-(-1500); -1/(0.017)  -1/(0.017)  0]
%lamda1=eig(E)        %symetric
%lamda2=eig((E+E')/2)  %non-symetric
 [R, P] = chol(E);     
 
sim('ELMNNSTRCECP',40);