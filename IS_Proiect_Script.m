%Prelucrare date primite

t=scope37(:,1);
u=scope37(:,2);
y1=scope37(:,3);
y2=scope37(:,4);

figure;
plot(t,u,t,y1,t,y2);title('Setul de date experimentale "scope37"');legend('u','y1','y2');

figure;
subplot(311);
plot(t,u);title('Semnalul u');legend('u');
subplot(312);
plot(t,y1,'red');title('Semnalul y1');legend('y1');

subplot(313);
plot(t,y2,'yellow');title('Semnalul y2');legend('y2');

%%
% Partea 1
% Identificarea neparametrica a semnalului y1

t=scope37(:,1);
u=scope37(:,2);
y1=scope37(:,3);
y2=scope37(:,4);

plot(t,u,t,y1);

t_u_min=492; 
t_u_max=480;
t_y1_min=496;
t_y1_max=485;

Mr1=(y1(t_y1_max)-y1(t_y1_min))/(u(t_u_max)-u(t_u_min)); 
zeta=sqrt(2-sqrt(4-4/Mr1^2))/2; 

K=mean(y1)/mean(u); 

index_t2=460;
index_t1=435;

T=t(index_t2)-t(index_t1); 
wr=2*pi/T;  
wn=wr/(sqrt(1-2*zeta^2)); 

H=tf(K*wn^2,[1 2*zeta*wn wn^2]);
A=[0 1; -wn^2 -2*zeta*wn];
B=[0; K*wn^2];
C=[1 0];
D=0;

Sys=ss(A,B,C,D);
[y1_cond]=lsim(Sys,u,t,[y1(1),(y1(2)-y1(1))/(t(2)-t(1))]);
figure; 
plot(t,[y1,y1_cond]);title('Compararea semnalului initial cu semnalul simulat')
legend('y1','y1 simulat')

J=1/sqrt(length(t))*norm(y1-y1_cond);            
eMPN=norm(y1-y1_cond)/norm(y1-mean(y1))*100;  
y1_initial=iddata(y1,u,T);
y1_simulat=iddata(y1_cond,u,T);
figure;
compare(y1_initial,y1_simulat);

%%
% Partea 2
% Estimarea Diagramei Bode pentru semnal y1


t=scope37(:,1);
u=scope37(:,2);
y1=scope37(:,3);
y2=scope37(:,4);

plot(t,u,t,y1);

%frecv joase
umin1=124;
ymin1=131;
umax1=64;
ymax1=67;

M1=(y1(ymax1)-y1(ymin1))/(u(umax1)-u(umin1));  
w1=pi/(t(ymin1)-t(ymax1)); 

delta_t1=t(umax1)-t(ymax1);
ph1=rad2deg(delta_t1*w1);  

umin2=289;
ymin2=294;
umax2=270;
ymax2=273;

M2=(y2(ymax2)-y1(ymin2))/(u(umax2)-u(umin2));
w2=pi/(t(ymin2)-t(ymax2));

delta_t2=t(umax2)-t(ymax2);
ph2=rad2deg(delta_t2*w2); 

umin3=326;
ymin3=330;
umax3=309;
ymax3=313;

M3=(y1(ymax3)-y1(ymin3))/(u(umax3)-u(umin3)); 
w3=pi/(t(ymin3)-t(ymax3)); 
delta_t3=t(umax3)-t(ymax3);
ph3=rad2deg(delta_t3*w3);

umin4=390;
ymin4=394;
umax4=376;
ymax4=379;

M4=(y1(ymax4)-y1(ymin4))/(u(umax4)-u(umin4)); 
w4=pi/(t(ymin4)-t(ymax4)); 
delta_t4=t(umin4)-t(ymin4); 
ph4=rad2deg(delta_t4*w4); 

%frecv medii

umin5=469;
ymin5=473;
umax5=456;
ymax5=460;

M5=(y1(ymax5)-y1(ymin5))/(u(umax5)-u(umin5));
w5=pi/(t(ymin5)-t(ymax5));
delta_t5=t(umax5)-t(ymax5);
ph5=rad2deg(delta_t5*w5); 

%rezonanta

t_u_min=492; 
t_u_max=480;
t_y1_min=496;
t_y1_max=485;

Mr1=(y1(t_y1_max)-y1(t_y1_min))/(u(t_u_max)-u(t_u_min)); 
wr1=2*pi/T ;   
delta_rez=t(t_u_min)-t(t_y1_min);  
phr1=rad2deg(delta_rez*wr1);  

umin6=556;
ymin6=561;
umax6=546;
ymax6=550;

M6=(y1(ymax6)-y1(ymin6))/(u(umax6)-u(umin6));
w6=pi/(t(ymin6)-t(ymax6)); 
delta_t6=t(umin6)-t(ymin6); 
ph6=rad2deg(delta_t6*w6); 

umin7=576;
ymin7=581;
umax7=566;
ymax7=571;

M7=(y1(ymax7)-y1(ymin7))/(u(umax7)-u(umin7));
w7=pi/(t(ymin7)-t(ymax7));
delta_t7=t(umax7)-t(ymax7);
ph7=rad2deg(delta_t7*w7);

%frecv inalte

umin8=667;
ymin8=671;
umax8=658;
ymax8=663;

M8=(y1(ymax8)-y1(ymin8))/(u(umax8)-u(umin8));
w8=pi/(t(ymin8)-t(ymax8));
delta_t8=t(umax8)-t(ymax8);
ph8=rad2deg(delta_t8*w8); 

umin9=979;
ymin9=983;
umax9=973;
ymax9=978;

M9=(y1(ymax9)-y1(ymin9))/(u(umax9)-u(umin9));
w9=pi/(t(ymin9)-t(ymax9));
delta_t9=t(umin9)-t(ymin9);
ph9=rad2deg(delta_t9*w9);

umin10=990;
ymin10=995;
umax10=984;
ymax10=989;

M10=(y1(ymax10)-y1(ymin10))/(u(umax10)-u(umin10));
w10=pi/(t(ymin10)-t(ymax10));

delta_t10=t(umin10)-t(ymin10);
ph10=rad2deg(delta_t10*w10);

umin11=817;
ymin11=822;
umax11=810;
ymax11=815;

M11=(y1(ymax11)-y1(ymin11))/(u(umax11)-u(umin11));
w11=pi/(t(ymin11)-t(ymax11));

delta_t11=t(umin11)-t(ymin11);
ph11=rad2deg(delta_t11*w11);

w=logspace(3,6);
[num,den]=tfdata(H,'v');
[M,ph]=bode(num,den,w);

vectorul_modulului=[M1,M2,M3,M4,M5,Mr1,M6,M7,M8,M9,M10,M11];
vectorul_fazei=[ph1,ph2,ph3,ph4,ph5,phr1,ph6,ph7,ph8,ph9,ph10,ph11];
vectorul_pulsatiei=[w1,w2,w3,w4,w5,wr1,w6,w7,w8,w9,w10,w11];

figure;
subplot(211);
semilogx(w,20*log10(M),vectorul_pulsatiei,20*log10(vectorul_modulului),"x");
grid on;ylabel('Amplitudine(dB)');title('Modul');hold on;

subplot(212);
semilogx(w, squeeze(ph),'color','green'); hold on; 

semilogx(vectorul_pulsatiei, vectorul_fazei, "x");
grid on;xlabel('Frecventa (rad/s)');ylabel('Faza(grade)');title('Faza');
 


Panta=20*log10(M11/M8)*10*(w11/w8); 

%%
%Partea 3
%Identificarea parametrica a semnalului y1

t=scope37(:,1);
u=scope37(:,2);
y1=scope37(:,3);
y2=scope37(:,4);

%semnalul y1
T=t(2)-t(1);
%%
dy1=iddata(y1,u,T);
M_arx_y1=arx(dy1,[2,2,0])

figure;
resid(dy1,M_arx_y1,'corr',5);
figure;
compare(dy1,M_arx_y1);

Hz=tf(M_arx_y1.B,M_arx_y1.A,T)
Hc=d2c(Hz)
%%
My1_pem=pem(dy1,M_arx_y1);
figure;
resid(dy1,My1_pem,7);
figure;
compare(dy1,My1_pem);

%%
M_armax_y1=armax(dy1,[2,2,2,0]);
figure;
resid(dy1,M_armax_y1);

figure;
compare(dy1,M_armax_y1);

Hz=tf(M_armax_y1.B,M_armax_y1.A,T)
Hc=d2c(Hz)

%%
M_iv4_y1=iv4(dy1,[2,2,0]);
figure;
resid(dy1,M_iv4_y1,5)

figure
compare(dy1,M_iv4_y1)
%%
M_oe_y1=oe(dy1,[2,2,0]);
figure;
resid(dy1,M_oe_y1,7)
figure;
compare(dy1,M_oe_y1)
Hz=tf(M_oe_y1.B,M_oe_y1.F,T)
Hc=d2c(Hz)

%%

dy1_n4=iddata(y1(50:end,1),u(50:end,1),T);
My1_n4sid=n4sid(dy1_n4,1:7); 

figure;
resid(dy1_n4,My1_n4sid,5);

figure;
compare(dy1_n4,My1_n4sid);

index_t2=460;
index_t1=435;

T=t(index_t2)-t(index_t1); 

y1c=dlsim(My1_n4sid.A,My1_n4sid.B,My1_n4sid.C(1,:),My1_n4sid.D(1,:),u);

[A,B,C,D]=d2cm(My1_n4sid.A,My1_n4sid.B,My1_n4sid.C(1,:),My1_n4sid.D(1,:),T);

plot(t(50:end,1),[y1(50:end,1),y1c(50:end,1)]);legend('y1','y1c');

[Ny1,Dy1]=ss2tf(A,B,C,D);
Hy1=tf(Ny1,Dy1);

wny1=sqrt(Dy1(1,3)); 
titay1=Dy1(1,2)/2/wny1;  
Ky1=Ny1(1,3)/Dy1(1,3); 

J=1/sqrt(length(t))*norm(y1-y1c);     
eMPN=norm(y1-y1c)/norm(y1-mean(y1))*100; 







 