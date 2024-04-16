clear all;
%% Шаг интегрирования
DT=1.E-3;   
INITLA;
%% Параметры режимов
TGLIDE=100.0;
%% Инициализация состояния самолета
% Требуемые параметры
HTR = 500.0;
VTR = 80.0;
THETATR = 0.0;
GAMTR = 0.0;
BETATR = 0.0;
% Начальное условие
H0 = HTR;
V0 = VTR;
THETA0 = THETATR;
BETA0 = BETATR;
GAM0 = GAMTR;
PSI0 = 0.0;
L0 = 0.0;
Z0 = 100.0;
W0 = [0 0 0];
% Балансировка 
[BAL,J]=BALANCING(HTR,VTR,THETATR,GAMTR,BETATR);
% Инициализация начального состояния ЛА
[X0]=INTSTATE(H0,V0,THETA0,BETA0,GAM0,PSI0,L0,Z0,W0,BAL);
%% Программное управления
% Для режима горизонтального полета
U0_HORIZ=[BAL(1);BAL(3);BAL(4);BAL(5)];
% Для режима захвата глиссады
HTR=15.0;
THETATR=-3.0*pi/180;
[BAL,J]=BALANCING(HTR,VTR,THETATR,GAMTR,BETATR);
U0_GLIDE=[BAL(1);BAL(3);BAL(4);BAL(5)];
HTR=500.0;
THETATR=0.0;
%% Инициализация ИВК
TET0=X0(9);
GAM0=X0(7);
PSI0=X0(8);
% СКО ошибок начальной выставки
SIGVYST=0.1*pi/180/3600.0;
% погрешность выставки в горизонт
DTET = (-1+2*rand())*SIGVYST;
DPH = (-1+2*rand())*SIGVYST;	
% погрешность выставки по азимуту
DPSI = (-1+2*rand())*SIGVYST;

%TETB=TET0+DTET;
%GAMB=GAM0+DPH;
%PSIB=PSI0+DPSI;

TETB=TET0;
GAMB=GAM0;
PSIB=PSI0;

XB0(1)=cos(PSIB)*cos(TETB);
XB0(2)=-cos(PSIB)*sin(TETB)*cos(GAMB)+sin(PSIB)*sin(GAMB);
XB0(3)=cos(PSIB)*sin(TETB)*sin(GAMB)+sin(PSIB)*cos(GAMB);

XB0(4)=sin(TETB);
XB0(5)=cos(TETB)*cos(GAMB);
XB0(6)=-cos(TETB)*sin(GAMB);

XB0(7)=-cos(TETB)*sin(PSIB);	
XB0(8)=cos(PSIB)*sin(GAMB)+sin(PSIB)*sin(TETB)*cos(GAMB);	
XB0(9)=cos(PSIB)*cos(GAMB)-sin(PSIB)*sin(TETB)*sin(GAMB);

CBN0=[XB0(1) XB0(2) XB0(3); XB0(4) XB0(5) XB0(6); XB0(7) XB0(8) XB0(9)];

Vx1=X0(1);  Vy1=X0(2);  Vz1=X0(3);

XB0(10)=Vx1*cos(PSI0)*cos(TET0)-Vy1*(cos(PSI0)*sin(TET0)*cos(GAM0)-sin(PSI0)*sin(GAM0))+Vz1*(sin(PSI0)*cos(GAM0)+cos(PSI0)*sin(TET0)*sin(GAM0));
XB0(11)=Vx1*sin(TET0)+Vy1*cos(TET0)*cos(GAM0)-Vz1*cos(TET0)*sin(GAM0);
XB0(12)=-Vx1*sin(PSI0)*cos(TET0)+Vy1*(cos(PSI0)*sin(GAM0)+sin(PSI0)*sin(TET0)*cos(GAM0))+Vz1*(cos(PSI0)*cos(GAM0)-sin(PSI0)*sin(TET0)*sin(GAM0));

XB0(13)=X0(10);
XB0(14)=X0(11);
XB0(15)=X0(12);

X01=[XB0(10);XB0(11);XB0(12);XB0(13);XB0(14);XB0(15)];
%% Характеристики ДУС
SIGW0=1.0*pi/180/ 3600.0;			% град/ч ->рад/с
SIGKMW=1.E-4;
SIGFIW=15*pi/180/3600.0;            % угл. сек. ->рад

% Смещения нуля
dw=[SIGW0*sign(randn());
    SIGW0*sign(randn());
    SIGW0*sign(randn())];
% Погрешности масш. коэфф.
kmw=[SIGKMW*sign(randn());
    SIGKMW*sign(randn());
    SIGKMW*sign(randn())];
% Неортогональности
fiw=[SIGFIW*sign(randn());
    SIGFIW*sign(randn());
    SIGFIW*sign(randn());
    SIGFIW*sign(randn());
    SIGFIW*sign(randn());
    SIGFIW*sign(randn())];
% Характеристики шумов измерения ДУС
PSD_GYRO=1.E-12;      % rad2s-1 (Paul D Groves)

%% Характеристики ДЛУ
SIGA0=0.001;			% м/с^2
SIGKMA=1.E-4;
SIGFIA=15*pi/180/3600.0;            % угл. сек. ->рад

% Смещения нуля
da0=[SIGA0*sign(randn());
     SIGA0*sign(randn());
     SIGA0*sign(randn())];
% Погрешности масш. коэфф.
kma=[SIGKMA*sign(randn());
     SIGKMA*sign(randn());
     SIGKMA*sign(randn())];
% Неортогональности
fia=[SIGFIA*sign(randn());
     SIGFIA*sign(randn());
     SIGFIA*sign(randn());
     SIGFIA*sign(randn());
     SIGFIA*sign(randn());
     SIGFIA*sign(randn())];
% Характеристики шумов измерения ДЛУ
PSD_ACEL=1.E-7;     % m2s-3(Paul D Groves)