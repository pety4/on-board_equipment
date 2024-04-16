%% ��������� ������������� ���������� ��������� ��
function [X0]=INTSTATE(H0,V0,THETA0,BETA0,GAM0,PSI0,L0,Z0,W0,BAL)
    V=V0;
    ALP=BAL(2);
    THETA = THETA0;
    BET = BETA0;
    %% ���������� ���������� ���������� ���������
    Vx1 = V*cos(BET)*cos(ALP);      % �������� � ���
    Vy1 =-V*cos(BET)*sin(ALP); 
    Vz1 = V*sin(BET);           
    Wx1 = W0(1);                    % ������� �������� � ���
    Wy1 = W0(2);      
    Wz1 = W0(3);      
    GAM = GAM0;                     % ���� ����������
    PSI = PSI0;          
    TET = ALP+THETA;
    H= H0;                          % ���������� � ���
    L= L0;
    Z= Z0;   % ��������� ���������� ������� �� ����� ���
    FI=0.0;
    LAMD=0.0;
    d=0.0;                          % ���������� �� ����� ��������
    Vwx1=0.0;                       % �������� ����� � ���
    Vwy1=0.0;
    Vwz1=0.0;
    deltaT=BAL(1);                  % ������ ����������
    deltaV=BAL(3);
    deltaN=BAL(4);
    deltaE=BAL(5);
    HFL=15.0;                       % �������� ���������� ������������
    EPS=0.0;                        % ���������� �� �������� ����������
    X1=0.0;                         % ��������������
    %% ������ ���������� ��������� ��
    X0=[Vx1,Vy1,Vz1,Wx1,Wy1,Wz1,GAM,PSI,TET,L,H,Z,FI,LAMD,d,deltaT,deltaV,deltaN,deltaE,Vwx1,Vwy1,Vwz1,HFL,EPS,X1];    
end
    


