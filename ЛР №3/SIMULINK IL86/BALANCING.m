%% �������� ������� ������ ������������ ������
function [BAL,J]=BALANCING(HTR,VTR,THETATR,GAMTR,BETATR)    
    %% �������� ������� ��������
    f = @(BAL)FCT(BAL,HTR,VTR,THETATR,GAMTR,BETATR);
    
    %% ��������� �������� ��������� ���������� ������������
    BAL0=[0.0,0.0,0.0,0.0,0.0];
%     options = optimset('fminunc');
%     options = optimset(options,'Display','iter','PlotFcns',@optimplotfval);
%     [BAL,J] = fminunc(f,BAL0,options);
    %% �����������
   [BAL,J] = fminunc(f,BAL0);  
end
