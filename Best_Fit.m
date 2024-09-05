clc
clear all
close all
%% Wing Drag Divergence
    Cl = 0.54;
    % Conventional
    DeltaM_div = -0.3496*Cl + 0.1921;

%% t/c from M_div
    M_div = 0.8007;
    % Conventional
    TC = -0.4673*M_div + 0.4854;

%% CL_Max
    cc = .06638;
    % Takeoff
    Clmax_takeoff = -15.7160*exp(-4.2422*cc) + 16.7680*exp(-2.9683*cc);
    % Landing
    Clmax_landing = 3.3287*exp(0.3480*cc) + -1.4142*exp(-13.6025*cc);

%% Fuel Ratio/All Out Range
    R_ao = 5546;
    % JT8D
    WF_WT = 0.3804*exp(0*R_ao) - 0.3789*exp(-0.0003*R_ao);

    %% Skin Friction
    Cf = 230.7517*Rn^-0.2891 + 1.0836;