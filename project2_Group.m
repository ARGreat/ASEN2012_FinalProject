%{
Author: Alex Godbout
Assignment: ASEN2012 Project 2
Creation Date: 11/21/2024
Inputs: Initial Bottle Parameters (Begins on line 21)
Outputs: Trajectory and Thrust of Bottle Rocket
Purpose: Devlope a model that can be used to find the ideal starting
parameters for a bottle rocket to acheive maximum distance traveled.
%}


clear;
clc;
close all;
%--------------------------%
%---Load Data/Constants----%
%--------------------------%
%The maximum and Minimum parameters
        %C_d  ,V_water, p_0 ,angle
bounds = [.425,  0    ,  0  , 0    ; %Min
           1,0.002    , 75  ,90    ];%Max


%State vectorvector contains any value from 0-1 that represents
%A fraction of the bounds
statevector_0 = [0.5,0.5,0.5,.5]; %STarting Position
statevector = statevector_0;

%---------------%
%---Run ode45---%
%---------------%

history = []; % To track every simulation preventing recalculations later on
stepSize = 0.1; %How far to move in rectilinear space after each step

%take "12" steps in the direction of the gradient
%Decrease step size each time to increase precision
dist = zeros(12,1);
for i = 1:12
    [statevector(i+1,:),dist(i+1),history] = calculateGradient(statevector(i,:),stepSize,bounds,history);
    stepSize=stepSize*(0.75);
end


%------------------%
%---Process Data---%
%------------------%

plot3(history(:,5),history(:,3),history(:,1));

%Dial In Pressure for Distance

%Begin by rounding down volume of water and angle to simpler number
idealLaunchParameters = [statevector(end,1),round(statevector(end,2),2),statevector(end,3),round(statevector(end,4),1)];
dist;
delta = 1;
%Get Fairly Close
while(abs(delta) > .1)
    dist = simulateLaunch(idealLaunchParameters, bounds,history);
    delta = 92-dist
    idealLaunchParameters = idealLaunchParameters +[0,0,.001,0]*delta/abs(delta);
end
%Get VERRRY CLOSE
while(abs(delta) > .01)
    dist = simulateLaunch(idealLaunchParameters, bounds,history);
    delta = 92-dist
    idealLaunchParameters = idealLaunchParameters +[0,0,.0001,0]*delta/abs(delta);
end

%Convert fractions to actual units
idealLaunchParameters = bounds(1,:)+idealLaunchParameters.*(bounds(2,:)-bounds(1,:));

%-----------%
%---Plots---%
%-----------%

%Coefficient of drag Analysis
sizeAnalysis = 50
steps(1,:) = linspace(0,1,sizeAnalysis);

for i=1:sizeAnalysis
    tempVector = unitsToFraction(idealLaunchParameters,bounds);
    tempVector(1) = steps(i);
    [cd(i),~] = simulateLaunch(tempVector,bounds,history);
    tempVector = unitsToFraction(idealLaunchParameters,bounds);
    tempVector(2) = steps(i);
    [v_w(i),~] = simulateLaunch(tempVector,bounds,history);
    tempVector = unitsToFraction(idealLaunchParameters,bounds);
    tempVector(3) = steps(i);
    [p_0(i),~] = simulateLaunch(tempVector,bounds,history);
    tempVector = unitsToFraction(idealLaunchParameters,bounds);
    tempVector(4) = steps(i);
    [alpha(i),~] = simulateLaunch(tempVector,bounds,history);
end
figure; hold on;
plot(bounds(1,1)+steps*(bounds(2,1)-bounds(1,1)),cd);
title("Coefficient of Drag vs. Distance");
figure;
plot(steps*bounds(2,2),v_w);
title("Initial Volume of Water vs. Distance");
figure;
plot(steps*bounds(2,3),p_0);
title("Initial Pressure vs. Distance");
figure;
plot(steps*bounds(2,4),alpha);
title("Launch Angle vs. Distance");
function [dist, history] = simulateLaunch(vars,bounds,history)
    vars = bounds(1,:) + (bounds(2,:)-bounds(1,:)).*vars;
    [row, ~] =size(history);
    history(row+1,2:5) = vars;
    tspan = [0, 10]; %
    const = constFunc(vars);
    statevector_0 = initialFunc(const);
    [~,StateVector] =ode45(@(t, stateVector) d_State_dt(t,stateVector, const), tspan, statevector_0);
    dist = StateVector(end,1);
    history(row+1,1) = dist;
end

% Function to model the state change
function [d_statevector_dt, linVars] = StateFunction(t, stateVector, const)
    % Extract state variables
    x = stateVector(1);
    v_x = stateVector(2);
    z = stateVector(3);
    v_z = stateVector(4);
    mTotal = stateVector(5);
    Vol_air = stateVector(6);
    massAir = stateVector(7);
    

    F_Thrust = 0;
    F_Grav = [0 -mTotal*const.g];
    F_Normal = 0;
    d_Vair = 0;
    d_mAir = 0;
    d_mTot = 0;
    
    v_exit = 0;
    rho_exit = 0;
    M = 0;
    %---On the Stand---%
    if x < const.l_s*cosd(const.theta_0) 

            h = [cosd(const.theta_0) sind(const.theta_0)];
            h = h/norm(h);

            %Didnt Need this lmao, I think it should be included tho
            %Account for Normal Force of Launcher
            %magNorm = abs(F_Grav(2)/cosd(const.theta_0));
            %F_Normal = [-magNorm*sind(const.theta_0) magNorm*cosd(const.theta_0)];

    else %Off Launcher
            h = [v_x v_z];
            h = h/norm(h);
    end

    %Always calculate so p_end can be used in Stage 2 if statement
    p_end = const.p_0*(const.vAir_0/const.volumeBottle)^(const.gamma);

    %---Stage 1---%
    if Vol_air < const.volumeBottle
        p = (const.p_0)*((const.vAir_0)/Vol_air).^const.gamma;
        if(p > const.atomosphericPressure)
            %Calculate the Change in Volume
            %Latex form p_{\text{air}}' = C_d \cdot A_t \cdot \sqrt{\frac{2}{\rho_w} \left( p_0 \left( \frac{V_{w,0}}{V_{\text{air}}} \right)^\lambda - p_{\text{atm}} \right)}
            d_Vair = (const.c_dis)*(const.A_Throat)*sqrt(2/(const.rhoWater)*(p-const.atomosphericPressure));
            d_mTot = -d_Vair*const.rhoWater;   
            %Calculate Net Force
            F_Thrust = 2*(const.c_dis)*(const.A_Throat)*(p-const.atomosphericPressure);
        else
            F_Thrust = 0;
        end
    %---Stage 2---%
    elseif p_end*(massAir/const.mAir_0)^(const.gamma) > const.atomosphericPressure


        %Compute State
        T_end = const.T_0*(const.vAir_0/const.volumeBottle)^(const.gamma-1);
        p = p_end*(massAir/const.mAir_0)^(const.gamma);
        rho_air = massAir/const.volumeBottle;
        T = p/(rho_air*const.Rair);
        %Type of Flow
        p_star = p*(2/(const.gamma+1)).^(const.gamma/(const.gamma-1));
        %Choked
        if p_star > const.atomosphericPressure 
            p_exit = p_star;
            T_exit = T*(2/(const.gamma+1));
            rho_exit = p_exit/(const.Rair*T_exit);
            v_exit = sqrt(const.gamma*const.Rair*T_exit);
        %Unchoked
        else
            p_exit = const.atomosphericPressure;
            M = sqrt( 2/(const.gamma-1) * ((p/const.atomosphericPressure)^((const.gamma-1)/(const.gamma))-1));
            T_exit = T/(1+(const.gamma-1)/2*M^2);
            rho_exit = const.atomosphericPressure/(const.Rair*T_exit);
            v_exit = M*sqrt(const.gamma*const.Rair*T_exit);
        end

        d_mAir = -const.c_dis*rho_exit*const.A_Throat*v_exit;
        F_Thrust = (-d_mAir*v_exit + (p_exit-const.atomosphericPressure)*const.A_Throat);
        d_mTot = d_mAir;    
    end

    %---Stage 3---%
    %Because the only force in stage 3 is the drag force and drag force is
    %always implemtned we do not have a stage 3 if statement

    %Always account for drag
    F_Drag = 1/2*const.rhoAir*norm([v_x v_z])^2*const.c_D*const.A_Bottle;
    
    %Add a Normal force to 'model' impact force
    if z <=0
        F_Normal = -h *(v_x^2 +v_z^2)*1000;
    end

    %Calculate Net Force and Acceleration (h breaks thrust and drag into x and z components)
    F_Net = h*(F_Thrust - F_Drag) + F_Normal + F_Grav;
    acceleration = F_Net/mTotal;
    
    
    linVars = [F_Thrust, v_exit , rho_exit,M];
    d_statevector_dt = [ v_x ;  acceleration(1)  ; v_z ; acceleration(2)   ; d_mTot;d_Vair ; d_mAir];
end

function const = constFunc(vars)
const.atomosphericPressure = 83426.563088; %N/m^2

%Control Variables
const.c_D = vars(1); 
const.vWater_0 = vars(2);%M^3
const.p_0 = (vars(3)+12.1)*6894.76; %N/m
const.theta_0 = vars(4); %Degrees

%Posible Adjustable Variables
const.A_Throat =(.021)^2/4*pi; %m
const.A_Bottle = (.105)^2/4*pi; %m
const.T_0 = 310;% K

%Position Vars
const.x_0 = 0; %m
const.z_0 = .25;%m
const.v_0 = 0;% m/s

%Launcher Vars
const.volumeBottle = 0.002;% m^3
const.l_s = 0.5;% m
const.massBottle = 0.15;% kg

%Constants
const.gamma = 1.4;
const.atomosphericPressure = 83426.563088; %N/m^2
const.rhoAir = 0.961;% kg/m^3
const.g = 9.81;% m/s^2`
const.rhoWater = 1000;% kg/m^3
const.Rair = 287;% J/(kg*K)
const.c_dis = 0.78;

%Helpful For Calculations
const.vAir_0 = const.volumeBottle-const.vWater_0;
const.mAir_0 = (const.vAir_0*const.p_0)/(const.T_0*const.Rair);
const.mass_0 = const.massBottle + const.vWater_0*const.rhoWater + const.mAir_0;
end

function initial = initialFunc(const)
initial = [const.x_0 const.v_0 const.z_0 const.v_0 const.mass_0 const.vAir_0 const.mAir_0];
end

%Function to seperate linear variables and derivatives for ode45
function d_statevector_dt=d_State_dt(t, stateVector, const)
    [d_statevector_dt,~] = StateFunction(t, stateVector, const);
end

function [statevector_0, distFinal, history]= calculateGradient(statevector_0,stepSize,bounds,history);
    
    d_statevector_0 = diag(zeros(1,4)+stepSize);
    [dist_0,history] = simulateLaunch(statevector_0,bounds, history);
    gradient = zeros(1,4);
    for i=1:4
        [dist_delta,~] = simulateLaunch(verifyStateVector(statevector_0 + d_statevector_0(i,:)),bounds,history);
        gradient(i) = (dist_delta-dist_0);
    end

    gradient = gradient/norm(gradient);

    statevector_0 = verifyStateVector(statevector_0 + gradient.*stepSize);
    [dist_delta,history] = simulateLaunch(statevector_0,bounds,history);
    while(dist_delta > dist_0)
        dist_0 = dist_delta;
        statevector_0 = verifyStateVector(statevector_0 + gradient.*stepSize);
        [dist_delta,history] = simulateLaunch(statevector_0,bounds,history);
    end
    distFinal = dist_0;
end

function statevector = verifyStateVector(statevector)
    for i=1:4
        statevector(i) = max(statevector(i),0);
        statevector(i) = min(statevector(i),1);
    end
end

function statevector=unitsToFraction(statevector,bounds)
    statevector=(statevector-bounds(1,:))./(bounds(2,:)-bounds(1,:));
end
