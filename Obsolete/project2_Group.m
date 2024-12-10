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



absMin = [.425,  .0005,  40, 20];
absMax = [   1,.002, 68,45];
newMin = absMin;
newMax = absMax;
%---------------%
%---Run ode45---%
%---------------%

%%[t, stateVectors] = ode45(@(t, stateVector) d_State_dt(t, stateVector, const), tspan, statevector_0);
numloops=20;
initialConds =   zeros(numloops,4);
maxDist = zeros(numloops,1);
for i=1:numloops
    dist = zeros(9,1);
    testValues = GenTestValues(newMin,newMax);
    testMatrix = GenTestMatrix(testValues);
    for j=1:9
        const = constFunc(testMatrix(j,:));
        dist(j) = simulateLaunch(const);
        "complete test " + j + " in loop " + i
    end
    maxDist(i) = max(dist);
    initialConds(i,:) = newInitialVals(dist,testValues);
    [newMin,newMax] = genMinValues(initialConds(i,:),absMin,absMax,i);
end

%------------------%
%---Process Data---%
%------------------%

function dist = simulateLaunch(const)
    tspan = [0, 10]; %
    statevector_0 = initialFunc(const);
    [~,StateVector] =ode45(@(t, stateVector) d_State_dt(t, stateVector, const), tspan, statevector_0);
    dist = max(StateVector(:,1));
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

function testValues = GenTestValues(min, max)
    cd = [min(1), (min(1)+max(1))/2 , max(1)];
    w = [min(2), (min(2)+max(2))/2 , max(2)];
    p = [min(3), (min(3)+max(3))/2 , max(3)];
    a = [min(4), (min(4)+max(4))/2 , max(4)];
    testValues = [cd;w;p;a];
end
function testMatrix= GenTestMatrix(testValues)
%coefficient of drag, amount of water, rocket air pressure, and launch angle
cd = testValues(1,:);
w = testValues(2,:);
p = testValues(3,:);
a = testValues(4,:);

%Tanguchi Matrix P = 4 L = 3
testMatrix =    [ ...
                 cd(1),w(1), p(1),a(1) ; ... %Exp1
                 cd(1),w(2), p(2),a(2) ; ... %Exp2
                 cd(1),w(3), p(3),a(3) ; ... %Exp3
                 cd(2),w(1), p(2),a(3) ; ... %Exp4
                 cd(2),w(2), p(3),a(1) ; ... %Exp5
                 cd(2),w(3), p(1),a(3) ; ... %Exp6
                 cd(3),w(1), p(3),a(2) ; ... %Exp7
                 cd(3),w(2), p(1),a(3) ; ... %Exp8
                 cd(3),w(3), p(2),a(1) ; ... %Exp9
                 ];
end

function [Min, Max] = genMinValues(PrevResults, absMin, absMax,j)
    diff = absMax - absMin;
    for i = 1:4
        % Update min and max symmetrically around PrevResults
        Min(i) = PrevResults(i) - diff(i)/(2^(j+1));
        Max(i) = PrevResults(i) + diff(i)/(2^(j+1));

        % Enforce bounds
        Min(i) = max([Min(i), absMin(i)]);
        Max(i) = min([Max(i), absMax(i)]);
    end
end

function varIndexMatrix=newInitialVals(dist,testValues)

    a_bar(1) = dist(1) + dist(2) + dist(3);
    a_bar(2) = dist(4) + dist(5) + dist(6);
    a_bar(3) = dist(7) + dist(8) + dist(9);

    b_bar(1) = dist(1) + dist(4) + dist(7);
    b_bar(2) = dist(2) + dist(5) + dist(8);
    b_bar(3) = dist(3) + dist(6) + dist(9);

    c_bar(1) = dist(1) + dist(6) + dist(8);
    c_bar(2) = dist(2) + dist(4) + dist(9);
    c_bar(3) = dist(3) + dist(5) + dist(7);

    d_bar(1) = dist(1) + dist(5) + dist(9);
    d_bar(2) = dist(2) + dist(6) + dist(7);  
    d_bar(3) = dist(3) + dist(4) + dist(8);

    [~, varIndexMatrix(1)] = max(a_bar);
    [~, varIndexMatrix(2)] = max(b_bar);
    [~, varIndexMatrix(3)] = max(c_bar);
    [~, varIndexMatrix(4)] = max(d_bar);

    varIndexMatrix(1) = testValues(1,varIndexMatrix(1));
    varIndexMatrix(2) = testValues(2,varIndexMatrix(2));
    varIndexMatrix(3) = testValues(3,varIndexMatrix(3));
    varIndexMatrix(4) = testValues(4,varIndexMatrix(4));
end