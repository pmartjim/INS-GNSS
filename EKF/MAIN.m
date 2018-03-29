clear; format short; clc; close all;

% Parameters
epochsInit= 5000; % Number of epochs used to calcualte bias
R= 0.001^2 * eye(3); % Measurement variance matrix
SWITCH_KF_UPDATE= 0;
f_GPS = 2;
dT= 1/f_GPS; % KF Update period
f_IMU = 2000; % IMU data frequency (1/DT)

% Convert from sensor-frame to body-frame (NED)
R_SENSOR2BODY = zeros(6);
R_SENSOR2BODY(4:6,4:6)= [-1 0 0;
                          0 1 0;
                          0 0 -1];
R_SENSOR2BODY(1:3,1:3)= [-1 0 0;
                          0 1 0;
                          0 0 -1];

% Tau for bias -- from manufacturer
tau= 1000;  %Actual manuf. value
BiasInstabW= deg2rad(0.5/3600);     % rad/s
BiasInstabA= 0.05 * 9.80665 / 1000; % m/s2

% % Measured white noise SD
% sigma_IMU= [0.019342455000080;  % ax
%             0.018660194431183;  % ay
%             0.019979928602605;  % az
%             0.001239820655196;  % wx [rad]
%             0.001353946031205;  % wy [rad]
%             0.001297525572326]; % wx [rad]

% Manufacturer Data
sigma_a = 0.06/60/sqrt(1/f_IMU);
sigma_w = deg2rad(0.15/60/sqrt(1/f_IMU));
sigma_IMU= [sigma_a;  % ax
            sigma_a;  % ay
            sigma_a;  % az
            sigma_w;  % wx [rad]
            sigma_w;  % wy [rad]
            sigma_w]; % wx [rad]
       
% Convert to PSD
sigma_IMU= sigma_IMU * sqrt(1/f_IMU);

% PSD for continuous model
Qw= diag([sigma_IMU', ...
    BiasInstabA, BiasInstabA, BiasInstabA, ...
    BiasInstabW, BiasInstabW, BiasInstabW]).^2;


testInit= 1; testLast = 10;
for testNum= testInit:testLast
    disp(testNum);
        
    file= strcat('../DATA_STATIC/24x120k/20180228_12', num2str(testNum),'.txt');
    [timeIMU, gyrox, gyroy, gyroz, accx, accy, accz,...
        incx, incy, incz, gyroSts, accSts, ~, ~, ~, ~]= DataRead(file); % rads
    u= [accx, accy, accz, gyrox, gyroy, gyroz]';
    dt= diff(timeIMU);
    timeIMU_final= timeIMU(end);
    N= length(timeIMU);
    
    % Rotate to body frame (NED)
    u= R_SENSOR2BODY * u;
    accx= u(1,:);  accy= u(2,:);  accz= u(3,:);
    gyrox= u(4,:); gyroy= u(5,:); gyroz= u(6,:);
    I = R_SENSOR2BODY(1:3,1:3) * [incx, incy, incz]';
    incx = I(1,:);    incy = I(2,:);    incz = I(3,:);
    
    % Initial values/biases
    ax0= mean(accx(1:epochsInit));
    ay0= mean(accy(1:epochsInit));
    az0= mean(accz(1:epochsInit));
    
    % Inclinometer for initialization 
    ix0= mean(incx(1:epochsInit));
    iy0= mean(incy(1:epochsInit));
    iz0= mean(incz(1:epochsInit));
    phi0= (atan2(-iy0,-iz0));
    theta0= (atan2(ix0, sqrt(iy0^2 + iz0^2)));
    clear incx incy incz ix0 iy0 iz0 I 
    
    % G estimation
    g= sqrt(ax0^2 + ay0^2 + az0^2);
    G= [0; 0; -g];
    
    % Gyro bias estimation
    wx0= mean(gyrox(1:epochsInit));
    wy0= mean(gyroy(1:epochsInit));
    wz0= mean(gyroz(1:epochsInit));  
    
    % Initial biases
    u0= [ax0; ay0; az0; wx0; wy0; wz0];
    
    % Allocate variables
    Var= zeros(15, ceil(timeIMU_final/dT)); %For robustness
    Pk= zeros(15);
    x= zeros(15,N);
    xcorrected= zeros(15,N);
    x_star= [0 0 0 0 0 0 phi0 theta0 0 0 0 0 0 0 0]';    % Initialization state
    x_dot_star= zeros(15,1);
    x(:,1)= x_star;
    xcorrected(:,1)= x_star;
    t_sum= 0;
    k_update= 1;
    
    % Initial linearization & discretization
    [Fn, F, Gu_tilda, Gw, H]=...
        Matrices(x(:,1), u0(1), u0(2), u0(3), u0(4), u0(5), u0(6), tau);
    [Phi, Gamma, Gammaw_W_Gammaw]= Matrices2Discrete( F, Gu_tilda, Gw, H, Qw, dT );
    
    
    for k= 1:N-1
        
        % Simple integration
        x_dot= x_dot_star + F * ( x(:,k) - x_star ) + Gu_tilda * ( u(:,k) - u0 );  %General form of discrete integration
        x(:,k+1)= x(:,k) + dt(k)*x_dot;
        
        % KF update
        if t_sum >= dT
            
            % System Matrices
            [Fn, F, Gu_tilda, Gw, H]=...
                Matrices(x(:,k), u(1,k), u(2,k), u(3,k), u(4,k), u(5,k), u(6,k), tau);
            
            % New Star state
            x_star = x(:,k);
            x_dot_star = x_dot;
            u0 = [u(1,k), u(2,k), u(3,k), u(4,k), u(5,k), u(6,k)]';
            
            % Discretize system for GPS upate time (only for variance calculations)
            [Phi, Gamma, Gammaw_W_Gammaw]= Matrices2Discrete( F, Gu_tilda, Gw, H, Qw, dT );

            Var(:,k_update) = diag(Pk);
            Pk= Phi*Pk*Phi' + Gammaw_W_Gammaw;
            
            if SWITCH_KF_UPDATE
                disp('---- KF UPDATE ----')
                L= Pk*H'/(H*Pk*H' + R);
                ymeasured= [0 0 0]';
                y= H*x(:,k+1);
                dy= ymeasured - y;
                xcorrected(:,k+1)= x(:,k+1) + L*dy;
                Pk= Pk - L*H*Pk;
                x(:,k+1)= xcorrected(:,k+1);
            end
            
            t_sum= 0;
            k_update= k_update+1;
        end
        t_sum= t_sum + dt(k);
                
    end  
    Var(:, ceil(timeIMU_final/dT))= diag(Pk); %For robustness
    
    % Plot errors for this run
    figure(1); hold on;
    
    subplot(3,3,1); hold on;
    plot(timeIMU,x(1,:));
    ylabel('x [m]');
    
    subplot(3,3,2); hold on;
    plot(timeIMU,x(2,:));
    ylabel('y [m]');
    
    subplot(3,3,3); hold on;
    plot(timeIMU,x(3,:))
    ylabel('z [m]');
    
    subplot(3,3,4); hold on;
    plot(timeIMU,x(4,:))
    ylabel('v_x [m/s]');
    
    subplot(3,3,5); hold on;
    plot(timeIMU,x(5,:))
    ylabel('v_y [m/s]');
    
    subplot(3,3,6); hold on;
    plot(timeIMU,x(6,:));
    ylabel('v_z [m/s]');
    
    subplot(3,3,7); hold on;
    plot(timeIMU,rad2deg(x(7,:)))
    ylabel('\phi [deg]');
    
    subplot(3,3,8); hold on;
    plot(timeIMU,rad2deg(x(8,:)))
    ylabel('\theta [deg]');
    
    subplot(3,3,9); hold on;
    plot(timeIMU,rad2deg(x(9,:)))
    ylabel('\psi [deg]');
    
    
end

TT = 0:dT:timeIMU_final;
deviation = sqrt(Var);  %1sigma value

% Plots
subplot(3,3,1); hold on;
plot(TT, deviation(1,:), 'r.')
plot(TT, -deviation(1,:), 'r.')

subplot(3,3,2); hold on;
plot(TT, deviation(2,:), 'r.')
plot(TT, -deviation(2,:), 'r.')

subplot(3,3,3); hold on;
plot(TT, deviation(3,:), 'r.')
plot(TT, -deviation(3,:), 'r.')

subplot(3,3,4); hold on;
plot(TT, deviation(4,:), 'r.')
plot(TT, -deviation(4,:), 'r.')

subplot(3,3,5); hold on;
plot(TT, deviation(5,:), 'r.')
plot(TT, -deviation(5,:), 'r.')

subplot(3,3,6); hold on;
plot(TT, deviation(6,:), 'r.')
plot(TT, -deviation(6,:), 'r.')

subplot(3,3,7); hold on;
plot(TT,rad2deg(deviation(7,:)), 'r.')
plot(TT,-rad2deg(deviation(7,:)), 'r.')

subplot(3,3,8); hold on;
plot(TT,rad2deg(deviation(8,:)), 'r.')
plot(TT,-rad2deg(deviation(8,:)), 'r.')

subplot(3,3,9); hold on;
plot(TT,rad2deg(deviation(9,:)), 'r.')
plot(TT,-rad2deg(deviation(9,:)), 'r.')




