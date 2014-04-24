% This script uses a linear strength vortex panel method to find the lift 
% and pressure coefficients on a two element wing. The points are taken 
% from 2008 formula SAE rear wing. This file generates the points for the 
% formula sae wing using the banana airfoil coordinates for both the main 
% element and the flap. 
clc 
clear all 
close all 
% points file goes naturally from leading edge clockwise. We want 
% clockwise, lower trailing edge. 
alpha = 2; %alpha is defined as angle of attack from main element chord line. 
AL = alpha / 57.2958; % get into radians 
ref_length = 1.33; % distance from leading edge of main element to trailing 
 % edge of flap 
 
 % points file, unit length 
x_flap = [1.3297217118151354, 1.2904572302954775, 1.1876614831263712, 1.0605989518274743, 0.95780320465836799, 0.9185387231387101, 0.95780320465836799, 1.0605989518274743, 1.1876614831263712, 1.2904572302954775,1.3297217118151354]
y_flap = [-0.34607115247081488, -0.31009396750223184, -0.22074564599458646, -0.12181764386367, -0.062350878985984071, -0.064121696391919336, -0.11754846462230117, -0.1941909870370957, -0.27092319524055669, -0.32720193772202183,-0.34607115247081488]
x_main = [1.0, 0.90450849718747373, 0.65450849718747373, 0.34549150281252633, 0.09549150281252633, 0.0, 0.095491502812526274, 0.34549150281252622, 0.65450849718747361, 0.90450849718747373,1.0]
y_main = [0.0, 0.013914300349052097, 0.040917400250269027, 0.059574699935282791, 0.046048900461581513, 0.0, -0.046048900461581499, -0.059574699935282797, -0.040917400250269034, -0.013914300349052097,0.0]
x_flap=x_flap';
y_flap=y_flap';
x_main=x_main';
y_main=y_main';
% Establish panel endpoints for each separate element 
 for i = 1:(length(x_flap)-1) 
 flap_PT1(i,1) = x_flap(i,1); 
 flap_PT2(i,1) = x_flap(i+1,1); 
 flap_PT1(i,2) = y_flap(i,1); 
 flap_PT2(i,2) = y_flap(i+1,1); 
end 
 
for i = 1:(length(x_main)-1) 
 main_PT1(i,1) = x_main(i,1); 
 main_PT2(i,1) = x_main(i+1,1); 
 main_PT1(i,2) = y_main(i,1); 
 main_PT2(i,2) = y_main(i+1,1); 
end 
 
M_flap =length(x_flap) - 1; % number of flap panels 
M_main = length(x_main)-1; % number of main element panels 
 
M = M_flap+ M_main; %total number of panels 
N = M+2; 
 
% Find panel angles theta (panel orientation angle) 
for i = 1:M_flap 
 DY_flap = flap_PT2(i,2) - flap_PT1(i,2); 
 DX_flap = flap_PT2(i,1) - flap_PT1(i,1); 
 TH_flap(i) = atan2(DY_flap,DX_flap); 
 DL_flap(i) = sqrt(DX_flap^2 + DY_flap^2); 
end 
 
for i = 1:M_main 
 DY_main = main_PT2(i,2) - main_PT1(i,2); 
 DX_main = main_PT2(i,1) - main_PT1(i,1); 
 TH_main(i) = atan2(DY_main,DX_main); 
 DL_main(i) = sqrt(DX_main^2 + DY_main^2); 
end 
 
% Establish Collocation Points 
for i = 1:M_flap 
 CO_flap(i,1) = (flap_PT2(i,1)-flap_PT1(i,1))/2 + flap_PT1(i,1); 
 CO_flap(i,2) = (flap_PT2(i,2)-flap_PT1(i,2))/2 + flap_PT1(i,2); 
end 
 
for i = 1:M_main 
 CO_main(i,1) = (main_PT2(i,1)-main_PT1(i,1))/2 + main_PT1(i,1); 
 CO_main(i,2) = (main_PT2(i,2)-main_PT1(i,2))/2 + main_PT1(i,2); 
end 
 
% merge all into one set 
% do flap then main 
CO = [CO_flap;CO_main]; 
TH = [TH_flap';TH_main']; 
DL = [DL_flap';DL_main']; 
 
for i = 1:M 
 % determine if we're on flap or main 
 
 if i<=M_flap 
 %we're dealing with a collocation point on the flap 
 for j = 1:M+1 
 if j <= M_flap 
 %we're dealing with both collocation points on flap 
 
 % we find the influence coefficient for a specific collocation 
 % point, i, on each of the pannels, j. Then we move on to the next 
 % collocation point. 
 
 % Convert Collocation Point To Local Panel Coords 
 XT = CO(i,1) - flap_PT1(j,1); 
 YT = CO(i,2) - flap_PT1(j,2); 
 X2T = flap_PT2(j,1) - flap_PT1(j,1); 
 Y2T = flap_PT2(j,2) - flap_PT1(j,2); 
 
 X = XT*cos(TH(j)) + YT*sin(TH(j)); % collocation point 
 Y = -XT*sin(TH(j)) + YT*cos(TH(j)); %collocation point 
 X1 = 0; 
 Y1 = 0; 
 X2 = X2T*cos(TH(j)) + Y2T*sin(TH(j)); 
 Y2 = 0; 
 
 % Find the length of r1,r2,theta1 and theta2 
 R1 = sqrt((X-X1)^2 + (Y-Y1)^2); %length from panel point 1 to collocation point, in panel coords
 R2 = sqrt((X-X2)^2 + (Y-Y2)^2); %length from panel point 2 to collocation point, in panel coords 
 TH1 = atan2(Y-Y1,X-X1); 
 TH2 = atan2(Y-Y2,X-X2); 
 
 if i == j 
 Y = 0; 
 TH1 = 0; 
 end 
 % Compute velocity components as functions of Gamma1 and gamma2. 
 % Velocity of panel j due to collocation point i 
 if i == j 
 U1L = -0.5*(X-X2)/(X2); 
 U2L = 0.5*(X)/(X2); 
 W1L = -0.15916; 
 W2L = 0.15916; 
 else 
 U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2); 
 U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2); 
 W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2); 
 W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2); 
 end 
 
 % Transform the local velocities into global velocity functions 
 U1 = U1L*cos(-TH(j)) + W1L*sin(-TH(j)); 
 U2 = U2L*cos(-TH(j)) + W2L*sin(-TH(j)); 
 W1 = -U1L*sin(-TH(j)) + W1L*cos(-TH(j)); 
 W2 = -U2L*sin(-TH(j)) + W2L*cos(-TH(j)); 
 
 % Compute the coefficients of gamma in the influence matrix. 
 if j == (1 | (M_flap+1)) 
 A(i,1) = -U1*sin(TH(i)) + W1*cos(TH(i)); 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,1) = U1*cos(TH(i)) + W1*sin(TH(i)); 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 elseif j== (M_flap) 
 A(i,M_flap) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 A(i,M_flap+1) = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M_flap) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 B(i,M_flap+1) = U2*cos(TH(i)) + W2*sin(TH(i)); 
 else 
 A(i,j) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,j) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 end 
 
 elseif j>(M_flap+1) 
 % collocation point i is on flap, collocation point j is on 
 % main 
 
 % we find the influence coefficient for a specific collocation 
 % point, i, on each of the pannels, j. Then we move on to the next 
 % collocation point. 
 
 % since there are two separate airfoils, point j=M_flap+1 
 % should not be solved for (this is the 2nd edge of the last panel). 
 % During the for loop, any j after 
 % this should be j-1, to reference the correct pannel, except when referring to matrix A. 
 
 % Convert Collocation Point To Local Panel Coords 
 XT = CO(i,1) - main_PT1(j-M_flap-1,1); 
 YT = CO(i,2) - main_PT1(j-M_flap-1,2); 
 X2T = main_PT2(j-M_flap-1,1) - main_PT1(j-M_flap-1,1); 
 Y2T = main_PT2(j-M_flap-1,2) - main_PT1(j-M_flap-1,2); 
 
 X = XT*cos(TH(j-1)) + YT*sin(TH(j-1)); % collocation point 
 Y = -XT*sin(TH(j-1)) + YT*cos(TH(j-1)); %collocation point 
 X1 = 0; 
 Y1 = 0; 
 X2 = X2T*cos(TH(j-1)) + Y2T*sin(TH(j-1)); 
 Y2 = 0; 
 
 % Find the length of r1,r2,theta1 and theta2 
 R1 = sqrt((X-X1)^2 + (Y-Y1)^2); %length from panel point 1 to collocation point, in panel coords 
 R2 = sqrt((X-X2)^2 + (Y-Y2)^2); %length from panel point 2 to collocation point, in panel coords 
 TH1 = atan2(Y-Y1,X-X1); 
 TH2 = atan2(Y-Y2,X-X2); 
 
 if i == j-1 
 Y = 0; 
 TH1 = 0; 
 end 
 % Compute velocity components as functions of Gamma1 and gamma2. 
 % Velocity of panel j due to collocation point i 
 if i == j-1 
 U1L = -0.5*(X-X2)/(X2); 
 U2L = 0.5*(X)/(X2); 
 W1L = -0.15916; 
 W2L = 0.15916; 
 else 
 U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2); 
 U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2); 
 W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2); 
 W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2); 
 end 
 
% Transform the local velocities into global velocity functions 
 U1 = U1L*cos(-TH(j-1)) + W1L*sin(-TH(j-1)); 
 U2 = U2L*cos(-TH(j-1)) + W2L*sin(-TH(j-1)); 
 W1 = -U1L*sin(-TH(j-1)) + W1L*cos(-TH(j-1)); 
 W2 = -U2L*sin(-TH(j-1)) + W2L*cos(-TH(j-1)); 
 
 % Compute the coefficients of gamma in the influence matrix. 
 if j == ( (M_flap+2)) 
 A(i,M_flap+2) = -U1*sin(TH(i)) + W1*cos(TH(i)); 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M_flap+2) = U1*cos(TH(i)) + W1*sin(TH(i)); 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 elseif j== ( M+1) 
 A(i,M+1) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 A(i,N) = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M+1) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 B(i,N) = U2*cos(TH(i)) + W2*sin(TH(i)); 
 else 
 A(i,j) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,j) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 end 
 else 
     end 
 end 
 A(i,N+1) = cos(AL)*sin(TH(i))-sin(AL)*cos(TH(i)); 
 
 else 
 %we're dealing with a collocation point on the main 
 for j = 1:M+1 
 if j <= M_flap 
 % CO point i is on main, CO point j is on flap 
 
 % we find the influence coefficient for a specific collocation 
 % point, i, on each of the pannels, j. Then we move on to the next 
 % collocation point. 
 
 % Convert Collocation Point To Local Panel Coords 
 XT = CO(i,1) - flap_PT1(j,1); 
 YT = CO(i,2) - flap_PT1(j,2); 
 X2T = flap_PT2(j,1) - flap_PT1(j,1); 
 Y2T = flap_PT2(j,2) - flap_PT1(j,2); 
 
 X = XT*cos(TH(j)) + YT*sin(TH(j)); % collocation point 
 Y = -XT*sin(TH(j)) + YT*cos(TH(j)); %collocation point 
 X1 = 0; 
 Y1 = 0; 
 X2 = X2T*cos(TH(j)) + Y2T*sin(TH(j)); 
 Y2 = 0; 
 
 % Find the length of r1,r2,theta1 and theta2 
 R1 = sqrt((X-X1)^2 + (Y-Y1)^2); %length from panel point 1 to collocation point, in panel coords 
 R2 = sqrt((X-X2)^2 + (Y-Y2)^2); %length from panel point 2 to collocation point, in panel coords 
 TH1 = atan2(Y-Y1,X-X1); 
 TH2 = atan2(Y-Y2,X-X2); 
 
 if i == j 
 Y = 0; 
 TH1 = 0; 
 end 
% Compute velocity components as functions of Gamma1 and gamma2. 
 % Velocity of panel j due to collocation point i 
 if i == j 
 U1L = -0.5*(X-X2)/(X2); 
 U2L = 0.5*(X)/(X2); 
 W1L = -0.15916; 
 W2L = 0.15916; 
 else 
 U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2); 
 U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2); 
 W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2); 
 W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2); 
 end 
 
 % Transform the local velocities into global velocity functions 
 U1 = U1L*cos(-TH(j)) + W1L*sin(-TH(j)); 
 U2 = U2L*cos(-TH(j)) + W2L*sin(-TH(j)); 
 W1 = -U1L*sin(-TH(j)) + W1L*cos(-TH(j)); 
 W2 = -U2L*sin(-TH(j)) + W2L*cos(-TH(j)); 
 
 % Compute the coefficients of gamma in the influence matrix. 
 if j == (1) 
 A(i,1) = -U1*sin(TH(i)) + W1*cos(TH(i)); 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,1) = U1*cos(TH(i)) + W1*sin(TH(i)); 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 elseif j== (M_flap) 
 A(i,M_flap) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 A(i,M_flap+1) = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M_flap) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 B(i,M_flap+1) = U2*cos(TH(i)) + W2*sin(TH(i)); 
 else 
 A(i,j) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,j) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 end 
 elseif j>(M_flap+1) 
 % collocation point i is on main, collocation point j is on 
 % main 
 
 % we find the influence coefficient for a specific collocation 
 % point, i, on each of the pannels, j. Then we move on to the next 
 % collocation point. 
 
 % since there are two separate airfoils, point j=M_flap+1 
 % should not be solved for (this is the 2nd edge of the last panel). 
 % During the for loop, any j after 
 % this should be j-1, to reference the correct pannel, except when referring to matrix A. 
 
 % Convert Collocation Point To Local Panel Coords 
 XT = CO(i,1) - main_PT1(j-M_flap-1,1); 
 YT = CO(i,2) - main_PT1(j-M_flap-1,2); 
 X2T = main_PT2(j-M_flap-1,1) - main_PT1(j-M_flap-1,1); 
 Y2T = main_PT2(j-M_flap-1,2) - main_PT1(j-M_flap-1,2); 
 
 X = XT*cos(TH(j-1)) + YT*sin(TH(j-1)); % collocation point 
 Y = -XT*sin(TH(j-1)) + YT*cos(TH(j-1)); %collocation point 
 X1 = 0; 
 Y1 = 0; 
 X2 = X2T*cos(TH(j-1)) + Y2T*sin(TH(j-1)); 
 Y2 = 0; 
 
 % Find the length of r1,r2,theta1 and theta2 
 R1 = sqrt((X-X1)^2 + (Y-Y1)^2); %length from panel point 1 to collocation point, in panel coords 
 R2 = sqrt((X-X2)^2 + (Y-Y2)^2); %length from panel point 2 to collocation point, in panel coords 
 TH1 = atan2(Y-Y1,X-X1); 
 TH2 = atan2(Y-Y2,X-X2); 
 
 if i == j-1 
 Y = 0; 
TH1 = 0; 
 end 
 % Compute velocity components as functions of Gamma1 and gamma2. 
 % Velocity of panel j due to collocation point i 
 if i == j-1 
 U1L = -0.5*(X-X2)/(X2); 
 U2L = 0.5*(X)/(X2); 
 W1L = -0.15916; 
 W2L = 0.15916; 
 else 
 U1L = -(Y*log(R2/R1)+X*(TH2-TH1)-X2*(TH2-TH1))/(6.28319*X2); 
 U2L = (Y*log(R2/R1) + X*(TH2-TH1))/(6.28319*X2); 
 W1L = -((X2-Y*(TH2-TH1)) - X*log(R1/R2) + X2*log(R1/R2))/(6.28319*X2); 
 W2L = ((X2 - Y*(TH2-TH1))-X*log(R1/R2))/(6.28319*X2); 
 end 
 
 % Transform the local velocities into global velocity functions 
 U1 = U1L*cos(-TH(j-1)) + W1L*sin(-TH(j-1)); 
 U2 = U2L*cos(-TH(j-1)) + W2L*sin(-TH(j-1)); 
 W1 = -U1L*sin(-TH(j-1)) + W1L*cos(-TH(j-1)); 
 W2 = -U2L*sin(-TH(j-1)) + W2L*cos(-TH(j-1)); 
 
 % Compute the coefficients of gamma in the influence matrix. 
 if j == ( (M_flap+2)) 
 A(i,M_flap+2) = -U1*sin(TH(i)) + W1*cos(TH(i)); 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M_flap+2) = U1*cos(TH(i)) + W1*sin(TH(i)); 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 elseif j== M+1 
 A(i,M+1) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 A(i,N) = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,M+1) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 B(i,N) = U2*cos(TH(i)) + W2*sin(TH(i)); 
 else 
 A(i,j) = -U1*sin(TH(i)) + W1*cos(TH(i)) + HOLDA; 
 HOLDA = -U2*sin(TH(i)) + W2*cos(TH(i)); 
 B(i,j) = U1*cos(TH(i)) + W1*sin(TH(i)) + HOLDB; 
 HOLDB = U2*cos(TH(i)) + W2*sin(TH(i)); 
 end 
 else 
 end 
 end 
 A(i,N+1) = cos(AL)*sin(TH(i))-sin(AL)*cos(TH(i)); 
 
 end 
end 
 
 % Add both kutta conditions. Be careful of where the ones are. 
 % matrix columns M_flap+1 and M+2 are the last edges of the airfoil. 
 
%flap 
A(N-1,1) = 1; 
A(N-1,M_flap+1) = 1; 
%main 
A(N,M_flap+2) = 1; 
A(N,N) = 1; 
 
R = rref(A); % solve 
G = R(:,N+1); 
CL = 0; 
% calculate variables of interest 
for i = 1:M 
 VEL = 0; 
 for j = 1:N 
 VEL = VEL + B(i,j)*G(j); 
 end 
 V = VEL + cos(AL)*cos(TH(i)) + sin(AL)*sin(TH(i)); 
 CP(i) = 1-V^2; 
 CL = CL + -1.*CP(i)*(cos(AL)*cos(TH(i)) + sin(AL)*sin(TH(i)))*DL(i); 
end 
CP = CP'; 
disp(CP);
% report values of interest 
CL = CL/ref_length % varies based on ref. length definition
unitxm = x_main(:,1)./ref_length; 
unitxf = x_flap(:,1)./ref_length; 
unitym = y_main(:,1)./ref_length; 
unityf = y_flap(:,1)./ref_length; 
figure 
plot(unitxm, unitym) 
title('Two Element Airfoil') 
hold on 
plot(unitxf,unityf) 
xlabel('x/c') 
ylabel('y/c') 
axis([0 1 -.3 .3]) 
 
figure 
plot(CO(:,1)./ref_length,CP,'o') 
title('Pressure Distribution on a Two Element Airfoil') 
xlabel('x/c') 
ylabel('pressure coefficient')