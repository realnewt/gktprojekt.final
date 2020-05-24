%% Konstanter
clear
clc
clf
CP=[1.39 0.3847 -1.846e-04 2.895e-08;
    16.05 0.2804 -1.091e-04 9.098e-09;
    27.14 0.009274 -1.3813e-05 7.645e-09; 
    32.24 0.001924 1.055e-05 -3.596e-09]; %matris med alla CP konstanter [J/mol/K]

CEPCI_Year_B=532.9;     %From kurs PM
CEPCI_Year_A=607.5;     %Average for 2019
E=1;     %svets verkningsgrad (= 1)
rhobed=1120;        %[kg/m3]

M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

%% Reaktor 1
options = optimset('Display', 'off');
e=0;
Wtot=1950;      %[kg] mass cat
XA_start=0;     %Omsättning
T0=950;     %[K]
FA0=78;    %[mol/s] inflöden
FB0=29.10;
FH0=11.66;
FW0=10*FA0;

HR=116.3e3;     %[J/mol] reaktionsentalpi
P=1;        %[bar]
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0,FH0),[0 Wtot],[XA_start T0]);

XA=Y(:,1); T=Y(:,2);

figure(1)
plot(W,XA)
xlabel('Mängd katalys(kg)','Fontsize',15)
ylabel('XA','Fontsize',15)
title('Omsättning mot mängd katalys, reaktor 1','Fontsize',15)
legend('XA')

figure(2)
plot(T,XA)
xlabel('Temp(K)','Fontsize',15)
ylabel('XA','Fontsize',15)
title('Omsättning mot temperaturfall, reaktor 1','Fontsize',15)
legend('Adiabat linje ')

%kostnad
Volym=Wtot/rhobed;          %[m3]
D=(2*Volym/pi)^(1/3);       %[m] diameter av reaktorkärl
L=2*D;      %[m] längd av reaktorkärl

Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)
t=(P*D)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 
V_inner=pi*(((D/2)-t)^2)*(L-2*t);     %[m^3] Volume of inner tank (air)
Volume_container=Volym-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_h_304=12800;      %konstanter för kostnadsberäkning
b_h_304=73;
n_h_304=0.85;
Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304; 
Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B); %[$] Accounts for inflation[SEK]
Cost_Dist_h_304=Cost_Year_A_h_304*10;        %[SEK]
Cost_Dist_h_304_with_kat_1=Cost_Dist_h_304*1.5;       %Cost with catalyst 
fprintf('REACTOR ONE\n')
fprintf('Inflow: %2.0f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\n',FA0,FB0,FH0,FW0)
fprintf('Dimensions: diameter is %0.3f m and length is %1.2f m. Volume is %0.2f m^3\n',D,L,Volym)
fprintf('Wall thickness is %0.3f mm and total mass of container is %2.0f kg. Amount of catalyst is %4.0f kg\n',t*1000,Shell_mass_dist,Wtot)
fprintf('Total cost is %6.0f SEK with a final conversion of %0.2f\n',Cost_Dist_h_304_with_kat_1,XA(end));
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))

%% Ugn 2
Tin=T(end); Tut=950; Tmedel=(Tin+Tut)/2;        %[K] definierar in och ut temperaturer

F_mass_ugn_1=[F_mas(FA0*(1-XA(end)),1);
              F_mas(FB0+XA(end)*FA0,2);
              F_mas(FH0+XA(end)*FA0,3);
              F_mas(FW0,4)];                   %[kg/s]

    CP_matrix_ugn_1=[Cp_new(Tmedel,1,1,1,1,1);
                     Cp_new(Tmedel,2,2,2,2,2);
                     Cp_new(Tmedel,3,3,3,3,3);
                     Cp_new(Tmedel,4,4,4,4,4)];     %[J/Kg*K] Cp each component

QBA=F_mass_ugn_1(1)*CP_matrix_ugn_1(1,1)*(Tut-Tin);
QBE=F_mass_ugn_1(2)*CP_matrix_ugn_1(2,1)*(Tut-Tin);
QH=F_mass_ugn_1(3)*CP_matrix_ugn_1(3,1)*(Tut-Tin);
QW=F_mass_ugn_1(4)*CP_matrix_ugn_1(4,1)*(Tut-Tin);      %[J/s]
Qtot=QBA+QBE+QH+QW;
         
chi=0.8;        %Verkningsgrad
Qheat=Qtot/chi;     %[J/s] Energi som krävs för uppvärmning 
Driftko_ugn_2=Qheat/1000*8000*0.2; %[SEK] Driftskostnad för 1 år, naturgas

%Cost Ugn: Cylindrical
a=80000;
b=109000;   %Kostnadsparametrar
n=0.8;
S_ugn_1=Qheat/(10^6);       %[kW]
Cost_Year_B=a+b*S_ugn_1.^n;
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_ugn_2=Cost_Year_A*10;

fprintf('OVEN ONE\n')
fprintf('Inflow: %2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))
fprintf('Effect needed is %4.0f kW and the operation cost per year is %8.0f SEK when heated with natural gas. Cost of investment is %7.0f SEK \n', Qheat/1000,Driftko_ugn_2,Cost_ugn_2)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour at a temperature of %3.0f K \n\n\n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,Tut)

%% Reaktor 2
fprintf('REAKTOR 2\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K\n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,Tut)
XA_start=XA(end);
Wtot=1950;
P=1;
[W,Y]=ode15s(@(W,Y)ode_func(W,Y,HR,P,CP,FA0,FB0,FW0,FH0),[0 Wtot],[XA_start T0]);
XA=Y(:,1); T=Y(:,2);

figure(3)
plot(W,XA)
xlabel('Mängd katalys(kg)','Fontsize',15)
ylabel('XA','Fontsize',15)
title('Omsättning mot mängd katalys, reaktor 2','Fontsize',15)

figure(4)
plot(T,XA)
xlabel('Temp(K)','Fontsize',15)
ylabel('XA','Fontsize',15)
title('Omsättning mot temperaturfall, reaktor 2','Fontsize',15)
legend('Adiabat linje ')
rhobed=1120; %kg/m3
Volym=Wtot/rhobed;
D=(2*Volym/pi)^(1/3);
L=2*D;

Density=8000;       %[kg/m^3] from kurs PM at 900 F=755 K (would be better with 950 K, yes?)
S_max=103.4*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)
t=(P*D)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 
V_inner=pi*(((D/2)-t)^2)*(L-2*t);     %[m^3] Volume of inner tank (air)
Volume_container=Volym-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

%Kostand horistontell
Cost_Year_B_h_304=a_h_304+b_h_304*Shell_mass_dist.^n_h_304;
Cost_Year_A_h_304=Cost_Year_B_h_304*(CEPCI_Year_A/CEPCI_Year_B); %[$] Accounts for inflation %[SEK]
Cost_Dist_h_304=Cost_Year_A_h_304*10;        %[SEK]
Cost_Dist_h_304_with_kat_2=Cost_Dist_h_304*1.5;       %Includes cost of catalyst

fprintf('Dimensions: diameter is %0.3f m and length is %1.2f. Volume is %0.2f m^3\n',D,L,Volym)
fprintf('Wall thickness is %0.3f mm and total mass of container is %2.0f kg. Amount of catalyst is %4.0f kg\n',t*1000,Shell_mass_dist,Wtot)
fprintf('Total cost is %f SEK with a final conversion of %0.2f\n',Cost_Dist_h_304_with_kat_2,XA(end));
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour at a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))

%% Cooler 1
format shortG
fprintf('COOLER 1\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))
%Data
%Area calculations
U_cooler_1=200;       %[W/m2] KVärmegenomgångstal from kurs PM
F_cooler_1=1;

%Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler_1=[(FA0*(1-XA(end)));
       (FB0+FA0*XA(end));
       (FH0+FA0*XA(end));
       (FW0)];       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass_cooler_1=[F_mas(F_mol_cooler_1(1,1),1);
        F_mas(F_mol_cooler_1(2,1),2);
        F_mas(F_mol_cooler_1(3,1),3);
        F_mas(F_mol_cooler_1(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O

% Setting the desired out temperatures and the current in temperatures
TH_in_cooler_1=T(end);      %[K] Initial temperature of mixture
TH_ut_cooler_1=319; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_cooler_1=273-20;       %Inlet temperature of coolant
TC_ut_cooler_1=273+20;        %Outlet temperature of coolant   

TH_medel_cooler_1=(TH_in_cooler_1+TH_ut_cooler_1)/2;   %an average to calculate Cp constants

CP_matrix_cooler_1=[Cp_new(TH_medel_cooler_1,1,1,1,1,1);
           Cp_new(TH_medel_cooler_1,2,2,2,2,2);
           Cp_new(TH_medel_cooler_1,3,3,3,3,3);
           Cp_new(TH_medel_cooler_1,4,4,4,4,4)];     %[J/Kg*K] Cp each component
    
       %Energy transfer calculations
q_matrix_cooler_1=F_mass_cooler_1(:,1).*CP_matrix_cooler_1(:,1).*(TH_in_cooler_1-TH_ut_cooler_1);      %Energy transfer matrix       
q_cooler_1=sum(q_matrix_cooler_1);      %[J/s]

delta_T_a_cooler_1=TH_in_cooler_1-TC_ut_cooler_1;
delta_T_b_cooler_1=TH_ut_cooler_1-TC_in_cooler_1;
delta_T_lm_cooler_1=(delta_T_a_cooler_1-delta_T_b_cooler_1)/(log(delta_T_a_cooler_1/delta_T_b_cooler_1)); %logarithmic temperature
Area_cooler_1=q_cooler_1/(U_cooler_1*delta_T_lm_cooler_1*F_cooler_1); %calculating the area needed 

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_cooler_1.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Cooler_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_cooler_1=q_cooler_1/1000*8000*1; %operation cost for cooler when cooled with a coolant at inlet T=Tc_in_coolant_1

fprintf('Total area needed is %4.0f m2 and investment cost is %7.0f SEK. Energy required is %5.0f kW and the operation cost per year is %9.0f \n', Area_cooler_1, Cost_Cooler_1, q_cooler_1/1000, Drift_cooler_1)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,TH_ut_cooler_1)

%% Condensor 1
format shortG
fprintf('CONDENSOR 1\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,T(end))
%Data
%Area calculations
U_cooler_1=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_cooler_1=1;

%Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler_1=[(FA0*(1-XA(end)));
       (FB0+FA0*XA(end));
       (FH0+FA0*XA(end));
       (FW0)];       %[mol/s] Butane-Butene-H2-H2O  

% Setting the desired out temperatures and the current in temperatures
TH_in_cooler_1=319;      %[K] Initial temperature of mixture
TH_ut_cooler_1=319; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_cooler_1=273-20;       %Inlet temperature of coolant   
TC_ut_cooler_1=273+34;        %Outlet temperature of coolant     

%Energy transfer calculations
Havg_cooler_1=F_mol_cooler_1(4,1)/sum(F_mol_cooler_1)*18.4*1000+(F_mol_cooler_1(1,1)+F_mol_cooler_1(2,1)+F_mol_cooler_1(3,1))/sum(F_mol_cooler_1)*22.3*1000; %finding the average enthalpy of vaporization for the mixture isobutane/water 

q_cooler_1=sum(F_mol_cooler_1)*Havg_cooler_1;        %[J/s] Energy transfer, effect needed
delta_T_a_cooler_1=TH_in_cooler_1-TC_ut_cooler_1;
delta_T_b_cooler_1=TH_ut_cooler_1-TC_in_cooler_1;
delta_T_lm_cooler_1=(delta_T_a_cooler_1-delta_T_b_cooler_1)/(log(delta_T_a_cooler_1/delta_T_b_cooler_1)); %logarithmic temperature

Area_cooler_1=q_cooler_1/(U_cooler_1*delta_T_lm_cooler_1*F_cooler_1); %area needed

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant
Cost_Year_B=a_s+b_s*Area_cooler_1.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Condensor_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_condensor_1=q_cooler_1/1000*8000*1;

fprintf('Total area needed is %3.0f m2 and total cost is %7.0f SEK. Energy required is %5.0f kW and the operation cost per year is %9.0f \n', Area_cooler_1, Cost_Condensor_1, q_cooler_1/1000,Drift_condensor_1)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,TH_ut_cooler_1)

%% Destillation 1, buten-vatten
clear x
clear y
clear i
clear n
fprintf('DESTILLATION 1\n')
fprintf('Inflow: %2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n',FA0*(1-XA(end)),FB0+FA0*XA(end),FH0+FA0*XA(end),FW0,TH_ut_cooler_1)
%1 = buten, 2 = vatten 

%Antoine constants for degC, mmHg, log10

A1=15.7564; B1=2132.42; C1=-33.15 ;%buten
A2=18.3036; B2=3816.44; C2=-46.13; %vatten
%total pressure
P =760;  %mmHg

tb1=B1/(A1-log(P))-C1;
tb2=B2/(A2-log(P))-C2;

x1 = linspace(0,1,1000);
Tstart=(tb1+tb2)/2;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tbnonideal(T,x1(i),x2,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = exp(A1-B1./(Tb(i)+C1));
    P1 = P01.*x1(i);
    y1(i) = P1./P;
end
figure(5)
hold on
axis([0 1 min(Tb) max(Tb)])
plot(x1, Tb)
plot(y1, Tb)
xlabel('x1,y1','Fontsize',15)
ylabel('T [K]','Fontsize',15)
title('Jämviktsdiagram för isobutan/vatten system','Fontsize',15)
legend('Mättad vätska','Mättad ånga')

hold off

figure(6)
hold on
plot(x1,y1)
plot(x1,x1,'red')
xlabel('x','Fontsize',15)
ylabel('y','Fontsize',15)
title('Ångfasssammansättnig mot vätskefassammansättning för isobutan/vatten system','Fontsize',12)
legend('Jämviktssammansättning','Diagonal')
% Molar flows 
F_cooler_1 = sum(F_mol_cooler_1); xF=(sum((F_mol_cooler_1(1:3))))/F_cooler_1; xD=0.99; xB=0.01; %KÃ¤nt inflÃ¶de F o xF, Ã¶nskade sammansÃ¤ttningar xD(toppen) och xB=azeotrop(botten)
A=[1 1 F_cooler_1; xD xB xF*F_cooler_1]; Flows=rref(A); %rÃ¤knar ut toppflÃ¶de och bottenflÃ¶de med total samt komponentbalans
D=Flows(1,3); B=Flows(2,3);
q_cooler_1 = 1; %previous cooler cooled the flow to a saturated liquid at current pressure and composition

%Calculate yF
Tbf = fsolve(@(T)find_Tbnonideal(T,xF,1-xF,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tbf+C1));
P1 = P01.*xF;
yF= P1./P; 

%Calculate R
LV = (xD - yF)/(xD - xF);
Rmin = LV/(1 - LV);
R = 1.5*Rmin;

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V+F_cooler_1*(q_cooler_1-1);
Lbar = L + F_cooler_1*q_cooler_1; 

%Initial temperature estimation and starting composition bottom of tower
Tstart = (tb1 + tb2)/2;
x(1) = xB;
Tb = fsolve(@(T)find_Tbnonideal(T,xB,1-xB,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tb+C1));
P1 = P01.*xB;
y(1)= P1./P;
% Botten
i = 1; 
while x<xF
    i = i + 1;
    x(i)=Vbar/Lbar*y(i-1) + B/Lbar*xB;
    y(i)=idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i));
    
end
% toppen
while y<xD
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F_cooler_1*xF);
    y(i)=idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i));
    i = i + 1;
end
n=i-1;
real=n/0.7;
% Torndimensioner

ts=0.45;
H=(real+1)*ts;

%Medelmolmassor
ML=x(1)*56.11+(1-x(1))*18.015; 
MV=y(1)*56.11+(1-y(1))*18.015;

%Temperatur i botten
Tbb = fsolve(@(T)find_Tbnonideal(T,x(1),1-x(1),1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
%Temperatur i toppen
Ttop = fsolve(@(T)find_Tbnonideal(T,xD,1-xD,1,1,A1,B1,C1,A2,B2,C2,P),Tstart, options);
%Medeldensiteter

MVrho=((ML/1000)*(P*133.322368))/(8.314*Tbb); %Kg/m^3
MLrho=(x(1))*588+(1-x(1))*997; %kg/mÂ³

FLV=(Lbar*ML)/(Vbar*MV)*sqrt(MVrho/MLrho);
CF=0.29; %flooding constant from diagram
sigma=72.8; %taken from some page water at 20Celcius
FST=(sigma/20)^0.2;
C=CF*FST;
Uf=C*sqrt(((MLrho-MVrho)/MVrho));
ada=0.1+(FLV-0.1)/9;
DT=sqrt((4*V*(MV/1000))/(0.8*(Uf/3.28)*pi*(1-ada)*MVrho));

% Totalkondensor och Ã¥terkokare
Hvap1=22.5*1000; %J/mol isobutene
Hvap2=44200; %J/mol vatten vid 20 grader C, 42911 vid 50 grader
Havgtop=xD*Hvap1+(1-xD)*Hvap2;
Havgbot=x(1)*Hvap1+(1-x(1))*Hvap2;

%condenser
Qc=D*(R+1)*Havgtop; %Joule/s

%reboiler
Qr=Vbar*Havgbot; %Joule/s

% COST DESTILLATIOn
%Pressure vesel

Density=7900;       %[kg/m^3] kolstål from kurs PM 

Dim=DT;     %[m] diameter of container 
L=H;        %[m] length of container 
S_max=118*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_k=11600; 
b_v_k=34; 
n_v_k=0.85; 

Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;    %[$] Cost calculations for different distillation trays
Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);              %[$] Accounts for inflation
Cost_Dist_v_k=Cost_Year_A_v_k*10;        %[SEK]

%Cost Trays   Valve trays=_v 
%Constants
a_v=210; 
b_v=400; 
n_v=1.9;
S_dist=Dim;      %[m] diameter of trays

Cost_Year_B_v=a_v+b_v*S_dist.^n_v;      %[$] 
Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation
Cost_Distillation_v=Cost_Year_A_v*10*ceil(real);       %[SEK] Cost of all trays
Cost_Dist_1=Cost_Dist_v_k+Cost_Distillation_v;      %Total cost
fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler_1(1)/(sum((F_mol_cooler_1(1:3))))*B*(xB),F_mol_cooler_1(2)/(sum((F_mol_cooler_1(1:3))))*B*(xB),F_mol_cooler_1(3)/(sum((F_mol_cooler_1(1:3))))*B*(xB),B*(1-xB),Tbb)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler_1(1)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(2)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(3)/(sum((F_mol_cooler_1(1:3))))*D*(xD),D*(1-xD),Ttop)
fprintf('Diameter of tower is %1.2f m and height is %1.2f m with %1.0f real trays. Total cost of tower is %6.0f SEK\n',Dim,H,ceil(real),Cost_Dist_1)
fprintf('Energy requirement in reboiler %1.2f MW. Energy requirement in condenser %1.2f MW\n\n\n',Qr*10^-6, Qc*10^-6)
m2=[(F_mol_cooler_1(1)/(sum((F_mol_cooler_1(1:3))))*B*(xB)); %butan
         (F_mol_cooler_1(2)/(sum((F_mol_cooler_1(1:3))))*B*(xB));  %buten
             (F_mol_cooler_1(3)/(sum((F_mol_cooler_1(1:3))))*B*(xB)); %H
             (B*(1-xB));];

%% Reboiler dest 1
format shortG
fprintf('REBOILER DEST 1\n')
%Data
%Area calculations
U_reboiler_1=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_reboiler_1=1;

% Setting the desired out temperatures and the current in temperatures
TH_in_reboiler_1=210+273;      %[K] Initial temperature of mixture
TH_ut_reboiler_1=273+150; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_reboiler_1=Tbb; 
TC_ut_reboiler_1=Tbb;

delta_T_a_reboiler_1=TH_in_reboiler_1-TC_ut_reboiler_1;
delta_T_b_reboiler_1=TH_ut_reboiler_1-TC_in_reboiler_1;
delta_T_lm_reboiler_1=(delta_T_a_reboiler_1-delta_T_b_reboiler_1)/(log(delta_T_a_reboiler_1/delta_T_b_reboiler_1));

Area_reboiler_1=Qr/(U_reboiler_1*delta_T_lm_reboiler_1*F_reboiler_1);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_reboiler_1.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_reboiler_dest_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_reboiler_dest_1=Qr/1000*8000*0.16;

fprintf('Total area needed is %3.0f m2 and investment cost is %6.0f SEK and the operation cost per year is %8.0f \n\n\n', Area_reboiler_1, Cost_reboiler_dest_1,Drift_reboiler_dest_1)

%% Condensor dest 1
format shortG
fprintf('CONDENSOR DEST 1\n')
%Data
%Area calculations
U_condensor_dest_1=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_condensor_dest_1=1;

% Setting the desired out temperatures and the current in temperatures
TH_in_condensor_dest_1=Ttop;      %[K] Initial temperature of mixture
TH_ut_condensor_dest_1=Ttop; 
TC_in_condensor_dest_1=273-30; 
TC_ut_condensor_dest_1=273-13;

delta_T_a_condensor_dest_1=TH_in_condensor_dest_1-TC_ut_condensor_dest_1;
delta_T_b_condensor_dest_1=TH_ut_condensor_dest_1-TC_in_condensor_dest_1;
delta_T_lm_condensor_dest_1=(delta_T_a_condensor_dest_1-delta_T_b_condensor_dest_1)/(log(delta_T_a_condensor_dest_1/delta_T_b_condensor_dest_1));

Area_condensor_dest_1=Qr/(U_condensor_dest_1*delta_T_lm_condensor_dest_1*F_condensor_dest_1);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_condensor_dest_1.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);

Cost_condensor_dest_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_condensor_dest_1=Qc/1000*8000*2.25;

fprintf('Total area needed is %3.0f m2 and investment cost is %f SEK and the operation cost per year is %f \n\n\n', Area_condensor_dest_1, Cost_condensor_dest_1,Drift_condensor_dest_1)

%% Reboiler 1
format shortG
fprintf('REBOILER 1\n')
%Data
%Area calculations
U_reboiler=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_reboiler=1;

% Setting the desired out temperatures and the current in temperatures
TH_in_reboiler=273+180;      %[K] Initial temperature of mixture
TH_ut_reboiler=273+100; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_reboiler=Ttop; 
TC_ut_reboiler=Ttop;

delta_T_a_reboiler=TH_in_reboiler-TC_ut_reboiler;
delta_T_b_reboiler=TH_ut_reboiler-TC_in_reboiler;
delta_T_lm_reboiler=(delta_T_a_reboiler-delta_T_b_reboiler)/(log(delta_T_a_reboiler/delta_T_b_reboiler));

Q=D*Havgtop; %effect needed based on the average enthalpy of vaporization for the mixture calculated in Destillation 1

Area_reboiler=Q/(U_reboiler*delta_T_lm_reboiler*F_reboiler);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_reboiler.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_reboiler_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_reboiler_1=Qr/1000*8000*0.16;

fprintf('Total area needed is %2.2f m2 and investment cost is %f SEK and the operation cost per year is %8.0f \n\n\n', Area_reboiler, Cost_reboiler_1,Drift_reboiler_1)

%% Compressor 
fprintf('COMPRESSOR\n')
fprintf('Inflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', F_mol_cooler_1(1)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(2)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(3)/(sum((F_mol_cooler_1(1:3))))*D*(xD),D*(1-xD),Ttop)
format shortG

T_in=Ttop; %inlet temperature

F_mol_compressor=[F_mol_cooler_1(1)/sum((F_mol_cooler_1(1:3)))*Flows(1,3)*xD;
       F_mol_cooler_1(2)/sum((F_mol_cooler_1(1:3)))*Flows(1,3)*xD;
       F_mol_cooler_1(3)/sum((F_mol_cooler_1(1:3)))*Flows(1,3)*xD;
       Flows(1,3)*(1-xD)];
   
F_mass_compressor=[F_mas(F_mol_compressor(1,1),1);
        F_mas(F_mol_compressor(2,1),2);
        F_mas(F_mol_compressor(3,1),3);
        F_mas(F_mol_compressor(4,1),4)]; %[kg/s]

CP_matrix_cooler_1=[Cp_new(T_in,1,1,1,1,1);
           Cp_new(T_in,2,2,2,2,2);
           Cp_new(T_in,3,3,3,3,3);
           Cp_new(T_in,4,4,4,4,4)];     %[J/Kg*K] Cp each component
  
% Calculations
C_matrix=F_mass_compressor.*CP_matrix_cooler_1;      
C_tot=sum(C_matrix);        %[W/K] Summan av m*cp för alla komponenter i flödet

P_in=1*10^5;     %[Pa] Ingående tryck till kompressorerna

P_ut=50*10^5;     %[Pa] Utgående tryck
eta_is=0.8;     %[] Isentropverkningsgrad
R=8.314;        %[J/mol*K]

R_matrix=[R/M(1,1);
          R/M(2,1);
          R/M(3,1);
          R/M(4,1)];        %Butane-Butene-H2-H2O in J/kg*K
  
Cv_matrix=[CP_matrix_cooler_1(1,1)-R_matrix(1,1);
           CP_matrix_cooler_1(2,1)-R_matrix(2,1);
           CP_matrix_cooler_1(3,1)-R_matrix(3,1);
           CP_matrix_cooler_1(4,1)-R_matrix(4,1)];       %%Butane-Butene-H2-H2O
       
kappa_matrix=[CP_matrix_cooler_1(1,1)/Cv_matrix(1,1);
              CP_matrix_cooler_1(2,1)/Cv_matrix(2,1);
              CP_matrix_cooler_1(3,1)/Cv_matrix(3,1);
              CP_matrix_cooler_1(4,1)/Cv_matrix(4,1)];
          

kappa=sum(kappa_matrix)./4;     %Divided by number of components (average)
       
%beräkningsfunktionen
%Pressure increase per step
P_step = (P_ut/P_in)^(1/3);  %[]
%Temperature out for each step in isentrope compression
T_ut_is = T_in*P_step^((kappa-1)/kappa);  %[K] 
%Actual temperature out from every compression step.
T_ut = T_in + (T_ut_is-T_in)/eta_is; %[K] 
%Required compressor effekt for one compresson step.
W = C_tot*(T_ut-T_in); %[W] 
%Total required compression effect (3 steg).
W_tot = 3*W; %[W] 
%Required cooling effekt for 1 cooler inbetween steps 
Q_kyl = C_tot*(T_ut-T_in);%[W] 
%Total required cooling effect for all coolers between steps (2)
Q_kyl_tot = 2*Q_kyl; %[W] 
%Coolant temperature
T_kv = 273.15-15; %[K] 
%Maximal temperature that the coolant may be heated to
T_kv_max = 273.15-10; %[K] 
%Logarithmic average temperature difference
deltaTlm = ((T_in-T_kv)-(T_ut-T_kv_max))/log((T_in-T_kv)/(T_ut-T_kv_max)); %[]
%U-value for coolerbetween steps (gas-liquid)
Ukyl = 200; %[W/(m2K)] 
%Heat exchange unit for 1 cooler between steps
A_kyl = Q_kyl/(Ukyl*deltaTlm); %[m2] 
%Total heat exchange area for coolers 
A_kyl_tot = 2*A_kyl; %[m2] 

% Utdata
W_tot;        %[W] Totalt effektbehov för kompressionen
Q_kyl;      %[W] Kylbehov i mellankylare
A_kyl_tot;     %[m2] Total värmeväxlararea för mellankylare
T_ut;    %[K] Utgående temperatur
% 3. COST: Compressor      (Centrifugal compressor)        %[kW]
a=580000;
b=20000;
n=0.6;
S_comp=W_tot/1000;       %[kW]
Cost_Year_B=a+b*S_comp.^n;
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Compressor_1=Cost_Year_A*10;

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*A_kyl_tot.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Cooler_comp_1=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_kost_komp=W_tot/1000*8000*0.3;
Drift_kost_kyl=Q_kyl/1000*8000*1;
fprintf('Outflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K and pressure %1.0f MPa\n', F_mol_cooler_1(1)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(2)/(sum((F_mol_cooler_1(1:3))))*D*(xD),F_mol_cooler_1(3)/(sum((F_mol_cooler_1(1:3))))*D*(xD),D*(1-xD),T_ut,P_ut*10^-6)
fprintf('Investment cost of compressor is %8.0f SEK and investment cost for coolers is %6.0f SEK. Operation cost for compressor %8.0f SEK and the operation cost for the coolers is %8.0f \n\n\n',Cost_Compressor_1,Cost_Cooler_comp_1,Drift_kost_komp,Drift_kost_kyl)

%% Cooler 2 
format shortG
fprintf('COOLER 2\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour\nat a temperature of %3.0f K \n',F_mol_compressor(1,1),F_mol_compressor(2,1),F_mol_compressor(3,1),F_mol_compressor(4,1),T_ut)
%Data
%Area calculations
U_cooler_2=200;       %[W/m2] KVärmegenomgångstal from kurs PM
F_cooler_2=1;

%Fetching and converting molar flow to mass flow
%Molar flow
F_mol_cooler_2=F_mol_compressor;       %[mol/s] Butane-Butene-H2-H2O  

%Mass flow 
F_mass_cooler_2=[F_mas(F_mol_cooler_2(1,1),1);
        F_mas(F_mol_cooler_2(2,1),2);
        F_mas(F_mol_cooler_2(3,1),3);
        F_mas(F_mol_cooler_2(4,1),4)];      %[kg/s] Butane-Butene-H2-H2O

% Setting the desired out temperatures and the current in temperatures
TH_in_cooler_2=T_ut;      %[K] Initial temperature of mixture
TH_ut_cooler_2=253; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_cooler_2=273-30;       %Inlet temperature of coolant
TC_ut_cooler_2=273-10;        %Outlet temperature of coolant

TH_medel_cooler_2=(TH_in_cooler_2+TH_ut_cooler_2)/2;       

CP_matrix_cooler_2=[Cp_new(TH_medel_cooler_2,1,1,1,1,1);
           Cp_new(TH_medel_cooler_2,2,2,2,2,2);
           Cp_new(TH_medel_cooler_2,3,3,3,3,3);
           Cp_new(TH_medel_cooler_2,4,4,4,4,4)];     %[J/Kg*K] Cp each component
     
%Energy transfer calculations
q_matrix_cooler_2=F_mass_cooler_2(:,1).*CP_matrix_cooler_2(:,1).*(TH_in_cooler_2-TH_ut_cooler_2);      %Energy transfer matrix       
q_cooler_2=sum(q_matrix_cooler_2);       %[J/s] Energy transfer

delta_T_a_cooler_2=TH_in_cooler_2-TC_ut_cooler_2;
delta_T_b_cooler_2=TH_ut_cooler_2-TC_in_cooler_2;
delta_T_lm_cooler_2=(delta_T_a_cooler_2-delta_T_b_cooler_2)/(log(delta_T_a_cooler_2/delta_T_b_cooler_2));

Area_cooler_2=q_cooler_2/(U_cooler_2*delta_T_lm_cooler_2*F_cooler_2);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_cooler_2.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Cooler_2=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_cooler_2=q_cooler_2/1000*8000*1;

fprintf('Total area needed is %3.0f m2 and investment cost is %6.0f SEK. Energy required is %4.0f kW and the operation cost per year is %8.0f \n', Area_cooler_2, Cost_Cooler_2, q_cooler_2/1000, Drift_cooler_2)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour at a temperature of %3.0f K \n \n \n',F_mol_compressor(1,1),F_mol_compressor(2,1),F_mol_compressor(3,1),F_mol_compressor(4,1),TH_ut_cooler_2)

%% Condensor 2 (partial)
format shortG
fprintf('CONDENSOR 2\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour\nat a temperature of %3.0f K \n',F_mol_compressor(1,1),F_mol_compressor(2,1),F_mol_compressor(3,1),F_mol_compressor(4,1),TH_ut_cooler_2)
%Data
%Area calculations
U_cooler_1=850;       %[W/m2] KVärmegenomgångstal from kurs PM
F_cooler_1=1;

%Molar flow
F_mol_cooler_1=F_mol_compressor;   %[mol/s] Butane-Butene-H2-H2O  

% Setting the desired out temperatures and the current in temperatures
TH_in_cooler_1=253;      %[K] Initial temperature of mixture
TH_ut_cooler_1=253; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_cooler_1=273-40;       %Inlet temperature of coolant     
TC_ut_cooler_1=273-30;        %Outlet temperature of coolant     

%Energy transfer calculations
A=[1 1 sum(F_mol_cooler_1);0.98707 0.085665 F_mol_cooler_1(3,1)];% solving flow balances and component balances based on a isobutane/H2 system to see how much of the flow is liquid at the current state
flows=rref(A);
L=flows(2,3);
Havg=0.085665*0.44936*1000+(1-0.085665)*22.6*1000; %finding average enthalpy of vaporization based on a isobutane/hydrogen system

q_condensor=L*Havg; %effect required in the condensor
delta_T_a_cooler_1=TH_in_cooler_1-TC_ut_cooler_1;
delta_T_b_cooler_1=TH_ut_cooler_1-TC_in_cooler_1;
delta_T_lm_cooler_1=(delta_T_a_cooler_1-delta_T_b_cooler_1)/(log(delta_T_a_cooler_1/delta_T_b_cooler_1));

Area_cooler_1=q_condensor/(U_cooler_1*delta_T_lm_cooler_1*F_cooler_1); %area needed

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_cooler_1.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_Condensor_2=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_condensor_2=q_condensor/1000*8000*1;

fprintf('Total area needed is %3.0f m2 and total cost is %7.0f SEK. Energy required is %4.0f kW and the operation cost per year is %8.0f \n', Area_cooler_1, Cost_Condensor_2, q_condensor/1000,Drift_condensor_2)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour\nat a temperature of %3.0f K \n \n \n',F_mol_compressor(1,1),F_mol_compressor(2,1),F_mol_compressor(3,1),F_mol_compressor(4,1),TH_ut_cooler_2)

%% Återflöde
fprintf('ÅTERFLÖDE\n')
fprintf('Inflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour at a temperature of %3.0f K \n',F_mol_compressor(1,1),F_mol_compressor(2,1),F_mol_compressor(3,1),F_mol_compressor(4,1),TH_ut_cooler_1)
isobutan=F_mol_cooler_2(1)/(sum((F_mol_cooler_2(1:3))))*D*(xD)+0.26161; %adding the recirculated amount 
isobuten=F_mol_cooler_2(2)/(sum((F_mol_cooler_2(1:3))))*D*(xD)+0.78615;
h2=F_mol_cooler_2(3)/(sum((F_mol_cooler_2(1:3))))*D*(xD)+82.35;
water=D*(1-xD)+0.016801;
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %0.2f mol/s water vapour at a temperature of %3.0f K \n \n \n',isobutan,isobuten,h2,water,TH_ut_cooler_2);

%% Flash
%Antoine constants
Butanin=isobutan; Butenin=isobuten; H2in=h2; Win=water;
fprintf('FLASH\n')
fprintf('Inflow:%1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', Butanin,Butenin,H2in,Win,TH_ut_cooler_2)
A1 = [13.6333, 164.90, 3.19]; %hydorgen
A2 = [15.5381, 2032.73, -33.15]; %butan

%total pressure
P = 50*760;  %mmHg

x1 = linspace(0, 1);
options = optimset('Display', 'off');
%For loop which will calculate the temperature which gives y1+y2=1 for each x1 for the
%system based on Raoults law and assumption of ideal mixture
for i=1:length(x1) 
x2=@(x1)(1-x1(i));

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3)))); %calculating Psat for (1) and (2)

y1=@(T)(P1sat(T).*x1(i)./P); y2=@(T)(P2sat(T).*x2(x1)./P); %using Raoults law to calculate y1 and y2

ysum=@(T)(y1(T)+y2(T)-1); %creating the function to be solved 

T(i)=fsolve(ysum,200,options);%finding the T which gives ysum=0

y11(i)=y1(T(i)); %calculating these ys and saving in a vector using the functions created in line 18 and 20
end

%plotting the temperature vs x-y for the system
figure(7)
plot(x1,T,'blue')
hold on
plot(y11,T,'red')
axis([0 1 min(T) max(T)])
xlabel('x','Fontsize',15)
ylabel('T [K]','Fontsize',15)
title('Temperatur mot sammansättning för isobutan/H_{2} system','Fontsize',15)
grid on
legend('Mättad vätska', 'Mättad ånga')

T=253; %temperature at which the flash will operate
fun=@(x1)(find_Tbflash(P,T,A1,A2,x1)-T);  %function to be solved
x11=fzero(fun,0.1); %solving the function to find the component fraction of the liquid phase
[T,y1]=find_Tbflash(P,T,A1,A2,x11); %using the found x-value to calculate the component fraction in the vapour phase

Ftot=Butanin+Butenin+H2in+Win; 
z1in=H2in/Ftot; %defining the inlet component fraction 
A=[1 1 Ftot;y1 x11 Ftot*z1in];%assembling the equations to be solved
flows=rref(A); %solving the equations by row reduction. First line gives flow for vapour, second gives liquid flow
hold on
line([z1in z1in],[0 T]); line([x11 y1], [T T]); line([x11 x11], [0 T], 'Linestyle', '--'); line([y1 y1], [0 T], 'Linestyle', '--') %drawing the flash operation into the diagram
legend('Liquid phase composition', 'Vapor phase composition')

% Dimensions

MV=y1*2*1.00784+(1-y1)*58.12; %average Mr of vapor g/mol
ML=x11*2*1.00784+(1-x11)*58.12; %average Mr of liquid g/mol

rhoV=((MV/1000)*(P*133.322368))/(8.314*T); %density of vapor kg/m3
rhoL=x11*70.85+(1-x11)*563; %density of liquid kg/m3

ut=0.07*sqrt((rhoL-rhoV)/rhoV); %m/s

L=flows(2,3)*(ML/1000)/rhoL; %liquid flow (mol/s)*(kg/mol)/(kg/m3)=m3/s
V=flows(1,3)*(MV/1000)/rhoV; %vapor flow (mol/s)*(kg/mol)/(kg/m3)=m3/s

DT=sqrt((4*V)/(pi*0.15*ut));

HL=(L*10*60)/((pi/4)*DT^2);

H=HL+1.5*DT;

%cost
Density=7900;       %[kg/m^3] kolstål from kurs PM 

Dim=DT;     %[m] diameter of container (1-3.5)
L=H;        %[m] length of container 
S_max=88.9*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
E=1;     %svets verkningsgrad (= 1)
P=50*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass_flash=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_k=11600;        
b_v_k=34;       
n_v_k=0.85;        

Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass_flash.^n_v_k;         %[$] Cost calculations for different distillation trays
Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);     %[$] Accounts for inflation
Cost_Flash_v_k_1=Cost_Year_A_v_k*10;        %[SEK]

fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n',Butanin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),Butenin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),x11*flows(2,3),Win/(Butanin+Butenin+Win)*(1-x11)*flows(2,3),253)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', Butanin/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),Butenin/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),y1*flows(1,3),Win/(Butanin+Butenin+Win)*(1-y1)*flows(1,3),253)
fprintf('Dimensions: diameter is %1.2f m and height is %2.1f m. The cost of the flash is %7.0f\n\n\n',Dim,H,Cost_Flash_v_k_1)

%% Strypventil
Put=3.3*760; 
Tut=T; %antar att temperaturen ej förändras i strypventilen då vi rör oss inom underkyldvätska området

%% Destillation 2 butan-buten
%1 = buten, 2 = butan
butanin=Butanin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3); butenin=Butenin/(Butanin+Butenin+Win)*(1-x11)*flows(2,3); h2in=x11*flows(2,3); win=Win/(Butanin+Butenin+Win)*(1-x11)*flows(2,3);
fprintf('DESTILLATION 2\n')
fprintf('Inflow:%1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', butanin,butenin,h2in,win,T)
clear x
clear y
clear i
clear n
%Wilson parameters
W12 = 0.48584; 
W21 = 1.64637;

%Antoine constants for degC, mmHg, log10

A2=15.7564; B2=2132.42; C2=-33.15 ;%buten
A1=15.5381; B1=2032.73; C1=-33.15;%butan
%total pressure
P =Put;  %mmHg

tb1=B1/(A1-log(P))-C1;
tb2=B2/(A2-log(P))-C2;

x1 = linspace(0,1,1000);
Tstart=(tb1+tb2)/2;  %temperature at which to start the search

for i = 1:length(x1)
    x2 = 1-x1(i);

    %activity coefficients at x1
    gamma2 = exp(-log(x1(i)+W12*x2)+x2.*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
    gamma1 = exp(-log(x2+W21*x1(i))-x1(i).*((W12./(x1(i)+W12*x2))-(W21./(W21*x1(i)+x2))));
   
    %use fsolve function to find bubble point temperature (Tb) for x1
    %find_Tb is a function we need to create that will check if a certain value of T satisfies y1+y2-1=0 
    %current value of x1 and other constants are passed to find_Tb
    options = optimset('Display', 'off');
    Tb(i) = fsolve(@(T)find_Tbnonideal(T,x1(i),x2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
    
    P01 = exp(A1-B1./(Tb(i)+C1));
    P1 = gamma1.*P01.*x1(i);
    y1(i) = P1./P;
end
figure(9)
hold on
axis([0 1 min(Tb) max(Tb)])
plot(x1, Tb)
plot(y1, Tb)
xlabel('x1,y1','Fontsize',15)
ylabel('T [K]','Fontsize',15)
title('Jämviktsdiagram för isobutan/vatten system','Fontsize',15)
legend('Mättad vätska','Mättad ånga')
hold off

figure(10)
hold on
plot(x1,y1)
plot(x1,x1,'red')
xlabel('x','Fontsize',15)
ylabel('y','Fontsize',15)
title('Ångfasssammansättnig mot vätskefassammansättning för isobutan/vatten system','Fontsize',12)
legend('Jämviktssammansättning','Diagonal')
% Flows
F = butanin+butenin+h2in+win; xF2=(butanin+h2in)/(butanin+butenin+h2in+win); xD2=0.55; xB2=0.06; %Känt inflöde F o xF, önskade sammansättningar xD(toppen) och xB=azeotrop(botten)
A=[1 1 F; xD2 xB2 xF2*F]; Flows=rref(A); %räknar ut toppflöde och bottenflöde med total samt komponentbalans
D=Flows(1,3); B=Flows(2,3);


%Calculate yF
gamma2 = exp(-log(xF2+W12*(1-xF2))+(1-xF2).*((W12./(xF2+W12*(1-xF2)))-(W21./(W21*xF2+(1-xF2)))));
gamma1 = exp(-log((1-xF2)+W21*xF2)-xF2.*((W12./(xF2+W12*(1-xF2)))-(W21./(W21*xF2+(1-xF2)))));
Tbf = fsolve(@(T)find_Tbnonideal(T,xF2,1-xF2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tbf+C1));
P1 = gamma1.*P01.*xF2;
yF= P1./P; 
lambda=21.5*1000;
q_cooler_1=(95.21*(Tbf-T)+(lambda))/lambda;

%Calculate R
LV = (xD2 - yF)/(xD2 - xF2);
Rmin = LV/(1 - LV);
R = 2*Rmin;

%Flows through tower
L = R*D;
V = D*(R+1);
Vbar = V+F*(q_cooler_1-1);
Lbar = L + F*q_cooler_1;

%Initial temperature estimation and starting composition bottom of tower
Tstart = (tb1 + tb2)/2;
x(1) = xB2;
gamma2 = exp(-log(xB2+W12*(1-xB2))+(1-xB2).*((W12./(xB2+W12*(1-xB2)))-(W21./(W21*xB2+(1-xB2)))));
gamma1 = exp(-log((1-xB2)+W21*xB2)-xB2.*((W12./(xB2+W12*(1-xB2)))-(W21./(W21*xB2+(1-xB2)))));
Tb = fsolve(@(T)find_Tbnonideal(T,xB2,1-xB2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);
P01 = exp(A1-B1./(Tb+C1));
P1 = gamma1.*P01.*xB2;
y(1)= P1./P;
% Botten
i = 1; 
while x<xF2
    i = i + 1;
    x(i)=Vbar/Lbar*y(i-1) + B/Lbar*xB2;
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21);
    
end
% toppen
while y<xD2
    x(i) = V/L*y(i - 1) + 1/L*(B*x(1)-F*xF2);
    y(i)=nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x(i), W12, W21);
    i = i + 1;
end
n=i-1;
real=n/0.7;
% Torndimensioner

ts=0.45;
H=(real+1)*ts;

%Medelmolmassor
ML=x(1)*58.12+(1-x(1))*56.11;
MV=y(1)*58.12+(1-y(1))*56.11;

%Temperatur i toppen
gamma2 = exp(-log(xD2+W12*(1-xD2))+(1-xD2).*((W12./(xD2+W12*(1-xD2)))-(W21./(W21*xD2+(1-xD2)))));
gamma1 = exp(-log((1-xD2)+W21*xD2)-xD2.*((W12./(xD2+W12*(1-xD2)))-(W21./(W21*xD2+(1-xD2)))));
Ttop = fsolve(@(T)find_Tbnonideal(T,xD2,1-xD2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P),Tstart, options);

%Medeldensiteter

MVrho=((ML/1000)*(P*133.322368))/(8.314*Tbb); %Kg/m^3
MLrho=(x(1))*563+(1-x(1))*588; %Kg/m^3

FLV=(Lbar*ML)/(Vbar*MV)*sqrt(MVrho/MLrho);
CF=0.28; %flooding constant from diagram
sigma=15.8; %taken from some sketchy page isobutene at 20Celcius
FST=(sigma/20)^0.2;
C=CF*FST;
Uf=C*sqrt(((MLrho-MVrho)/MVrho));
ada=0.1+(FLV-0.1)/9;
DT=sqrt((4*V*(MV/1000))/(0.8*(Uf/3.28)*pi*(1-ada)*MVrho));

% Totalkondensor och återkokare
Hvap1=21.5*1000; %J/mol isobutan
Hvap2=22.5*1000; %J/mol isobutene
Havgtop=xD2*Hvap1+(1-xD2)*Hvap2;
Havgbot=x(1)*Hvap1+(1-x(1))*Hvap2;

%condenser

Qc=D*(R+1)*Havgtop; %Joule/s

%reboiler
Qr=Vbar*Havgbot; %Joule/s

%Cost calculations
%Cost Pressure vesel
Density=7900;       %[kg/m^3] from kurs PM 

Dim=DT;     %[m] diameter of container (1-3.5)
L=H;        %[m] length of container 
S_max=88.9*10^6;     %[N/m^2]maximalt tillåtna materialspänning 
P=1.1*1*10^5;     %konstructionstrycket (10% större än arbetstryck) (Pa)

t=(P*Dim)/(2*S_max*E-1.2*P);        %[m] Reactor wall thickness 

V_full=pi*(((Dim/2)+t)^2)*(L+2*t);      %[m^3]Volume of full tank
V_inner=pi*(((Dim/2))^2)*(L);     %[m^3] Volume of inner tank (air)
Volume_container=V_full-V_inner;        %[m^3] Volume of shell
Shell_mass_dist=Density*Volume_container;        %[kg] Mass of shell (S)

a_v_k=17400;
b_v_k=79;
n_v_k=0.85;

Cost_Year_B_v_k=a_v_k+b_v_k*Shell_mass.^n_v_k;      %[$] Cost calculations for different distillation trays
Cost_Year_A_v_k=Cost_Year_B_v_k*(CEPCI_Year_A/CEPCI_Year_B);       %[$] Accounts for inflation
Cost_Dist_v_k=Cost_Year_A_v_k*10;        %[SEK]

%Cost Trays: Valve trays=_v 
%Constants
a_v=210;
b_v=400;
n_v=1.9;
S_dist=Dim;      %[m] diameter of trays
Cost_Year_B_v=a_v+b_v*S_dist.^n_v;     %[$] Cost calculations for different distillation trays

Cost_Year_A_v=Cost_Year_B_v*(CEPCI_Year_A/CEPCI_Year_B);        %[$] Accounts for inflation
Cost_Distillation_v=Cost_Year_A_v*10*ceil(real);       %[SEK]
Cost_Dist_2=Cost_Dist_v_k+Cost_Distillation_v;

fprintf('Flows at bottom: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n',butanin/(butanin+h2in)*xB2*B,butenin/(butenin+win)*(1-xB2)*B,h2in/(butanin+h2in)*(xB2)*B,win/(butenin+win)*(1-xB2)*B, Tb)
fprintf('Flows at top: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.0f K\n', butanin/(butanin+h2in)*D*xD2, butenin/(butenin+win)*(1-xD2)*D,h2in/(butanin+h2in)*(xD2)*D,win/(butenin+win)*(1-xD2)*D, Ttop)
fprintf('Diameter of tower is %1.2f m and height is %1.2f m with %1.0f real trays. Cost of tower is %f SEK\n',Dim,H,ceil(real),Cost_Dist_2)
fprintf('Energy requirement in reboiler %1.2f MW. Energy requirement in condenser %1.2f MW\n\n\n',Qr*10^-6, Qc*10^-6)

%% Reboiler dest 2
format shortG
fprintf('REBOILER DEST 2\n')
%Data
%Area calculations
U_reboiler_2=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_reboiler_2=1;

% Setting the desired out temperatures and the current in temperatures
TH_in_reboiler_2=273+180;      %[K] Initial temperature of mixture
TH_ut_reboiler_2=273+150; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_reboiler_2=Tb; 
TC_ut_reboiler_2=Tb;

delta_T_a_reboiler_2=TH_in_reboiler_2-TC_ut_reboiler_2;
delta_T_b_reboiler_2=TH_ut_reboiler_2-TC_in_reboiler_2;
delta_T_lm_reboiler_2=(delta_T_a_reboiler_2-delta_T_b_reboiler_2)/(log(delta_T_a_reboiler_2/delta_T_b_reboiler_2));

Area_reboiler_2=Qr/(U_reboiler_2*delta_T_lm_reboiler_2*F_reboiler_2);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_reboiler_2.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_reboiler_dest_2=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_reboiler_dest_2=Qr/1000*8000*0.16;

fprintf('Total area needed is %3.0f m2 and investment cost is %6.0f SEK and the operation cost per year is %8.0f \n\n\n', Area_reboiler_2, Cost_reboiler_dest_2,Drift_reboiler_dest_2)

%% Condensor dest 2
format shortG
fprintf('CONDENSOR DEST 2\n')
%Data
%Area calculations
U_condensor_dest_2=1000;       %[W/m2] KVärmegenomgångstal from kurs PM
F_condensor_dest_2=1;

% Setting the desired out temperatures and the current in temperatures
TH_in_condensor_dest_2=Ttop;      %[K] Initial temperature of mixture
TH_ut_condensor_dest_2=Ttop; %fsolve(@(T)find_Tbnonideal(T,(F_mol_cooler(2,1))/(sum(F_mol_cooler)),1-((F_mol_cooler(2,1))/(sum(F_mol_cooler))),1,1,A1,B1,C1,A2,B2,C2,760),320, options);      %[K] Temperature after heat final exchange
TC_in_condensor_dest_2=273-5; 
TC_ut_condensor_dest_2=273;

delta_T_a_condensor_dest_2=TH_in_condensor_dest_2-TC_ut_condensor_dest_2;
delta_T_b_condensor_dest_2=TH_ut_condensor_dest_2-TC_in_condensor_dest_2;
delta_T_lm_condensor_dest_2=(delta_T_a_condensor_dest_2-delta_T_b_condensor_dest_2)/(log(delta_T_a_condensor_dest_2/delta_T_b_condensor_dest_2));

Area_condensor_dest_2=Qr/(U_condensor_dest_2*delta_T_lm_condensor_dest_2*F_condensor_dest_2);

% Cost of heat exchange units (uses Cooler.m)
% 1. COST: Cooler        (Floating head, 10 m^2 < size < 1000 m^2)       %[m^2]

a_s=32000;        %Cost constant
b_s=70;       %Cost constant
n_s=1.2;      %Equipment constant

Cost_Year_B=a_s+b_s*Area_condensor_dest_2.^n_s;      %[USD $ in 2010]  Total cost of heat exchange units (convert to SEK in 2020)
% Conversion for cooler cost
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_condensor_dest_2=Cost_Year_A*10; %Conversion from Dollars to SEK
Drift_condensor_dest_2=Qc/1000*8000*1;

fprintf('Total area needed is %3.0f m2 and investment cost is %7.0f SEK and the operation cost per year is %9.0f \n\n\n', Area_condensor_dest_2, Cost_condensor_dest_2,Drift_condensor_dest_2)

%% Strypning 
fprintf('STRYPNING\n')
fprintf('Inflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.3f K and pressure of %1.1f bar\n', butanin/(butanin+h2in)*D*xD, butenin/(butenin+win)*(1-xD)*D,h2in/(butanin+h2in)*(xD)*D,win/(butenin+win)*(1-xD)*D, Ttop,3.3)
R=8.314;
Tc=408.1; Pc=36.48; w=0.177; %kritiska parameterar
T1=Ttop;   %temp in (K)
P1=Put;    % tryck in (bar)
Tr=T1/Tc;
Pr=P1/Pc;
T2=300;     %gissning av ut temp
Tmedel=(T1+T2)/2;
cpBA=@(Tmedel)CP(1,1)+CP(1,2).*Tmedel+CP(1,3).*Tmedel.^2+CP(1,4).*Tmedel.^3;     cpBA=cpBA(Tmedel); %cp-beräkning
B0=0.083-(0.422/(Tr^1.6));
B1=0.139-(0.172/(Tr^4.2));
H=R*Tc*Pr*(B0-(Tc*(0.422*1.6/(Tr^2.6)))+(w*(B1-(Tr*(0.172*4.2/(Tr^5.2))))))/1000000; %beräkning av entalpi innan strypning 
T2=T1+(H/cpBA);
fprintf('Outflow: %1.2f mol/s isobutane, %1.2f mol/s isobutene, %1.2f mol/s H2, %0.2f mol/s water at a temperature of %3.3f K and pressure of %1.0f bar\n\n\n', butanin/(butanin+h2in)*D*xD, butenin/(butenin+win)*(1-xD)*D,h2in/(butanin+h2in)*(xD)*D,win/(butenin+win)*(1-xD)*D, T2,1)

%% Mix 1 destströmmarna vid 1 bar
fprintf('Inflow: The bottom flow from dest 1 and the top flow from dest 2 are mixed.\n')
m1=[(butanin/(butanin+h2in)*xD2*D); %butan                                   
         (butenin/(butenin+win)*(1-xD2)*D);  %buten
         (h2in/(butanin+h2in)*(xD2)*D); %H
         (win/(butenin+win)*(1-xD2)*D)];    %water
         
cpA=@(T)CP(1,1)+CP(1,2).*T+CP(1,3).*T.^2+CP(1,4).*T.^3;
cpB=@(T)CP(2,1)+CP(2,2).*T+CP(2,3).*T.^2+CP(2,4).*T.^3;
cpH=@(T)CP(3,1)+CP(3,2).*T+CP(3,3).*T.^2+CP(3,4).*T.^3;
cpW=@(T)CP(4,1)+CP(4,2).*T+CP(4,3).*T.^2+CP(4,4).*T.^3;  

y=@(T) m2(1,1)*cpA(T)*(Tbb-T)+m2(2,1)*cpB(T)*(Tbb-T)+m2(3,1)*cpH(T)*(Tbb-T)+m2(4,1)*cpW(T)*(Tbb-T)-m1(1,1)*cpA(T)*(T-Ttop)-m1(2,1)*cpB(T)*(T-Ttop)-m1(3,1)*cpH(T)*(T-Ttop)-m1(4,1)*cpW(T)*(T-Ttop);
Tut=fsolve(y,(Ttop+Tbb)/2,options); %uttempereratur för flödet

m=m1+m2;
fprintf('Outflow: %2.2f mol/s isobutane, %2.2f mol/s isobutene, %2.2f mol/s H2, %4.0f mol/s water at a temperature of %3.0f K \n\n\n',m(1),m(2),m(3),m(4),Tut)

%% Mix 2 färskt mixas med destströmmarna påväg in i startugnen
fprintf('Inflow: the flow from mix 1 is mixed with the fresh feed containing isobutane and water at 283 K.\n')
m4=[(FA0-m(1,1)); %färskt inflöde med butan och vatten
    0;
    0;
    FW0-m(4,1)];
y=@(t) m(1,1).*cpA(t).*(Tut-t)+m(2,1).*cpB(t).*(Tut-t)+m(3,1).*cpH(t).*(Tut-t)+m(4,1).*cpW(t).*(Tut-t)-m4(1,1).*cpA(t).*(t-283)-m4(4,1).*cpW(t).*(t-283);      
T_ut=fsolve(y,300,options);
m0=m4+m;
fprintf('Outflow: %2.2f mol/s isobutane, %2.2f mol/s isobutene, %2.2f mol/s H2, %4.0f mol/s water at a temperature of %3.0f K \n\n\n',m0(1),m0(2),m0(3),m0(4),T_ut)

%% Ugn 0 
Tin=T_ut; Tut=950; Tmedel=(Tin+Tut)/2;        %[K]

F_mass_ugn_1=[F_mas(m0(1),1);
              F_mas(m0(2),2);
              F_mas(m0(3),3);
              F_mas(m0(4),4)];  
          
    CP_matrix_ugn_1=[Cp_new(Tmedel,1,1,1,1,1);
                     Cp_new(Tmedel,2,2,2,2,2);
                     Cp_new(Tmedel,3,3,3,3,3);
                     Cp_new(Tmedel,4,4,4,4,4)];     %[J/Kg*K] Cp each component

QBA=F_mass_ugn_1(1)*CP_matrix_ugn_1(1,1)*(Tut-Tin);
QBE=F_mass_ugn_1(2)*CP_matrix_ugn_1(2,1)*(Tut-Tin);
QH=F_mass_ugn_1(3)*CP_matrix_ugn_1(3,1)*(Tut-Tin);
QW=F_mass_ugn_1(4)*CP_matrix_ugn_1(4,1)*(Tut-Tin);      %[J/s]
Qtot=QBA+QBE+QH+QW;
         
chi=0.8;        %Verkningsgrad
Qheat=Qtot/0.8;     %[J/s] Energi som krävs för uppvärmning 
Driftko_ugn_1=Qheat/1000*8000*0.2; %[SEK] Driftskostnad för 1 år, naturgas

%Cost Ugn: Cylindrical

a=43000;
b=111000;
n=0.8;
S_ugn_1=Qheat/(10^6); %[kW]
Cost_Year_B=a+b*S_ugn_1.^n;
Cost_Year_A=Cost_Year_B*(CEPCI_Year_A/CEPCI_Year_B);
Cost_ugn_1=Cost_Year_A*10;
fprintf('OVEN 0 \n')
fprintf('Inflow: %2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour at a temperature of %3.0f K \n',m0(1),m0(2),m0(3),m0(4),T_ut)
fprintf('Effect needed is %4.0f kW and the cost per year is %8.0f SEK då naturgas används. Investeringskostnaden blir %7.0f SEK \n', Qheat/1000,Driftko_ugn_1,Cost_ugn_1)
fprintf('Outflow:%2.2f mols/s isobutane, %0.2f mol/s isobutene, %0.2f mol/s hydrogen gas, %3.0f mol/s water vapour\nat a temperature of %3.0f K \n\n\n ',m0(1),m0(2),m0(3),m0(4),Tut)

%% Total Cost
% Investment cost
I_Cost=Cost_Dist_h_304_with_kat_1+Cost_ugn_2+Cost_Dist_h_304_with_kat_2+Cost_Cooler_1+Cost_Condensor_1+Cost_Dist_1+Cost_reboiler_dest_1+Cost_condensor_dest_1+Cost_reboiler_1+Cost_Compressor_1+Cost_Cooler_comp_1+Cost_Cooler_2+Cost_Condensor_2+Cost_Flash_v_k_1+Cost_Dist_2+Cost_reboiler_dest_2+Cost_condensor_dest_2+Cost_ugn_1;
% Operation cost
D_Cost=Driftko_ugn_2+Drift_cooler_1+Drift_condensor_1+Drift_reboiler_dest_1+Drift_condensor_dest_1+Drift_reboiler_1+Drift_kost_komp+Drift_kost_kyl+Drift_cooler_2+Drift_condensor_2+Drift_reboiler_dest_2+Drift_condensor_dest_2+Driftko_ugn_1;
fprintf('Total Cost \n')
fprintf('The total investment cost for the process is %f SEK and the operation cost is %f SEK per year \n\n\n',I_Cost,D_Cost)

close all

%% Functions
function dYdW=ode_func(W,Y,HR,P,CP,FA0,FB0,FW0,FH0)
%function file containing differential equations

XA=Y(1); T=Y(2); %plockar ut Xa och T från inputmatris

R=8.314; %J/K/mol
K1=22.90; %bar^-1.5
K2=7.56; %/bar
Ea=141e3; %J/mol

FA=(FA0)*(1-XA); FB=FB0+FA0*XA; FH=FH0+FA0*XA; FW=FW0; %räknar nya flöden
Ftot=FA+FB+FH+FW; %beräknar totala flödet
PA=FA/Ftot*P; PB=FB/Ftot*P; PH=FH/Ftot*P; %beräknar partial trycken

%beräknar nya reaktionskonstanter
k=0.0596*exp((Ea/R).*(1./(550+273)-1./T)); %mol/(kg cat)/s/bar
Ke=2.1*10^7*exp(-122/(R*T)); %bar

%definierar funktioner för att beräkna värmekapacitet med användning av de
%konstanter som matades in i matrisen CP
cpA=@(T)CP(1,1)+CP(1,2).*T+CP(1,3).*T.^2+CP(1,4).*T.^3;
cpB=@(T)CP(2,1)+CP(2,2).*T+CP(2,3).*T.^2+CP(2,4).*T.^3;
cpH=@(T)CP(3,1)+CP(3,2).*T+CP(3,3).*T.^2+CP(3,4).*T.^3;
cpW=@(T)CP(4,1)+CP(4,2).*T+CP(4,3).*T.^2+CP(4,4).*T.^3;

%korrekterar reaktionsentalpin för den rådande temperaturen
delHr=HR-integral(cpA,298,T)+integral(cpB,298,T)+integral(cpH,298,T);

%beräknar ny reaktionshastighet
r=(k*(PA-(PB*PH)/Ke))/(1+K1*PB*PH^0.5+(K2*PH)^0.5);

%samlar ekvationerna som ska lösas
dYdW=[r/FA0 %mole balance
    (r*(-delHr))/(FA*cpA(T)+FB*cpB(T)+FH*cpH(T)+FW*cpW(T))]; %heat balance
end
function F_mass=F_mas(F_mol,m)

%Molar mass
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

F_mass=(F_mol.*M(m,1));       %Butane-Butene-H2-H2O in kg/s
end
function Cp=Cp_new(T,a,b,c,d,m)      %Calculates Cp for each component

%Molar mass
M=[58.12*10^-3;
   56.1*10^-3;
   2*1.00784*10^-3;
   18.01528*10^-3];     %Butane-Butene-H2-H2O in kg/mol

%Specific heat capacities 
Cp_initial=[1.39 0.3847 -1.846e-04 2.895e-08;
            16.05 0.2804 -1.091e-04 9.098e-09;
            27.14 0.009274 -1.3813e-05 7.645e-09;
            32.24 0.001924 1.055e-05 -3.596e-09];       %Butane-Butene-H2-H2O in J/mol/K

Cp=(Cp_initial(a,1)+Cp_initial(b,2).*T+Cp_initial(c,3).*T^2+Cp_initial(d,4).*T^3)./M(m,1);   %based on c=(a+b*T+c*T^2+d*T^3)/M;        
end
function res = find_Tbnonideal(T,x1,x2,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P)

%Use T,x1,gamma1,gamma2,A1,B1,C1,A2,B2,C2,P
%to calculate y1 and y2

P01 = exp(A1-B1./(T+C1));
P02 = exp(A2-B2./(T+C2));

P1 = gamma1.*P01.*x1; 
P2 = gamma2.*P02.*x2;

y1 = P1./P;
y2 = P2./P;

res = y1+y2-1;
end
function y1 = idealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x1)

x2=@(x1)(1-x1);

%gamma2= 1; %@(x1)(exp(-log(x1+W12.*x2(x1))+(x2(x1).*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));
%gamma1= 1;  %@(x1)(exp(-log(x2(x1)+W21.*x1)-(x1.*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));

P1sat=@(T)(exp(A1-B1./(T+C1))); P2sat=@(T)(exp(A2-B2/(T+C2)));

y1=@(T)(P1sat(T).*x1./P); y2=@(T)(P2sat(T).*x2(x1)./P);

sumy=@(T)(y1(T)+y2(T)-1);

T=fsolve(sumy,Tstart,optimset('Diagnostics','off', 'Display','off'));
y1=y1(T);
y2=y2(T);
end
function [T,y1]=find_Tbflash(P,T0,A1,A2,x1)
%function which finds T at which y1+y2=1 
x2=@(x1)(1-x1);

P1sat=@(T)(exp(A1(1)-A1(2)./(T+A1(3)))); P2sat=@(T)(exp(A2(1)-A2(2)/(T+A2(3))));

y1=@(T)(P1sat(T).*x1./P); y2=@(T)(P2sat(T).*x2(x1)./P);

sumy=@(T)(y1(T)+y2(T)-1);

T=fsolve(sumy,T0,optimset('Diagnostics','off', 'Display','off'));
y1=y1(T);
y2=y2(T);
end
function y1 = nonidealTb(P,Tstart,A1,B1,C1,A2,B2,C2,x1, W12, W21)

x2=@(x1)(1-x1);

gamma2=@(x1)(exp(-log(x1+W12.*x2(x1))+(x2(x1).*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));
gamma1=@(x1)(exp(-log(x2(x1)+W21.*x1)-(x1.*(W12./(x1+W12.*x2(x1))-W21./(W21.*x1+x2(x1))))));

P1sat=@(T)(exp(A1-B1./(T+C1))); P2sat=@(T)(exp(A2-B2/(T+C2)));

y1=@(T)(gamma1(x1).*P1sat(T).*x1./P); y2=@(T)(gamma2(x1).*P2sat(T).*x2(x1)./P);

sumy=@(T)(y1(T)+y2(T)-1);

T=fsolve(sumy,Tstart,optimset('Diagnostics','off', 'Display','off'));
y1=y1(T);
y2=y2(T);
end

