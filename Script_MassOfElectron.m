%/////////////////////////////////////////////////////////
% By: Lee Allers                                         /
%For: Numerical Computation, 2016                        /
%     University of New Mexico                           /
%NOTE: None of my scripts are built to be robust, they   /
%      are merely an implementation of a given set of    /
%      data or instructions!                             /
%/////////////////////////////////////////////////////////
clc; close all; clear all

Erest = 511; %KeV
Background = xlsread('data\Background.csv', 'A21:G1044'); Background(:,2) = [];
Ba133 = xlsread('data\Ba133.csv', 'A21:G1044'); Ba133(:,2) = []; Ba133 = [Ba133(:,1) Ba133(:,2) - Background(:,2)]; Ba133(Ba133 < 0) = 0; Ba133(:,2) = Ba133(:,2).*10^(-1);
Co60 = xlsread('data\Co60.csv', 'A21:G1044'); Co60(:,2) = []; Co60 = [Co60(:,1) Co60(:,2) - Background(:,2)]; Co60(Co60 < 0) = 0;
Cs137 = xlsread('data\Cs137.csv', 'A21:G1044'); Cs137(:,2) = []; Cs137 = [Cs137(:,1) Cs137(:,2) - Background(:,2)]; Cs137(Cs137 < 0) = 0;
Na22 = xlsread('data\Na22.csv', 'A21:G1044'); Na22(:,2) = []; Na22 = [Na22(:,1) Na22(:,2) - Background(:,2)]; Na22(Na22 < 0) = 0;


%Cs137 Backscatter Peak at
%CsBy = max(Cs137([70 140],2)); %y value of backscatter peak
%CsBx = Cs137(Cs137(:,2) == CsBy); %x value of backscatterpeak
%CsBE = CsBx * Eco; %Converted x value to energy
iX = 310; iX2 = 365;
ResY = max(Cs137([iX:iX2],2)); Cs137_x_energy = Cs137(Cs137(:,2) == ResY);
iX = 550; iX2 = 600;
ResY = max(Co60([iX:iX2],2)); Co60_x_energyA = Co60(Co60(:,2) == ResY);
iX = 650; iX2 = 700;
ResY = max(Co60([iX:iX2],2)); Co60_x_energyB = max(Co60(Co60(:,2) == ResY));
iX = 200; iX2 = 300;
ResY = max(Na22([iX:iX2],2)); Na22_x_energyA = Na22(Na22(:,2) == ResY);
iX = 600; iX2 = 700;
ResY = max(Na22([iX:iX2],2)); Na22_x_energyB = Na22(Na22(:,2) == ResY);
iX = 120; iX2 = 220;
ResY = max(Ba133([iX:iX2],2)); Ba133_x_energy = Ba133(Ba133(:,2) == ResY);

Ec = 1.97; %Mean of the Energy Coefficients

Cs137_E = 662 / Cs137_x_energy;
Co60_EA = 1170 / Co60_x_energyA;
Co60_EB = 1330 / Co60_x_energyB;
Na22_EA = 511 / Na22_x_energyA;
Na22_EB = 1275 / Na22_x_energyB;
Ba133_E = 356 / Ba133_x_energy;


Defined_Energies = [ 356 , 511, 662, 1170, 1275, 1330 ];
y = [Ba133_x_energy Na22_x_energyA Cs137_x_energy Co60_x_energyA Na22_x_energyB Co60_x_energyB ];
P = polyfit(Defined_Energies,y,1);
Y = polyval(P,Defined_Energies);

string1 = sprintf('Ba133: %0.00f',y(1));
string2 = sprintf('Na22_{A}: %0.00f',y(2));
string3 = sprintf('Cs137: %0.00f',y(3));
string4 = sprintf('Co60_{A}: %0.00f',y(4));
string5 = sprintf('Na22_{B}: %0.00f',y(5));
string6 = sprintf('Co60_{B}: %0.00f',y(6));
string7 = sprintf('Fit: %.02fx + %.00f',P); 
%Plots the E versus channel
figure(1)
plot(Defined_Energies(1),y(1),'xr',Defined_Energies(2),y(2),'dg',Defined_Energies(3),y(3),'xr',Defined_Energies(4),y(4),'ok',Defined_Energies(5),y(5),'dg',Defined_Energies(6),y(6),'ok');hold on
plot(Defined_Energies,Y);
errorbar(Defined_Energies(1),y(1),12)
errorbar(Defined_Energies(2),y(2),18)
errorbar(Defined_Energies(3),y(3),6)
errorbar(Defined_Energies(4),y(4),13)
errorbar(Defined_Energies(5),y(5),11)
errorbar(Defined_Energies(6),y(6),9)

hold off
xlabel('Energy KeV'); ylabel('Channel #'); title('Calibration Parameters')
legend1 = legend(string1,string2,string3,string4,string5,string6,string7);
set(legend1,...
    'Position',[0.695828274646815 0.172812079185113 0.187120287524104 0.259999992630699]);


%clear val ResY idx iX iX2 pval y Y legend1 string1 string2 string3 string4 string5 string6 string7

Comp =@(x1,x2) ((x1 + x2) / 2) * Ec;
eRestMass =@(loc,comp) 2*Defined_Energies(loc)./((1/(1 - comp/Defined_Energies(loc))) - 1); %Rest mass of electron measured
%eRestMass =@(loc,comp) 2*Defined_Energies(loc)/( (Defined_Energies(loc)/comp) -1 );
eRestMass2 =@(loc,loc2,comp) (Defined_Energies(loc) + Defined_Energies(loc2))./((1/(1 - comp/Defined_Energies(loc))) - 1); %Rest mass of electron measured (averaged peaks of Co
EGamma_exp =@(loc,erest) Defined_Energies(loc)/(1 + (erest/(2*Defined_Energies(loc)))); %Experimental Backscatter Peak
EGamma_theory = Defined_Energies./(1 + (511./(2*Defined_Energies))); %Theory Backscatter Peak


%Mass of Electron

%//////////Ces133
idx = 240; idx2 = 275; loc = 3;
Cs137_Comp = Comp(idx,idx2);
Cs137_Erest = eRestMass(loc,Cs137_Comp);
idx = 300 ; idx2 = 400;
Cs137_fwhm = 100*fwhm(Cs137(idx:idx2,1),Cs137(idx:idx2,2))/662;

Cs137_error_Edge = 100 - EGamma_theory(3)/Cs137_Comp * 100;
Cs137_error_Mass = 100 - Cs137_Erest/Erest * 100;

%//////////Co60 (Averaged energies of 2 peaks)
idx = 450; idx2 = 570; loc = 4; loc2 = 6;
Co60_Comp = Comp(idx,idx2);
Co60_Erest = eRestMass2(loc,loc2,Co60_Comp);
%********Peak(s) FWHM
idx = 550 ; idx2 = 630;
Co60_1_fwhm = 100*fwhm(Co60(idx:idx2,1),Co60(idx:idx2,2))/1170;
idx = 630 ; idx2 = 720;
Co60_2_fwhm = 100*fwhm(Co60(idx:idx2,1),Co60(idx:idx2,2))/1330;
Co60_fwhm = (Co60_1_fwhm + Co60_2_fwhm)/2;

Co_error_Edge = 100 - Co60_Comp/((EGamma_theory(4) + EGamma_theory(6))/2) * 100;
Co_error_Mass = 100 - Co60_Erest/Erest * 100;

%//////////Na22 
idx = 140; idx2 = 220; loc = 2;
Na22_Comp = Comp(idx,idx2);
Na22_Erest = eRestMass(loc,Na22_Comp);

idx = 520; idx2 = 580; loc = 5;
Na22_CompB = Comp(idx,idx2);
Na22_ErestB = eRestMass(loc,Na22_CompB);
%********Peak(s) FWHM
idx = 230 ; idx2 = 292;
Na22_fwhm = 100*fwhm(Na22(idx:idx2,1),Na22(idx:idx2,2))/511;

Na22_error_Edge = 100 - EGamma_theory(2)/Na22_Comp * 100;
Na22_error_Mass = 100 - Na22_Erest/Erest * 100;
Na22_B_error_Edge = 100 - EGamma_theory(5)/Na22_CompB * 100;
Na22_B_error_Mass = 100 - Na22_ErestB/Erest * 100;


%//////////Ba133
idx = 95; idx2 = 125; loc = 1;
Ba133_Comp = Comp(idx,idx2);
Ba133_Erest = eRestMass(loc,Ba133_Comp);
%********Peak(s) FWHM
idx = 150 ; idx2 = 250;
Ba133_fwhm = 100*fwhm(Ba133(idx:idx2,1),Ba133(idx:idx2,2))/356;
Ba133_error_Edge = 100 - EGamma_theory(1)/Ba133_Comp * 100;
Ba133_error_Mass = 100 - Ba133_Erest/Erest * 100;


clear idx idx2


 BS =@(E) E/(2 + Erest/E);

Peaks = Ec*[BS(Cs137_x_energy),BS(Ba133_x_energy),BS(Co60_x_energyA),BS(Co60_x_energyB),BS(Na22_x_energyA),BS(Na22_x_energyB)];


xx = Cs137(:,1)*Ec;
figure ('units','normalized','position',[0 0 .6 .6]); 

%$$$$$ PLOT CESIUM-137 $$$$$
% subplot(4,2,1)
% plot(xx,Cs137(:,2))
% xlabel('Energy (KeV)')
% title('Cs137')

string = sprintf('FWHM: %.00f%%',Cs137_fwhm);

subplot(2,2,1)
plot(xx,Cs137(:,2))

tx = text(400,2000,string);
tx.HorizontalAlignment = 'center';    % set horizontal alignment to center
tx.VerticalAlignment = 'top';         % set vertical alignment to top
title('Cs137')
xlabel('Energy (KeV)')
xlim([0 ceil(400*Ec)])
l1 = line([EGamma_theory(3) EGamma_theory(3)],[0 1000],'Color', 'red');
l2 = line([Cs137_Comp Cs137_Comp],[0 1000],'Color', 'magenta');
l3 = line([Peaks(1) Peaks(1)],[0,1000],'Color','green');

string1 = sprintf('Comp Theory: %.02f KeV',EGamma_theory(3));
string2 = sprintf('Comp Experimental: %.02f KeV',Cs137_Comp);
string3 = sprintf('BS Theory: %.02f KeV',Peaks(1));
legend1 = legend([l1 l2 l3],{string1,string2,string3});
set(legend1,...
    'Position',[0.358217596214402 0.906121400648674 0.157118051933746 0.0563271590221075]);


%$$$$$ PLOT COBALT-60 $$$$$
% subplot(4,2,4)
% plot(xx,Co60(:,2))
% xlabel('Energy (KeV)')
% title('Co60')
string = sprintf('FWHM: %.00f%%',Co60_fwhm);
T = (EGamma_theory(4) + EGamma_theory(6))/2;
subplot(2,2,2)
plot(xx,Co60(:,2))
tx = text(670,600,string);
tx.HorizontalAlignment = 'center';    % set horizontal alignment to center
tx.VerticalAlignment = 'top';         % set vertical alignment to top
title('Co60')
xlabel('Energy (KeV)')
xlim([0 ceil(700*Ec)])
l1 = line([T T],[0 300],'Color', 'red');
l3 =  line([Co60_Comp Co60_Comp],[0 300],'Color', 'magenta');
l4 = line([Peaks(3) Peaks(3)],[0,300],'Color','green');
l5 = line([Peaks(4) Peaks(4)],[0,300],'Color','yellow');
string1 = sprintf('Theory: %.02f KeV',T);
string3 = sprintf('Experimental: %.02f KeV',Co60_Comp);
string4 = sprintf('BS Theory: %.02f KeV',Peaks(4));
string5 = sprintf('BS Theory: %.02f KeV',Peaks(3));

legend1 = legend([l1 l3 l4 l5],{string1, string3, string4, string5});
set(legend1,...
    'Position',[0.813946762881069 0.893004117433917 0.157118051933746 0.0817901212492107]);
clear T

%$$$$$ PLOT SODIUM-22 $$$$$
% subplot(4,2,7)
% plot(xx,Na22(:,2))
% xlabel('Energy (KeV)')
% title('Na22')
string = sprintf('FWHM: %.00f%%',Na22_fwhm);

subplot(2,2,3)
plot(xx,Na22(:,2))
tx = text(670,2000,string);
tx.HorizontalAlignment = 'center';    % set horizontal alignment to center
tx.VerticalAlignment = 'top';         % set vertical alignment to top
title('Na22')
xlabel('Energy (KeV)')
xlim([0 ceil(700*Ec)])
l1 = line([EGamma_theory(2) EGamma_theory(2)],[0 1000],'Color', 'red');
l2 = line([EGamma_theory(5) EGamma_theory(5)],[0 1000],'Color', 'red');
l3 = line([Na22_Comp Na22_Comp],[0 1000],'Color', 'magenta');
l4 = line([Na22_CompB Na22_CompB],[0 1000],'Color', 'magenta');
l5 = line([Peaks(5) Peaks(5)],[0,1000],'Color','green');
l6 = line([Peaks(6) Peaks(6)],[0,1000],'Color','yellow');

string1 = sprintf('Theory: %.02f KeV',EGamma_theory(2));
string2 = sprintf('Theory: %.02f KeV',EGamma_theory(5));
string3 = sprintf('Experimental: %.02f KeV',Na22_Comp); 
string4 = sprintf('Experimental: %.02f KeV',Na22_CompB); 
string5 = sprintf('BS Theory: %.02f KeV',Peaks(5));
string6 = sprintf('BS Theory: %.02f KeV',Peaks(6));

legend1 = legend([l1 l2 l3 l4 l5 l6],{string1,string2,string3,string4, string5,string6});
set(legend1,...
    'Position',[0.355613429728826 0.412294241626567 0.163194440641544 0.107253083476314]);


%$$$$$ PLOT Barium-133 $$$$$
% subplot(4,2,1)
% plot(xx,Ba133(:,2))
% xlabel('Energy (KeV)')
% title('Ba133')
string = sprintf('FWHM: %.00f%%',Ba133_fwhm);

subplot(2,2,4)
plot(xx,Ba133(:,2))
tx = text(225,7000,string);
tx.HorizontalAlignment = 'center';    % set horizontal alignment to center
tx.VerticalAlignment = 'top';         % set vertical alignment to top
title('Ba133')
xlabel('Energy (KeV)')
xlim([0 ceil(250*Ec)])
l1 = line([EGamma_theory(1) EGamma_theory(1)],[0 1000],'Color', 'red');
l2 = line([Ba133_Comp Ba133_Comp],[0 1000],'Color', 'magenta');
l3 = line([Peaks(2) Peaks(2)],[0,1000],'Color','green');

string1 = sprintf('Theory: %.02f KeV',EGamma_theory(1));
string2 = sprintf('Experimental: %.02f KeV',Ba133_Comp);
string3 = sprintf('BS Theory: %.02f KeV',Peaks(2));

legend1 = legend([l1 l2 l3],{string1,string2, string3});
set(legend1,...
    'Position',[0.830439818436625 0.426183129043736 0.157118051933746 0.0563271590221076]);

M_c = [510.1 506.9 510.2 511.9 520.6 508.5];
M_b = [506.4 513.5 510.2 512.9 520.6 510.9]; 
O_e = 1:length(M_c);
M_c_err = std(M_c);
M_b_err = std(M_b);
figure(3)
hold on
scatter(O_e,M_c)

scatter(O_e,M_b,'x')

l1 = line([1 O_e(end)],[511 511],'Color','red');
legend('Rest Mass Estimate (Compton)','Rest Mass Estiamte (Backscatter)','Actual Electron Mass')

%clear loc loc2 legend1 P tx xx string string1 string2 string3 string4 l1 l2 l3 l4 eRestMass eRestMass2 Ec  Erest G EGamma_exp Comp %EGamma_theory

eRestMass =@(g,comp) 2*g./((1/(1 - comp/g)) - 1); %Rest mass of electron measured
B = @(b,g) 2.*(b .* g)./(g - b) .* 1.97;


figure ('units','normalized','position',[0 0 .6 .6]); 
subplot(4,1,1)
plot(Cs137(:,1),Cs137(:,2))
subplot(4,1,2)
plot(Co60(:,1),Co60(:,2))
subplot(4,1,3)
plot(Na22(:,1),Na22(:,2))
subplot(4,1,4)
plot(Ba133(:,1),Ba133(:,2))



figure;
hold on
plot(Cs137(:,1),Cs137(:,2),'r')
plot(Ba133(:,1),Ba133(:,2),'g')
plot(Co60(:,1),Co60(:,2),'b')
plot(Na22(:,1),Na22(:,2),'m')
plot(Background(:,1),Background(:,2),'k')
whitebg(1,'k')
ylim([0 2000]);xlim([0 700])
xlabel('Channel')
ylabel('Counts')
legend(sprintf('^{137}Cs'),sprintf('^{133}Ba'),sprintf('^{60}Co'),sprintf('^{22}Na'),'Background Noise')
hold off
%/////////////////////////////////

Background = xlsread('data\Background.csv', 'A21:G1044'); Background(:,2) = [];
Ba133 = xlsread('data\Ba133.csv', 'A21:G1044'); Ba133(:,2) = []; Ba133 = [Ba133(:,1) Ba133(:,2) - Background(:,2)]; Ba133(Ba133 < 0) = 0;
Co60 = xlsread('data\Co60.csv', 'A21:G1044'); Co60(:,2) = []; Co60 = [Co60(:,1) Co60(:,2) - Background(:,2)]; Co60(Co60 < 0) = 0;
Cs137 = xlsread('data\Cs137.csv', 'A21:G1044'); Cs137(:,2) = []; Cs137 = [Cs137(:,1) Cs137(:,2) - Background(:,2)]; Cs137(Cs137 < 0) = 0;
Na22 = xlsread('data\Na22.csv', 'A21:G1044'); Na22(:,2) = []; Na22 = [Na22(:,1) Na22(:,2) - Background(:,2)]; Na22(Na22 < 0) = 0;
Ba133 = [Ec*Ba133(:,1) , Ba133(:,2)];
Co60 = [Ec*Co60(:,1) , Co60(:,2)];
Cs137 = [Ec*Cs137(:,1) , Cs137(:,2)];
Na22 = [Ec*Na22(:,1) , Na22(:,2)];




figure ('units','normalized','position',[0 0 .6 .6]); 
subplot(4,1,1)
plot(Cs137(:,1),Cs137(:,2))
subplot(4,1,2)
plot(Co60(:,1),Co60(:,2))
subplot(4,1,3)
plot(Na22(:,1),Na22(:,2))
subplot(4,1,4)
plot(Ba133(:,1),Ba133(:,2))


mco =@(E,E2) 1/( E/E2 - 1)/(2*E);
