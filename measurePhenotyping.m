clear all

%% Set variables 

D=27.65;     % effective D of node-pore[um]
L=6750;     % channel length [um]
w_np=22;    % channel width at node-pore
wc=10.5;    % width of the contraction channel
h=22.07;    % channel height
f=60;     %[kHz] sampling frequency

%% readout pulse
% load MeasureOut matrix from save file
% need to set correct FILE NUMBER
load('MeasureOut_trial066.mat') 

for i=1:size(MeasureOut,2)
    I=MeasureOut(3,i);    %baseline current
    dI=MeasureOut(4,i);   %current drop at node-pore
    dIc=MeasureOut(5,i);  %current drop at the contraction channel
    dT=MeasureOut(6,i);   %transit time at pore  
    dTc=MeasureOut(7,i);  %transit time at the contraction channel
    Tr=MeasureOut(8,i);   %recovery time
    
    d=(dI*L*(D^2)/(I+dI*L*0.8/D))^(1/3); %cell diameter
    e=(d-wc)./d;            %strain
    D2=D.*(wc./w_np).^0.5;  %effective D of cont. channel[um]
    
    dD=(dIc*L*(D2^2)/(I+dIc*L*0.8/D2))^(1/3);
    dD=pi.*(dD.^2)./(4*wc); %deformed diameter
    
    Phenotyping(:,i)=[MeasureOut(1,i); d; dT; dTc; dD; Tr]; %file #, cell diameter, transit time at pore[ms], transit time at cont. channel[ms], deformed diameter, recovery time
end

    Phenotyping=transpose(Phenotyping);
    save Phenotyping.mat

