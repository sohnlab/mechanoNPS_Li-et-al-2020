clear all

%% Set variables 

D=27.65;     %effective D of node-pore[um]
L=6750;     %channel length [um]
w_np=22;    %channel width at node-pore
wc=10.5;
h=22.07;    %channel height
f=50;     %[kHz] sampling frequency

%% readout pulse
load('MeasureOut_trial066.mat')

for i=1:size(MeasureOut,2)
    I=MeasureOut(3,i);
    dI=MeasureOut(4,i);  
    dIc=MeasureOut(5,i);
    dT=MeasureOut(6,i);
    dTc=MeasureOut(7,i);
    Tr=MeasureOut(8,i);
    
    d=(dI*L*(D^2)/(I+dI*L*0.8/D))^(1/3); %cell diameter
    e=(d-wc)./d;            %strain
    D2=D.*(wc./w_np).^0.5;  %effective D of cont. channel[um]
    
    dD=(dIc*L*(D2^2)/(I+dIc*L*0.8/D2))^(1/3);
    dD=pi.*(dD.^2)./(4*wc); %deformed diameter
    
    
    Phenotyping_Ctrl(:,i)=[MeasureOut(1,i); d; dT; dTc; dD; Tr];


end

    Phenotyping_Ctrl=transpose(Phenotyping_Ctrl);
    save Phenotyping_Ctrl.mat

