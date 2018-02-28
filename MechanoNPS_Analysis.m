
Nstart =1;
FiletoRead=50;

file_num = [Nstart,FiletoRead]; % file numbers, 1 is the first #, 2 is the number of files to read after
Fs=60;  %sampling frequency [kHz]
Ndown=10;    %down sampling period
buffer_input = 8000;  %buffer to measure baseline current at the both end of signal
MaxSignalSize= 12000; %Maximum length of signal to set matrix
gap = 4000;
error=0.08;
Thr_in=1e-4;      %initial value for threshold finding

%% read a new data file
load ('/Users/junghyun/Documents/Research/Collaboration/Mark Lab/Exp Data/180227_TEST8/mcf7-houseAir-dev2.mat','data');
data_temp=data(2,:);
mydata=data_temp(Nstart*Fs*1000:(Nstart+FiletoRead)*Fs*1000);

clear data_temp data_whole data 

y=mydata;
%y=cat(1,mydata{:}); % concatenate data
y_smoothed=fastsmooth(y',100,1,1); % perform smoothing
y_detrend=detrend(y_smoothed); % remove trend

clear myfilename mydata i k file_name file_id FiletoRead Nstart thresholds,

%% find outliers && group by clumps

out = isoutlier(y_detrend);
out_ind = 1;
group = zeros(2,length(out));
group_ind = 1;
while out_ind < length(out)
    if out(out_ind) == 1
        group(1,group_ind) = out_ind; % first number in the group
        group(2,group_ind) = 1; %counter for consecutive 1s
        out_ind = out_ind + 1;
            while out_ind < length(out) && out(out_ind) == 1 %checks if the following data is 0
                group(2,group_ind) = group(2,group_ind) + 1; %adds to the count
                out_ind = out_ind + 1;
            end
    end
    out_ind = out_ind + 1;
    group_ind = group_ind + 1;
end
group(:, all(group==0))=[]; % removes 0, 0 data

clear out out_len out_ind,

%% combines signals if the gap btwn the signals is small

group_ind = 1;
while group_ind < length(group)
    last_index = group(1,group_ind) + group(2,group_ind); %last point of group
    if group(1, group_ind + 1) - last_index < gap % checks if the gap btwn groups is small. If so, we combine the groups
        group(2, group_ind) = group(1, group_ind + 1) + group(2, group_ind + 1) - group(1, group_ind); %last index of 2nd group - first index of 1st = length
        group(:,group_ind + 1) = [];
    else
        group_ind = group_ind + 1;
    end
end

clear gap group_ind

%% creates "signals" matrix
% takes each group and puts it in matrix "signals" with each column being a cell signal data

signal_detrend = zeros(MaxSignalSize,length(group)); % creates a matrix with size (group x (number of data points + buffer)) 
signal_smooth = zeros(MaxSignalSize,length(group));
buffer=buffer_input;

for xx=1:length(group)
    if group(2,xx)<MaxSignalSize

        for i = 1:length(group)
            start_signal = group(1,i) - (buffer_input/2);

            if start_signal < 1
                buffer = buffer - start_signal;
                start_signal = 1;
            end
            end_signal = start_signal + group(2,i) + buffer;
            
            if (end_signal-start_signal)< MaxSignalSize && end_signal < length(y_smoothed)
                data1 = y_smoothed(start_signal:end_signal);
                data2 = y_detrend(start_signal:end_signal);
            end
            
            signal_smooth(1:length(data1),i) = data1;
            signal_detrend(1:length(data2),i) = data2;
            signal_down=downsample(signal_detrend,Ndown);
            signal_length(1,i)=size(signal_down(:,i),1);
            buffer = buffer_input;
        end
    end
end
clear start_signal end_signal buffer gap last_index group_ind data1 data2 signal_max xx i y y_detrend y_smoothed,

%% derivative of signal
global ii
for ii=1:length(group)
    
    ThrX=Thr_in; %derivative threshold value for node-pore
    
    signal_temp=signal_down(1:signal_length(1,ii),ii);
    signal_temp(signal_temp==0)=[];
    signal_diff=diff(signal_temp);
    clear signal_temp
    
    pks=pksFinder(signal_diff,ThrX); %Find peaks in signal_diff
    
    hFig=figure(ii); %Prinf figure to check peak finding
    set(hFig,'Position', [0 300 1200 400])
    
        subplot(1,2,1)
        plot(signal_detrend(:,ii));
        set(gca,'FontSize',15)
        
        subplot(1,2,2)
        plot(signal_diff)
        hold on
        plot(pks(:,1),pks(:,2),'ob');
        hold off
        set(gca,'FontSize',15)
    
    if isempty(pks)~=1&&length(pks)>3
        
        fprintf([ 'ENTER will return peaks,\n' ... 
                        'S will reture 0\n' ...
                        'N will set a new initial threshold \n' ...
                        'Or, just enter new threshold for node-pore \n' ...
                        ]);
        OK = input('---\n','s');

        switch OK
            case []
                Thr(ii)=ThrX;

            case{'s','S'}
                Thr(ii)=[0];

            case{'n','N'}
                Thr(ii)=[100];
                a=str2double(OK);
                clear a
                fprintf(['Now enter new threshold for the contraction channel \n']);
                ThrX=input('---\n');
                pks=pksFinder(signal_diff,ThrX); %Find peaks in signal_diff

                Thr(ii)=ThrX;

                hFig=figure(ii); %Prinf figure to check peak finding
                set(hFig,'Position', [0 300 1200 400])
                    subplot(1,2,1)
                    plot(signal_detrend(:,ii));
                    set(gca,'FontSize',15)

                    subplot(1,2,2)
                    plot(signal_diff)
                    hold on
                    plot(pks(:,1),pks(:,2),'ob');
                    hold off
                    set(gca,'FontSize',15)

                    fprintf(['Press ENTER \n']);
                    a=input('---\n'); 

            otherwise
                Thr1=str2double(OK);
                Thr(ii)=100;
                fprintf(['Now enter new threshold for the contraction channel \n']);
                Thr2=input('---\n');
                Thresholds=[Thr1, Thr2];

                signal_diff(abs(signal_diff) < Thresholds(1)) = 0;
                k=0;
                for i=1:length(signal_diff)
                    if signal_diff(i)<-Thresholds(2)
                       k=i+1;
                       for j=k:length(signal_diff)
                            if  abs(signal_diff(j))<Thresholds(2)
                                signal_diff(j)=0;
                            elseif signal_diff(j)>Thresholds(2)
                                break;
                            end
                       end
                    end
                end

                for i=1:size(signal_diff,1)
                    a(i)=i;
                    b(i)=signal_diff(i);
                end
                clear pks
                aa=transpose(a);
                bb=transpose(b);
                pks_temp=[aa,bb];
                pks1=pks_temp;
                pks1(any(pks_temp==0,2),:) = [];
                pks=pks1;
                %pks=pks1(1:14,:);
                clear pks_temp pks1 a b aa bb

                                set(hFig,'Position', [0 300 1200 400])
                                subplot(1,2,1)
                                plot(signal_detrend(:,ii));
                                set(gca,'FontSize',15)

                                subplot(1,2,2)
                                plot(signal_diff)
                                hold on
                                plot(pks(:,1),pks(:,2),'ob');
                                hold off
                                set(gca,'FontSize',15)

                                fprintf(['Press ENTER \n']);
                                a=input('---\n'); 

        end

        close all

        locs_temp=pks(:,1);
        pks_temp=pks(:,2);
        pks_temp(find(isnan(pks_temp)))=[];
        locs_temp(find(isnan(locs_temp)))=[];

        PKS=[locs_temp pks_temp];

        OUT(:,ii)=measureSignal(signal_smooth(:,ii), signal_detrend(:,ii), signal_diff, PKS, Ndown, file_num, Fs, Thr(ii), error);
        clear pks_temp locs_temp pks PKS
    end
    close all
    clear pks_temp locs_temp pks PKS
end

clear a i j k Thresholds ThrX Thr1 Thr2 OK signal_detrend signal_diff signal_down signal_length signal_smooth

global fj fi

if isempty(fj)
    fj=1;
end

savename=sprintf('MeasureOut_trial%03d.mat',fj);
fj=fj+1;

if isempty(fi)
    fi=1;
end

if exist('OUT')
OUT( :, all( ~any( OUT ), 1 ) ) = [];
    for i=1:size(OUT,2) 
            MeasureOut(:,fi)=OUT(:,i);
            save(savename,'MeasureOut')
            fi=fi+1;
    end
end

close all
clear buffer_input error group i MaxSignalSize Out Ndown Thr OUT
