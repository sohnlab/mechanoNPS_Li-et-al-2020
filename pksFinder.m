function pks=pksFinder(signal_diff,thr0)

    
    Thr0=thr0; %derivative threshold value for node-pore            
    [pks1,locs1]=findpeaks(signal_diff,'MinPeakProminence',Thr0); %find positive peaks
    [pks2,locs2]=findpeaks(-signal_diff,'MinPeakProminence',Thr0); %find negative peaks
    
    pks=[pks1;-pks2];   %combine peaks
    locs=[locs1;locs2]; %combine locations
    pks_temp=[locs, pks];
    
    [~,idx] = sort(pks_temp(:,1)); % sort just the first column
    sortedpks = pks_temp(idx,:);   % sort the whole matrix using the sort indices
    
    clear pks locs pks_temp idx
    
    TF1=isoutlier(sortedpks(:,2));  %find peaks at the cont. channel
    TF1 = double(TF1);
    k=0;
    for j=1:length(TF1)             %remove noise peaks at the cont. channel
        if TF1(j)==1
           k=j+1;
           for i=k:length(TF1)
                if TF1(i)==0
                     sortedpks(i,1)=NaN;
                     sortedpks(i,2)=NaN;
                elseif TF1(i)==1
                        break;
                end               
           end
           break;
        end
    end
    pks=sortedpks;
    
end