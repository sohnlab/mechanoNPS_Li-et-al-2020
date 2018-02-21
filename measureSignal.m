function [OUT] = measureSignal (y_smoothed, y_detrend, y_diff, pks, Ndown, file_num, Fs, flag, error)
 
 global ii;
 
 if flag~=0
    backset = 30;
    A = ones(3,length(pks)); % preallocation
    k=0;
    for i=1:length(pks)
        if abs(pks(i,2))>0
            k=k+1;
            A(1,k)=pks(i,1)*Ndown; % array index of nonzero
            A(2,k)=pks(i,2); % nonzero value
            A(3,k)=(pks(i,1)-backset)*Ndown; % backward offset for baseline current
        end
    end
    clear tmp backset k 
%% remove error from A
    k=1;
    while (k < size(A,2))
        i=1;
        while (i < size(A,2))

            % Case 1: current and next both positive && next > current
            if A(2,i) > 0 && A(2,i+1) > 0 && ...
                    A(2,i) < A(2,i+1)
                % move next into current
                A(:,i)=A(:,i+1);

            % Case 2: current and next both positive && current > next
            elseif A(2,i) > 0 && A(2,i+1) > 0 && ...
                    A(2,i) > A(2,i+1)
                % move current into next
                A(:,i+1)=A(:,i);

            % Case 3: current and next both negative && current > next
            elseif A(2,i) < 0 && A(2,i+1) < 0 && ...
                    A(2,i) > A(2,i+1)
                % move next into current
                A(:,i)=A(:,i+1);

            % Case 4: current and next both negative && next > current
            elseif A(2,i) < 0 && A(2,i+1) < 0 && ...
                    A(2,i) < A(2,i+1)
                % move current into next
                A(:,i+1)=A(:,i);
            end

            i=i+1;
        end
        k=k+1;
    end

%% remove repeats in A

    stop = 10;
    i = 1;
if size(A,2)>3
    while (i < stop)
       if (A(:,i) == A(:,i+1))
           A(:,i) = [];
           stop = stop - 1;      
       else
           i = i + 1;
       end
    end
end

  %% rectangularize pulses
    y_down=downsample(y_detrend,Ndown);
    
    k=1;
    while (k <= 50)
        i = 1;
        while (i < length(A))

            if A(2,i) < 0 && A(2,i+1) > 0 % look for sign change in differences
                y_down(A(1,i)/Ndown:A(1,i+1)/Ndown) = ...
                    mean(y_down(A(1,i)/Ndown:A(1,i+1)/Ndown)); 
                % replace all values in between with mean
            end

            i=i+1;
        end
        k=k+1;
    end
    
%% Detect NPS pulses, build pulse matrix P
    i=1;
    k = 0;
    while (i < size(A,2))
        if A(2,i) < 0 && A(2,i+1) > 0 && A(3,i) > 0 % check non-repeats
            k=k+1;
            P(k,1) = A(1,i);  % Start index
            P(k,2) = y_detrend(A(3,i)); % normalized baseline current
            P(k,3) = mean(y_detrend(A(1,i)+1:A(1,i+1)-1)); % avg current between pulses
            P(k,4) = A(1,i+1)-A(1,i); % change in index
            P(k,5) = y_smoothed(A(3,i)); % baseline current, no detrend
            P(k,6) = A(1,i+1); % End index
            %P(k,7) = y_detrend(A(1,i+1)+backset); %normalized baseline current after pulse
        end
        
            i=i+1;
    end
    
    clear backset; 
    
    if exist('P')==0
            OUT=zeros(8,1);
            return
    end
%% 
    N=1;
    out(:,N)=[0; 0; 0; 0; 0; 0; 0; 0];
    for k=1:size(P,1)-5
        start_index = P(k,1); % starting index of pulse
        fileN = floor(start_index/20000)+file_num(1); % compute which file # the pulse is in
        I = P(k,5); % compute baseline current
        dI = (abs(P(k,2) - P(k,3))+abs(P(k,2) - P(k+1,3)))/2; % average NP current drop
        dI1 = P(k,2) - P(k+2,3); % squeeze current drop 
        dT = (P(k,6) - P(k,1))/Fs; % NP transit time in ms
        dT1 = (P(k+2,6)-P(k+2,1))/Fs; % squeeze transit time in ms
        
       % post-squeeze NP current drops
        dI4 = abs(P(k,2) - P(k+3,3)); 
        dI5 = abs(P(k,2) - P(k+4,3));
        dI6 = abs(P(k,2) - P(k+5,3));

        % recovery time determined when post-squeeze NP current drop reaches
        % pre-squeeze NP current drop
        if abs(dI6-dI)/dI < error % last pulse is recovered
            Tr = (P(k+4,1) - P(k+2,1))/Fs; % recovery time in ms

        elseif abs(dI5-dI)/dI < error % 2nd to last pulse is recovered
            Tr = (P(k+3,1) - P(k+2,1))/Fs; % recovery time in ms

        elseif dI4 < dI*(1-error) && dI5 < dI*(1-error) && dI6 < dI*(1-error) % never recovers
            Tr = inf;

        else
            Tr = 0; 
        end
        
        if dI>0 && dI1>0
        out(:,N)=[fileN; start_index; I; dI; dI1; dT; dT1; Tr];
        N=N+1;
        end   
    end

    OUT=out(:,1);

    if out(1,1)~=0
        hFig=figure(ii);
        set(hFig,'Position', [0 300 1200 400])

        subplot(1,2,1)
        plot(y_down,'r-');
        xlabel('Time')
        ylabel('I')
        set(gca,'FontSize',10)

        subplot(1,2,2)
        plot(y_diff)
        hold on
        %plot(pks(:,1),pks(:,2),'ob');    
        plot(A(1,:)./Ndown,A(2,:),'ro')
        hold off
        title('Original signal (y)')
        ylabel('dI/dt')
        xlabel('Time')
        set(gca,'FontSize',10)
    end
    
    fprintf([ 'ENTER will return data,\n' ... 
                        'S will reture 0\n' ...
                        ]);
    OK = input('---\n','s');
                
    switch OK
        case []
            OUT=out(:,1);
        case{'s','S'}
            OUT=zeros(8,1);
    end
    
 elseif flag==0
         OUT=zeros(8,1);
 end
end