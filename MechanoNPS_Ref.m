clear all;
close all;
%
Nstart =210;
Nfiles =Nstart+5;
th=1e-4;
th2=1e-3;
error=0.08;
dt=100;
N=20;       %down sampling period

%load DATA--------------------------------------

mydata = cell(1,Nfiles);
for k = Nstart:Nfiles
  myfilename    = sprintf('trial1_%05d.txt', k);
  mydata{k}     = textread(myfilename);
end
  
y=cat(1,mydata{:});
y=fastsmooth(y,100,1,1);
y_detrend=detrend(y);
y_down=downsample(y_detrend,N); 
ym=transpose(y_down);

%parameters
t = 0:1:size(ym,2)-1; 

for i=1:size(ym,2)-1
    slp(i)=ym(i+1)-ym(i);
    if abs(slp(i))<th
       slp(i)=0;
    end
end

for i=1:size(ym,2)-1
    if slp(i)<-th2
       k=i+1;
       for j=k:size(ym,2)-1
            if  abs(slp(j))<th2
                slp(j)=0;
            elseif slp(j)>th2
                break;
            end
       end
    end
end

tmp=0;
for i=1:size(slp,2)
    if abs(slp(i))>0
        tmp=tmp+1;
        A(1,tmp)=i*N;
        A(2,tmp)=slp(i);
        A(3,tmp)=(i-dt)*N;
    end
end

%remove error in slp
k=1;
while k<=4;
tmpr=0;
i=1;
while i<=size(A,2)-1
    if A(2,i)>0 && A(2,i+1)>0 && A(2,i) < A(2,i+1)
       A(1,i)=A(1,i+1);
       A(2,i)=A(2,i+1);
       A(3,i)=A(3,i+1);
    elseif A(2,i)>0 && A(2,i+1)>0 && A(2,i) > A(2,i+1)
           A(1,i+1)=A(1,i);
           A(2,i+1)=A(2,i);
           A(3,i+1)=A(3,i);
    elseif A(2,i)<0 && A(2,i+1)<0 && A(2,i) > A(2,i+1)
           A(1,i)=A(1,i+1);
           A(2,i)=A(2,i+1);
           A(3,i)=A(3,i+1);
     elseif A(2,i)<0 && A(2,i+1)<0 && A(2,i) < A(2,i+1)
            A(1,i+1)=A(1,i);
            A(2,i+1)=A(2,i);
            A(3,i+1)=A(3,i);
    end
    i=i+1;
end
k=k+1;
end

ym2=ym;

%make rectangular signal
ym2=ym;
k=1;
while k<=50;
i=1;
while i<=size(A,2)-1
    if A(2,i)<0 && A(2,i+1)>0
        ym2(A(1,i)/N:A(1,i+1)/N)=mean(ym(A(1,i)/N:A(1,i+1)/N));
    end
    i=i+1;
end
k=k+1;
end

%Figures
figure('units', 'pixels', 'pos',[1200 1200 800 1000])
subplot(3,1,1)
plot(slp,'b-');
title('Differenciation of I')
set(gca,'FontSize',20)
subplot(3,1,2)
plot(y_detrend,'k-');
title('y_{detrend}')
set(gca,'FontSize',20)
subplot(3,1,3)
plot(ym2);
title('ym2')
set(gca,'FontSize',20)

%detect signal
i=1;
while i<=size(A,2)-1
    if A(2,i)<0 && A(2,i+1)>0 && A(3,i)>0
            tmpr=tmpr+1;
            p(tmpr,1)=A(1,i);  %Starttime
            p(tmpr,2)=y_detrend(A(3,i)); %I
            p(tmpr,3)=(mean(y_detrend(A(1,i)+1:A(1,i+1)-1))); %I'
            p(tmpr,4)=A(1,i+1)-A(1,i); %dT
            p(tmpr,5)=y(A(3,i)); %I w/o detrend
            p(tmpr,6)=A(1,i+1); %End time
            p(tmpr,7)=y_detrend(A(1,i+1)+200); %I after squeezing
    end
    i=i+1;
end


N=1;
for k=1:length(p)-6,
        sT=p(k,1);
        fileN=floor(sT/20000)+Nstart;
        I=p(k,5);
        dI=abs(p(k,2)-(p(k,3)+p(k+1,3)+p(k+2,3))/3);
        
        pp(N,:)=[N sT fileN I dI];
        N=N+1;       
end



