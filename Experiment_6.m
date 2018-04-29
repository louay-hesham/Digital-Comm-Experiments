n = 10; % However many numbers you want.
inbits1 = randi([0, 1], [1, n]); %Generate random bits of zeros and 1's
inbits = [0 1 0 1 1 0 0 1 0 ];
prompt = 'please enter the amplitude= ';
A=input(prompt);
prompt = 'please enter the time period for each bit= ';
Tb =input(prompt);

Rb = 1/Tb; %---Bit rate
Fs = 4*Rb;
N = length(inbits);   %---Bit Length of input bits
tTb = linspace(0,Tb,100); %---interval of bit time period

%%%%%%%% POLAR NON RETURN TO ZERO
NPRZ=[];
     for k = 1:N
            P=ones(1,length(tTb));
           
            NPRZ = [NPRZ ((-1)^(inbits(k) + 1))*A/2*P];
     end
figure(1);
t = linspace(0,N*Tb,length(NPRZ)) %---Time vector for n bits
subplot(6,1,1);
plot(t,NPRZ,'LineWidth',2);grid on;
title(['polar NRZ line coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(NPRZ)-A,max(NPRZ)+A]);


%%%%%%%% POLAR RETURN TO ZERO
PRZ=[];
 for k = 1:N
            
            c = ones(1,length(tTb)/2);
            b = zeros(1,length(tTb)/2);
            p = [c b]
            PRZ = [PRZ ((-1)^(inbits(k)+1))*(A/2)*p];
 end
       
subplot(6,1,2)  
plot(t,PRZ,'LineWidth',2);grid on;
title(['Polar RZ line coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(PRZ)-A,max(PRZ)+A]);

%%% MANCHESTER
Mchstr=[];
   for k = 1:N
            c = ones(1,length(tTb)/2);
            b = -1*ones(1,length(tTb)/2);
            p = [c b];
            Mchstr = [Mchstr ((-1)^(inbits(k)+1))*A/2*p];
   end

%figure(2);        
subplot(6,1,3);plot(t,Mchstr,'LineWidth',2);grid on;
title(['Manchester line coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(Mchstr)-A,max(Mchstr)+A]);



%%% AMI
AMI=[];
count=0
for k = 1:N
    if (inbits(k)==1)
        
     if(count ==0)   
        y =1;
        count =1;
     else
         y =0;
        count =0;
     end
     P=ones(1,length(tTb));
     AMI = [AMI ((-1)^(y + 1))*A/2*P];
    else
     p= zeros(1,length(tTb));   
     AMI = [AMI (p*(A/2))];
    end
    
    
end

%figure(2);        
subplot(6,1,4);plot(t,AMI,'LineWidth',2);grid on;
title(['AMI coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(AMI)-A,max(AMI)+A]);


%%% MULTI LEVEL

MLVL=0
count = -1
count1 = 0
m=zeros(1,length(tTb));
for  k = 1:N
  if(inbits(k)==1)
       if(count1==1&& count==0)%%
           count =count1;
           count1=0;
           p= zeros(1,length(tTb));
           m = p*(A/2)
           MLVL = [MLVL (m)];
       
       
   elseif(count1==0&& count==1)%%
            count =count1;
            count1=-1;
            P=ones(1,length(tTb));
            m=-1*A/2*P;
            MLVL = [MLVL (m)]; 
        
        
        elseif(count1==-1&& count==0)%%
           count =count1;
           count1=0;
           p= zeros(1,length(tTb));  
           m=p*(A/2);
           MLVL = [MLVL (m)];
           
        
        elseif(count1==0&& count==-1)%%
           count =count1;
           count1=1;
           P=ones(1,length(tTb));
           m=A/2*P;
           MLVL = [MLVL (m)]; 
           
        end
   else
      MLVL=[MLVL (m)]; 
    end
   
end
t = linspace(0,N*Tb,length(MLVL)) ;
        
subplot(6,1,5);
plot(t,MLVL,'LineWidth',2);grid on;
title(['MULTIlevel coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(MLVL)-A,max(MLVL)+A]);


%%% NRZ-I
NRZI=[];
count=1;
m=ones(1,length(tTb))*A/2;
for  k = 1:N
    
    if(inbits(k)==0)
        NRZI=[NRZI (m)]; 
        
    else
         if(count ==0)   
        y =1;
        count =1;
     else
         y =0;
        count =0;
     end
     P=ones(1,length(tTb));
     m=((-1)^(y + 1))*A/2*P;
      NRZI = [NRZI (m)] ;
     end
    
end
t = linspace(0,N*Tb,length(NRZI)) ;
subplot(6,1,6);
plot(t,NRZI,'LineWidth',2);grid on;
title(['NRZI coding for ', num2str(inbits)]);
%xlabel('Time in sec');
%ylabel('Amplitude');
axis([0,max(t),min(NRZI)-A,max(NRZI)+A]);

figure(3)
 [pyy,fy]=psd(NPRZ);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'r');
set(plotHandle,'LineWidth',2);

hold on;

[pyy,fy]=psd(PRZ);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'b');
set(plotHandle,'LineWidth',2);
hold on;

[pyy,fy]=psd(AMI);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'m');
set(plotHandle,'LineWidth',2);
hold on;

[pyy,fy]=psd(MLVL);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'g');
set(plotHandle,'LineWidth',2);
hold on;

[pyy,fy]=psd(Mchstr);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'c');
set(plotHandle,'LineWidth',2);
hold on;

[pyy,fy]=psd(NRZI);
plotHandle=plot(fy*Fs/2,10*log10((pyy)),'k');
set(plotHandle,'LineWidth',2);
legend('NPRZ','PRZ','AMI','MLVL','Manchester','NRZI');
title('PSD of Line Codes');
grid on;
hold off;

