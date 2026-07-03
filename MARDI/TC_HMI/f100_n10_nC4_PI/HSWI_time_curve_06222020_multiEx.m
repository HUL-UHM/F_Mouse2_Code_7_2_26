clear all
close all
savePath1='H:\Harmonic\ST-HMI\L74\f100_n10_nC4_PI';
nPulse=1;
pushFreq=4.0;
pushFreq=125/round(125/pushFreq)*1e6
%% multi excitaiton factor
Dmax=100e-6;
Dmin=30e-6;
funFreq=100;
totalHarm=10;
harmFac=1;
totalFreq=1:harmFac:totalHarm;
freqL=length(totalFreq);

Ncycle=4;flagLoad=1; % 2=equal distant push; 1 = previous location;
                 % else point sel
folName=['f',num2str(funFreq),'_n',num2str(totalHarm),'_nC',num2str(Ncycle),'_PI'];
savePath=fullfile(savePath1,folName);
if ~exist(savePath,'dir')
    mkdir(savePath)
end
AmFac=1.25;
sF=2;
harmPower=4;

filePath='H:\Harmonic\ST-HMI\L74\f100_n10_nC4_PI';
fileName='f1Fac_2_HP_4_3_5_nC_4.mat';
%%
prfTemp=10e3; % 10 kHZ PRF
if rem(prfTemp,funFreq)~=0
PRF=prfTemp + (funFreq - rem(prfTemp,funFreq));
else
   PRF=prfTemp;
end
if flagLoad==2
    samInt=10*1/PRF;
    iniSam=10;
end
PRFcom=round(PRF/nPulse);
PRF=PRFcom;
t=(0:1/PRF:(Ncycle/funFreq));
tAll=repmat(t,freqL,1);
% if harmPower==1
if harmPower~=1
    An=(1:freqL).^harmPower;
else
    An=ones(1,totalHarm);
end
% elseif harmPower==2
%     An=ones(1,totalHarm).^2;
% end

AnAll=repmat(An',1,length(t));
selFreq=totalFreq*funFreq;
freqAll=repmat(selFreq',1,length(t));
thetaAll=zeros(freqL,length(t));
thetaAll(1:2:freqL,:)=pi;
yAll=AnAll.*cos(2*pi*freqAll.*tAll+thetaAll);
yt=sum(yAll);
Am=-AmFac*min(yt);
ytAm=sqrt((Am+yt).^sF);
% ytAm = rescale(yt,Dmin,Dmax,'InputMin',min(yt),'InputMax',max(yt));
%%
figure;plot(t,yAll)
figure;plot(t*1e3,yt/max(yt))
hold all;plot(t*1e3,ytAm/max(ytAm))
%%
L=length(t);

Yf=fft(yt);YfAm=fft(ytAm);
P2 = abs(Yf/L);P2Am = abs(YfAm/L);
P1 = P2(1:round(L/2)+1);P1Am = P2Am(1:round(L/2)+1);
P1(2:end-1) = 2*P1(2:end-1);P1Am(2:end-1) = 2*P1Am(2:end-1);
freqVec=PRF*(0:round(L/2))/L;
h1=figure;plot(freqVec,log10(P1/max(P1)));xlim([0 2e3])
hold all;plot(freqVec,log10(P1Am/max(P1Am)))
%%
ytD=Dmax*ytAm/max(ytAm);
h2=figure;plot(t,ytD,'r');xlim([0 1/funFreq])
YfD=fft(ytD);
P2D = abs(YfD/L);
P1D = P2D(1:round(L/2)+1);
P1D(2:end-1) = 2*P1D(2:end-1);
figure(h1);plot(freqVec,log10(P1D/max(P1D)))
if flagLoad==1
load(fullfile(filePath,fileName),'pushLocT')
elseif flagLoad==2
    pushLocT=(iniSam/PRF):samInt:((1/funFreq)-(iniSam/PRF));
else
pushLocT=getpts(h2);
end
%%
pushLoc=zeros(size(pushLocT));
nPush=length(pushLoc);
ytS=zeros(size(ytD));
t=t';
for ii=1:nPush
    pushLoc(ii)=t(knnsearch(t,pushLocT(ii)));
    ytS(knnsearch(t,pushLocT(ii)):(PRF/funFreq):(Ncycle*PRF/funFreq))=ytD(knnsearch(t,pushLoc(ii)):(PRF/funFreq):(Ncycle*PRF/funFreq));
end
figure;plot(t*1e3,ytD)
hold all
plot(t*1e3,ytS,'ro')

YfS=abs(fft(ytS));
P2S = abs(YfS/L);
P1S = P2S(1:round(L/2)+1);
P1S(2:end-1) = 2*P1S(2:end-1);
figure;plot(freqVec,log10(P1S/max(P1S)));xlim([0 2e3])
title('FFT of ytS')
%% 
ft=PRF*linspace(0,1,L);
ftT=ft(knnsearch(ft',selFreq(1)-50):knnsearch(ft',selFreq(end)+50));
YfST=YfS(knnsearch(ft',selFreq(1)-50):knnsearch(ft',selFreq(end)+50));
[fftAmpFreq,freqLoc]=findpeaks(YfST,ftT,'NPeaks',totalHarm,'MinPeakDistance',funFreq-50);
fftAmpFreq./fftAmpFreq(1)
%%
ind=pushLoc(1)*1e4;

relaxationTime=1e-3;



on_duration=t(end);            
total_time=on_duration+relaxationTime;%[s] % in this time there will be no push
sampling_freq=PRF;
t_total = (0:1/sampling_freq:total_time)';
t_on=(0:1/sampling_freq:on_duration)';
d_final=zeros(size(t_total));
d_final(1:length(ytS))=ytS;
y_final=double(d_final>0);

NFFT=5e5;
indT=knnsearch(t_total,t_on(end));

figure;plot(t_total(1:indT)*1e3,d_final(1:indT),'ro');
xlabel('Sequence Time (ms)')
ylabel ('Force Durtation (\mu s)')
figure;plot(sampling_freq*linspace(0,1,NFFT),abs(fft(d_final(1:indT)-mean(d_final(1:indT)),NFFT)))
xlim([0 2e3])
%%
% ft=sampling_freq*linspace(0,1,NFFT);
% ftT1=ft(knnsearch(ft',selFreq(1)-50):knnsearch(ft',selFreq(end)+50));
% YfS1=abs(fft(d_final(1:indT)-mean(d_final(1:indT)),NFFT));
% YfST1=YfS1(knnsearch(ft',selFreq(1)-50):knnsearch(ft',selFreq(end)+50));
% [fftAmpFreq1,freqLoc1]=findpeaks(YfST1,ftT1,'NPeaks',totalHarm,'MinPeakDistance',funFreq-50)
% figure;plot(ftT1,YfST1)
% hold all;plot(freqLoc1,fftAmpFreq1,'ro')
%%
tnew=t_total(1):1/PRF:t_total(end);
d_final_interp=interp1(t_total,d_final,tnew);
y_final_interp=interp1(t_total,y_final,tnew);
indPush=find(y_final_interp==1);
pushDuration_all=d_final_interp(1,indPush);
if indPush(1)==1 % start with push
    track_btwn_num=diff(indPush)-1;
    if indT-indPush(end)>0
    track_btwn_num=[track_btwn_num indT-indPush(end)]; % track after last push
    end
    indPush_mod=indPush;
else % 1st push t>0.0
    indPush_mod=[1 indPush];
    track_btwn_num=diff(indPush_mod)-1;
    track_btwn_num(1)=track_btwn_num(1)+1;
    if indT-indPush(end)>0
        track_btwn_num=[track_btwn_num indT-indPush(end)]; % track after last push
    end
end
track_btwn_num(logical(mod(track_btwn_num,2)))=track_btwn_num(logical(mod(track_btwn_num,2)))-1
sum(track_btwn_num<3)
            
           
if ind~=0
    Duty_C=100*sum(pushDuration_all(1:nPush)*funFreq);
%     Duty_C_All(ind_prf,ind+1)=100*sum(pushDuration_all(1:nPush)/push_period);
else
    Duty_C=100*sum(pushDuration_all)/(Ncycle/funFreq);
%     Duty_C_All(ind_prf,ind+1)=100*sum(pushDuration_all)/(NCycle*funFreq);
end
push_cycle= round(pushDuration_all'.*pushFreq);
relaxationTimeTrack=round(relaxationTime*PRF);
%%
iniSamN=indPush(1);
save_name=['f',num2str(flagLoad),'Fac_',num2str(sF),'_HP_',num2str(harmPower),'_',num2str(iniSamN),'_',num2str(nPush),'_nC_',num2str(Ncycle)];
note=['push_start= ',num2str(pushLoc(1)*1e3),'; nPush= ',num2str(nPush),';  Cycle =',num2str(Ncycle)];
saveName=fullfile(savePath,save_name);
save(saveName,'note','indPush_mod','indPush','nPush','push_cycle','funFreq',...
        'track_btwn_num','pushFreq','pushDuration_all','PRFcom','PRF','relaxationTimeTrack','total_time',...
        'on_duration','relaxationTime','Duty_C','Ncycle','pushLocT','fftAmpFreq','freqLoc','Dmin','Dmax','nPulse')



