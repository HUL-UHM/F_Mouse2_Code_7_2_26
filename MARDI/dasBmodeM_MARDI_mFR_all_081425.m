clear all
% close all
clc
addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\mouse3';
flagPI=0; % 1 = process PI, 0 = Process Fundamental
xDucer='L115';
if flagPI
    bName='BmodeMat_PI';
else
    bName='BmodeMat_Fun';
end
fSh='HARF';
rSh='Bmode';
trFNumAll=[0.75;1.5];fNAll{1,1}='F0p75';fNAll{2,1}='F1p5';%fNAll{3,1}='F4p0';
LogThreshold=-60;
nPul=2;
nchan=128;
numRays = 128;
ARFlocs=1:nchan;
numpushlocs=length(ARFlocs);
lengthseq = numRays*nPul;
% buffNum=2;
for mm=1:size(main_path_all,1)

    main_path=main_path_all{mm,1};
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;

    file_name_all=[];

for kk=dummy_var:length(temp)

    file_name=subdir(fullfile(main_path,temp(kk).name,['*_',fSh,'*_Ac*.mat']));
    if ~isempty(file_name)
        
        file_name_cell=struct2cell(file_name);
        if kk==dummy_var && done_already==0
             done_already=1;
            file_name_all=file_name_cell(1,1:size(file_name_cell,2));
        else
            file_name_all=[file_name_all,file_name_cell(1,1:size(file_name_cell,2))];
        end
    else
       dummy_var=kk+1;
    end
   
end
%%
for bb=1:size(file_name_all,2)
    disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
    [filePath,fileName]=fileparts(file_name_all{1,bb});
    savePath=fullfile(filePath,bName);
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    trFNum=trFNumAll(2,:);
    fN=fNAll{2,1};
    saveName=strrep(fileName,fSh,[rSh,fN]);
    if ~exist(fullfile(savePath,[saveName,'.mat']),'file') && contains(filePath,'MARDI')
        clear Receive P aimg B_data
        load(file_name_all{1,bb});  
    
        if ~exist('Receive','var')
            cd(filePath)
            fL=dir('RcvBf_*.mat');
            fnNew=fL(size(fL,1)).name;
            load(fullfile(filePath,fnNew),'P','Receive','aimg');
            if ~exist('Receive','var')
                fnNew=fL(2).name;
                 load(fullfile(filePath,fnNew),'P','Receive','aimg');
            end
        end
        for uu=1:size(trFNumAll,1)
            trFNum=trFNumAll(uu,:);
            fN=fNAll{uu,1};
            
            axLen=length(Receive(1).startSample:Receive(1).endSample);tic
            ind=find( cellfun(@(x)isequal(x,1),{Receive.bufnum}),1,'first'); 
            frameTotal=size(B_data,3);% use middle frames
            RF_bmode_all=cell(frameTotal,nPul);
        %         fs=aimg.samplingRateMHz*1e6;
            fc =aimg.frequencyMHz*1e6;
            for frameNum=1:frameTotal
                temp=B_data(Receive(ind).startSample:Receive(ind+lengthseq-1).endSample,1:nchan,frameNum);
                temp=reshape(temp,[axLen,nPul,numRays,nchan]);
                temp=double(permute(temp,[1 4 3 2]));
                
                aimg.trackFNum=trFNum;
                if contains(filePath,'\PI')
                    aimg.bpf=0;aimg.bpfVal=[0 0];
                elseif contains(filePath,'\Fun')
                    aimg.bpf=1;aimg.bpfVal=[4 12.5];
                end
                aimg.apGrowth=1;
                aimg.track_time=1; % in ms
                aimg.figuresOn=0;
                aimg.iswiper=0;
                aimg.ARF_loc=1:128;
                aimg.StartDepthMm=P.startDepthMm; 
                for ii=1:nPul
                    aimg.data=squeeze(temp(:,:,:,ii));
                    aimg.scatSum = DAS_HARF_0526(aimg,0,0);
                    RF_bmode_all{frameNum,ii}=squeeze(aimg.scatSum);
                end
            end
            frameSelect=1;
            blat1=aimg.XMTspacingMM*((1:size(RF_bmode_all{frameSelect,1},2))-mean(1:size(RF_bmode_all{frameSelect,1},2)))'; %[mm]
            fs=  aimg.samplingRateMHz;
            c=aimg.c;  
            baxial=P.startDepthMm+(1000*c/(2*fs*1e6))*(0:size(RF_bmode_all{frameSelect,1},1)-1);
            shiftVal=1;
            baxial1=baxial(shiftVal:end);
            RF_bmodeAll=nan([size(RF_bmode_all{1,1}),nPul,frameTotal]);
            for frameSelect=1:frameTotal
                for ii=1:nPul
                    RF_bmodeAll(:,:,ii,frameSelect)=RF_bmode_all{frameSelect,ii};
                end
            end
    
            %%
            if contains(filePath,'\PI')
                if flagPI
                    RF_bmode=squeeze(sum(RF_bmodeAll,3));
                else
                    RF_bmode=squeeze(RF_bmodeAll(:,:,1,:)-RF_bmodeAll(:,:,2,:));
                end
            else
                 RF_bmode=squeeze(sum(RF_bmodeAll,3));
            end
            RF_bmode=squeeze(nanmean(RF_bmode,3));
            if contains(xDucer,'L74')
                if contains(filePath,'\PI')
                    if flagPI
                        if fs < 21
                            N = 10; Wn = [6.5 10]/(fs/2);
                        else
                            N = 10; Wn = [6.5 10]/(fs/2);
                        end
                    else
                        N = 10; Wn = [2 8]/(fs/2);
                    end
    %             elseif contains(filePath,'Fun')
    %                 N = 10; Wn = [2 8]/(fs/2);
                end
            elseif contains(xDucer,'L115')
                if contains(filePath,'\PI')
                    if flagPI
                        if fs < 32
                            N = 10; Wn = [7. 12.5]/(fs/2);
                        else
                            N = 10; Wn = [6.5 14.2]/(fs/2);
                        end
                    else
                        N = 10; Wn = [2 12]/(fs/2);
                    end
    %             elseif contains(filePath,'Fun')
    %                 N = 10; Wn = [2 12]/(fs/2);
                end
            elseif contains(xDucer,'L2214')
                N = 10; Wn = [13 24]/(fs/2);
            end
            if contains(filePath,'\PI')
                [b, a] = butter(N, Wn);
                RF_bmode = filtfilt(b, a, RF_bmode);
            end
            blat=blat1(1):0.01:blat1(end);
            baxial=baxial1(1):0.01:baxial1(end);
            [x,y]=meshgrid(blat1,baxial1);
            [X,Y]=meshgrid(blat,baxial);
            RF_bmodeIp=interp2(x,y,RF_bmode,X,Y,'spline');
    
    
            bmodeOrg=abs(hilbert(RF_bmodeIp(shiftVal:end,:)));
            bmode = (abs(hilbert(RF_bmodeIp(shiftVal:end,:))));
            bmodeLog= 20*log10(bmode/max(bmode(:)));
            bmodeLog(bmodeLog<LogThreshold)=LogThreshold;
    
            [len, width]=size(bmode);
            [counts,centers]=hist(bmode(10:end,round(size(bmode,2)/2)),100);
            bmode_cl=bmode(:);
            ind=find(counts(5:length(counts))>0.2*max(counts)); % skip first few bins due to waterpath
            cutVal=centers(min(ind)+9);
            [index, val]=find(bmode_cl<cutVal);
            bmode_cl(index)=cutVal;
            bmode_norm=reshape(bmode_cl,len,width);
            bmodeLogNorm= 20*log10(bmode_norm/max(bmode_norm(:)));
            bmodeLogNorm(bmodeLogNorm<LogThreshold)=LogThreshold;
    %%
            bmodeHist = adapthisteq(bmode_norm/max(bmode_norm(:)));
            % bmodeHist = adapthisteq(bmode/max(bmode(:)));
            bmodeLogHist= 20*log10(bmodeHist/max(bmodeHist(:)));
            bmodeLogHist(bmodeLogHist<LogThreshold)=LogThreshold;
            savefast(fullfile(savePath,saveName),...
                    'blat','baxial','bmode','bmodeOrg','bmode_norm','bmodeHist','RF_bmode','RF_bmodeAll','RF_bmodeIp','blat1','baxial1')
        end
    
            
    %     end
    
        
    end
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
end
end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
%% plot frequence of RF  data
% L=5*size(RF_bmodeAll,1);
% latInd=1:10; indT=round(mean(latInd))
% YfAl_l=zeros(L,numel(latInd));
% YfAl_2=zeros(L,numel(latInd));
% for nn=1:numel(latInd)
%     dd1=squeeze(RF_bmodeAll(:,latInd(nn),1,1));
%     YfAl_l(:,nn)=abs(fft(dd1,L)/L);
%     dd1=squeeze(RF_bmodeAll(:,latInd(nn),1,1))+squeeze(RF_bmodeAll(:,latInd(nn),2,1));
%     YfAl_2(:,nn)=abs(fft(dd1,L)/L);
% end        
% fs=  aimg.samplingRateMHz;
% freqVec=fs*(0:(L/2))/L;
% Yf=movmean(mean(YfAl_l,2),300);
% YfT=Yf(1:L/2+1);
% YfT(2:end-1)=2*YfT(2:end-1);
% YfLog_1=20*log10(YfT/max(YfT));
% figure;plot(baxial1,squeeze(RF_bmodeAll(:,indT,1,1)),'lineWidth',2)
% hold all;plot(baxial1,squeeze(RF_bmodeAll(:,indT,2,1)),'lineWidth',2)
% leg=legend('+ pulse','- pulse');
% set(leg,'location','best')
% figure;plot(freqVec,YfLog_1,'r');
% hold all;
% Yf=movmean(mean(YfAl_2,2),300);
% YfT=Yf(1:L/2+1);
% YfT(2:end-1)=2*YfT(2:end-1);
% YfLog_2=20*log10(YfT/max(YfT));
% plot(freqVec,YfLog_2,'b');
% plot(freqVec,-6*ones(size(freqVec)),'k');
% leg=legend('Fund','Harm');
% set(leg,'location','best')
% % % figure;imagesc(blat,baxial,bmodeLogNorm,[-55 0]);colormap gray;ylim([5 35])
%%
LogThreshold=-60;
bmodeLogNorm= 20*log10(bmode/max(bmode(:)));
bmodeLogNorm(bmodeLogNorm<LogThreshold)=LogThreshold;
figure;imagesc(blat,baxial,bmodeLogNorm,[LogThreshold 0]);colormap gray;ylim([15 35]);xlim([-15 15])