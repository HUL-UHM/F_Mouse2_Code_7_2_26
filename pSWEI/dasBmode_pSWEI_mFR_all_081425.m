clear all
% close all
clc
addpath(genpath('D:\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='D:\Zac\2025.12.1';

flagParfor=10;
fSh='pRF';
rSh='Bmode';
trFNumAll=[1.0;3.0];fNAll{1,1}='F1p0';fNAll{2,1}='F3p0';%fNAll{3,1}='F6p0';
LogThreshold=-60;
nChan=128;
numRays = 128;
ARFlocs=1:nChan;
numpushlocs=length(ARFlocs);
lengthseq = numRays;

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
% parpool(flagParfor)
for bb=1:size(file_name_all,2)
    for uu=1:size(trFNumAll,1)
        trFNum=trFNumAll(uu,:);
        fN=fNAll{uu,1};
        disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
        [filePath,fileName]=fileparts(file_name_all{1,bb});
        savePath=fullfile(filePath,'BmodeMat');
        if ~exist(savePath,'dir')
            mkdir(savePath)
        end
        saveName=strrep(fileName,fSh,[rSh,fN]);
        if ~exist(fullfile(savePath,[saveName,'.mat']),'file') && ~contains(filePath,'f100_n10_18_nC5_PRF15')
            clear Receive 
            load(file_name_all{1,bb})

            tempPath=filePath;

            indTemp=strfind(fileName,'_Pos');
            tempName=strrep(fileName,['_Pos',fileName(indTemp+4)],'_Pos2');
            tempName=strrep(tempName,tempName(1:2),'55');
            if contains(tempName,'_Ac2')
                 tempName=strrep(tempName,'Ac2','Ac1');
            end
            if ~exist('Receive','var')
                load(fullfile(tempPath,tempName),'Receive')
            end

            aimgB=aimg;
            aimgB.trackFNum=trFNum;
            aimgB.bpf=1;aimgB.bpfVal=[0 0];   
            aimgB.apGrowth=1;

            aimgB.figuresOn=0;
            aimgB.iswiper=0;
            aimgB.ARF_loc=1:128;
    %         aimgB.StartDepthMm=P.startDepthMm;   
            numpushlocs=length(ARFlocs);
            numRays = aimg.numRays;

            nAng=aimg.nAngB;
            aimgB.track_time=1:(nAng*aimg.numFramesBmode);
            lengthseq=numel(aimgB.track_time);
            axLen=length(Receive(1).startSample:Receive(1).endSample);
            nChan=128;

            tic

            ind=find( cellfun(@(x)isequal(x,1),{Receive.bufnum}),1,'first'); 
            tic
            temp=B_data(Receive(ind).startSample:Receive(ind+lengthseq-1).endSample,1:nChan,:);
            temp=permute(reshape(temp,[axLen,nAng,nChan,aimg.numFramesBmode,]),[1 3 2 4]);
            aimgB.data=single(temp);
            aimgB.nAng=aimg.nAngB;
            aimgB.angles=aimg.anglesB;
            RF_bmodeAll=DAS_SWEI_1202(aimgB,flagParfor); 

            fs=aimgB.samplingRateMHz;
            blat1=aimgB.XMTspacingMM*((1:size(RF_bmodeAll,2))-mean(1:size(RF_bmodeAll,2)))'; %[mm]
            c=aimgB.c;  
            baxial=aimgB.StartDepthMm+(1000*c/(2*fs*1e6))*(0:size(RF_bmodeAll,1)-1);
            shiftVal=1;
            baxial1=baxial(shiftVal:end);

            toc
            %%
            RF_bmode=squeeze(nanmean(nanmean(RF_bmodeAll,3),4));
%             N = 10; Wn = [3.5 8.5]/(fs/2);
%             N = 10; Wn = [4 12]/(fs/2);
% N = 10; Wn = [13 23]/(fs/2);
%             [b, a] = butter(N, Wn);
%             RF_bmode = filtfilt(b, a, RF_bmode);


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
    end
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
end
%  delete(gcp('nocreate'))
end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
%% plot frequence of RF  data
L=5*size(RF_bmode,1);
latInd=63:65;
YfAll=zeros(L,numel(latInd));
for nn=1:numel(latInd)
    dd1=squeeze(RF_bmode(:,latInd(nn),1));
    YfAll(:,nn)=abs(fft(dd1,L)/L);
end        
fs=  aimg.samplingRateMHz;
freqVec=fs*(0:(L/2))/L;
 Yf=movmean(mean(YfAll,2),300);
YfT=Yf(1:L/2+1);
YfT(2:end-1)=2*YfT(2:end-1);
YfLog=20*log10(YfT/max(YfT));
figure;plot(baxial1,squeeze(RF_bmode(:,64,1)))
figure;plot(freqVec,YfLog)
% figure;imagesc(blat,baxial,bmodeLogNorm,[-55 0]);colormap gray;ylim([5 35])
%%
LogThreshold=-60;
bmodeLogNorm= 20*log10(bmode/max(bmode(:)));
bmodeLogNorm(bmodeLogNorm<LogThreshold)=LogThreshold;
figure;imagesc(blat,baxial,bmodeLogNorm,[-60 0]);colormap gray;ylim([15 45]);xlim([-15 15])