clear all
% close all
clc
addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\mouse3';
xDucer='L115';
flagPI=0; % 1 = process PI, 0 = Process Fundamental
if flagPI
    bName='Process_PI';
    rmFr=0;
else
    bName='Process_Fun';
    rmFr=0;
end
coreNum=6;
trFn=1.5;trFnT='1p5';
for mm=1:size(main_path_all,1)

        buffNum=2;

    main_path=main_path_all{mm,1}
cd(main_path)
temp=dir;
temp(~[temp.isdir]) = [];  %remove non-directories
dummy_var=3;
done_already=0;

file_name_all=[];

for kk=dummy_var:length(temp)

    file_name=subdir(fullfile(main_path,temp(kk).name,'RcvBf_*.mat'));
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
    knLen =4;

    saveName=strrep(fileName,'RcvBf',['bf',trFnT,'rF',num2str(rmFr)]);

    if ~exist(fullfile(savePath,[saveName,'.mat']),'file') && contains(filePath,'MARDI')
        clear Receive P aimg
        load(file_name_all{1,bb});  
        dP=pwd;
        if ~exist('Receive','var')
            if  contains(filePath,'_DoA') || contains(filePath,'mm')
                cd(filePath)
                cd ..
                cd ..
                tempDir=dir;
                tempDir(~[tempDir.isdir]) = [];  %remove non-directories
                cd(tempDir(end).name)
                tempDir=dir;
                indAcq=[];
                for qq=3:size(tempDir,1)
                    indTemp=strfind(tempDir(qq).name,'_A');
                    indAcq(qq)=str2double(tempDir(qq).name(4:indTemp-1));
                end
                [~,lastAcq]=max(indAcq);
                cd(tempDir(lastAcq).name)
                 fL=dir('RcvBf_*.mat');
                load(fL(1).name,'P','Receive','aimg');

            else
                cd(filePath)
                fL=dir('RcvBf_*.mat');
                fnNew=fL(size(fL,1)).name;
                load(fullfile(filePath,fnNew),'P','Receive','aimg');
                if ~exist('Receive','var')
                    fnNew=fL(2).name;
                     load(fullfile(filePath,fnNew),'P','Receive','aimg');
                end
                

            end
        end
        cd(dP)


    aimg.maxDisplacement =60e-6; 
    tic
    aimg.Receive=Receive;
    aimg.figuresOn = 0;
%     aimg.ImgDataP = ImgDataP;
    % clearvars -EXCEPT aimg RcvData Receive Event
    if isfield(aimg,'nchan')
        nchan=aimg.nchan;
    else
        nchan=128;
    end
%     aimg.frequencyMHz=6.097; % tracking frequency 
%     Receive=aimg.Receive;
    lengthseq=aimg.numTrackLines;
    ARFlocs=aimg.ARF_loc;
    numpushlocs=length(ARFlocs);
    numRays = aimg.numRays;
    aimg.StartDepthMm=P.startDepthMm;
%     numpar = aimg.numParLines;
    aimg.rmFrame=rmFr;
    aimg.nAng =1;
    aimgIni=aimg;
            
    % Reshape the RcvData array into one that has dimensions of [axial x
    % lateral x rayline*acquisition] including the parallel receives, then
    % circshift the data to move each rayline to the center to prepare the data
    % for DnS.
    %%


%%  rearrange the channel data

    axLen=length(Receive(1).startSample:Receive(1).endSample);
    tic
    temp=newRData(:,1:nchan,:);
    try 
        temp=permute(reshape(temp,[axLen,length(aimg.track_time),nchan,length(ARFlocs)]),[1 3 4 2]);
    catch 
       warning('dimension is wrong'); 
       temp=permute(reshape(temp,[axLen,aimg.numTrackLines,nchan,length(ARFlocs)]),[1 3 4 2]);   
       temp(:,:,:,1)=temp(:,:,:,3);
       temp(:,:,:,2)=temp(:,:,:,4);
       temp(:,:,:,end)=[];
    end
    aimg.data=single(temp);
    clear temp
    if contains(filePath,'\PI')
        aimg.bpf=0;aimg.bpfVal=[3.0 8.5];
    elseif contains(filePath,'\Fun')
        aimg.bpf=1;aimg.bpfVal=[5.5 11];
    end
    aimg.trackFNum=trFn;
    % beam form
    scatSum = DAS_HARF_0526(aimg,0,coreNum);
    aimg=rmfield(aimg,'data');
    PRF_interval=median(diff(aimg.track_time));
    diffT=round(diff(aimg.track_time)*100)/100;
    if aimg.rmFrame==1
        indRm=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
    elseif aimg.rmFrame==2
        indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
        indRm_2=circshift(indRm_1,1); % 2nd frame just after push
        indRm=logical(indRm_1.*indRm_2);
    elseif aimg.rmFrame==3
        indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
        indRm_2=circshift(indRm_1,1); % 2nd frame just after push
        indRm_3=circshift(indRm_2,1); % 2nd frame just after push
        indRm=logical(indRm_1.*indRm_2.*indRm_3);  
    elseif aimg.rmFrame==4
        indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
        indRm_2=circshift(indRm_1,1); % 2nd frame just after push
        indRm_3=circshift(indRm_2,1); % 2nd frame just after push
        indRm_4=circshift(indRm_3,1); % 2nd frame just after push
        indRm=logical(indRm_1.*indRm_2.*indRm_3.*indRm_4);
    end
    if aimg.rmFrame>0
        aimg.track_time=aimg.track_time(indRm);
        aimgIni.track_time=aimgIni.track_time(indRm);
        scatSum=scatSum(:,:,indRm);
    end

    if contains(filePath,'\PI')
        if flagPI
            unFiltSum=scatSum+circshift(scatSum,1,3);
        else
            unFiltSum=scatSum-circshift(scatSum,1,3);
        end
        unFiltSum(:,:,1)=unFiltSum(:,:,2);
    elseif contains(filePath,'\Fun')
        aimg.scatSum=scatSum;
    end
    aimgIni.track_timeCom=aimgIni.track_time;
    aimg.track_timeCom=aimg.track_time;
    fs=  aimg.samplingRateMHz;
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
%         elseif contains(filePath,'Fun')
%             N = 10; Wn = [2 8]/(fs/2);
        end
    elseif contains(xDucer,'L115')
        if contains(filePath,'\PI')
            if flagPI
                if fs < 32
                    N = 10; Wn = [7.0 12.5]/(fs/2);
                else
                    N = 10; Wn = [7 14.2]/(fs/2);
                end
            else
                N = 10; Wn = [2 12]/(fs/2);
            end
%         elseif contains(filePath,'Fun')
%             N = 10; Wn = [2 12]/(fs/2);
        end
    elseif contains(xDucer,'L2214')
        N = 10; Wn = [13 24]/(fs/2);
    end

    
    if contains(filePath,'\PI')
        [b, a] = butter(N, Wn);
        aimg.scatSum= filtfilt(b,a,unFiltSum);                                                                         
    end

    savefast(fullfile(savePath,saveName),'aimg','aimgIni')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                close all
    toc
    end

end
    delete(gcp('nocreate'))
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
run2M_dispOnly_MARDI_mFR_all_081425
run3M_p2pdOnly_MARDI_mFR_all_081425
run4M_slope_MARDI_mFR_all_081425
dasBmodeM_MARDI_mFR_all_081425
%% 
% arfidata=dispData{1,1};
% latLoc=[-5 0 5];
% for ii=1:numel(latLoc)
%     figure(ii);hold all
%      plot(t(1:end),squeeze(arfidata(knnsearch(axial,32),knnsearch(lat,latLoc(ii)),:)))
% %      plot(t(1:end-1),squeeze(data1(knnsearch(axial*1e3,32),knnsearch(lat*1e3,latLoc(ii)),:)))
% %     plot(squeeze(arfidata(knnsearch(axial,25),knnsearch(lat,latLoc(ii)),:)),'r')
% end
%%
% LogThreshold=-60;
% bmodeLogNorm= 20*log10(bmode/max(bmode(:)));
% bmodeLogNorm(bmodeLogNorm<LogThreshold)=LogThreshold;
% figure;imagesc(blat,baxial,bmodeLogNorm,[-60 0]);colormap gray;ylim([6 20]);xlim([-8 8])