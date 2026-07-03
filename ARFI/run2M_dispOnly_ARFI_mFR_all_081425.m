%clear all
% close all
clc

addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\Mouse2\2026.03.03';
axL=[8 26];
flagDisp=3; % 1 = Loup, 2 = Loup inter Frame, 3= Xcorr, 4 = BSS, 5 = BSSiF
fSh='bf1p5rF0';
rmF=1;
if flagDisp==1
    rSh='LouL';
    flagCum=1;
elseif flagDisp ==2
    rSh='LouIfL';
    lag=2;
    flagCum=0;
elseif flagDisp ==3
    rSh='xCorL';
    flagCum=1;
elseif flagDisp ==4
    rSh='Bss';
    flagCum=1;
    latKern=1;
elseif flagDisp ==5
    rSh='BssIf';
    lag=2;
    flagCum=0;
    latKern=2;
end
knLen =4;
interpFac=4;
coreNumDisp=4;
for mm=1:size(main_path_all,1)
 % if mm==1
    
 % else
     % axL=[20 45];
 % end
   main_path=main_path_all{mm,1}
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;
    
    file_name_all=[];

for kk=dummy_var:length(temp)

    file_name=subdir(fullfile(main_path,temp(kk).name,[fSh,'*.mat']));
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
% parpool(coreNumDisp)
for bb=1:size(file_name_all,2)

    disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
    [filePath,fileName]=fileparts(file_name_all{1,bb});
    
    if contains(rSh,'LouIf') 
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),fSh(end-5:end),'Lg',num2str(lag),'rF',num2str(rmF),]);
    elseif contains(rSh,'Bss') && flagCum==1
         saveName=strrep(fileName,fSh,[rSh,num2str(knLen),num2str(latKern),fSh(end-5:end),'rF',num2str(rmF)]);
    elseif contains(rSh,'BssIf') && flagCum~=1
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),num2str(latKern),fSh(end-5:end),'Lg',num2str(lag),'rF',num2str(rmF),]);
    else
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),fSh(end-5:end),'rmF',num2str(rmF)]);
    end
    if ~exist(fullfile(filePath,[saveName,'.mat']),'file') && contains(filePath,'ARFI') && contains(filePath,'mm')
         tic
        tempData=load(file_name_all{1,bb});
        aimg=tempData.aimg;
        aimgIni=tempData.aimgIni;      
        aimg.kernelLengthWavs = knLen;       
        aimg.maxDisplacement =80e-6; 
        aimg.axL=axL;       
        aimgIni.axL=axL;
        aimg.rmFrRF=rmF;
        aimgIni.rmFrRF=rmF;
        PRF_interval=median(diff(aimg.track_time));

        % needed if we want to use the matlab version of the code
        ConfigParam.PRE_INTERP_FACTOR = 4;
        ConfigParam.field_sample_freq = aimg.samplingRateMHz*1e6;
        ConfigParam.c = aimg.c;
        ConfigParam.fo = aimg.frequencyMHz*1e6;
        ConfigParam.WAVELENGTHS = aimg.kernelLengthWavs;
        ConfigParam.MAX_SEARCH = aimg.maxDisplacement;
        ConfigParam.axDownFactor = ConfigParam.PRE_INTERP_FACTOR ;
        ConfigParam.latSpace = median(diff(sort(aimg.ARF_loc)))*aimg.elementSpacingMM;
        ConfigParam.CoreNum=coreNumDisp;
        ConfigParam.nAng=aimg.nAng;
        ConfigParam.iswiper=aimg.iswiper;   
        ConfigParam.axL=axL;
        ConfigParam.XMTspacingMM=aimg.XMTspacingMM;
        ConfigParam.numElementsTotal=aimg.numElementsTotal;
        ConfigParam.ARF_loc=aimg.ARF_loc;
        ConfigParam.StartDepthMm=aimg.StartDepthMm;

        diffT=round(diff(aimg.track_time)*100)/100;
        if aimg.rmFrRF==1
            indRm=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
        elseif aimg.rmFrRF==2
            indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
            indRm_2=circshift(indRm_1,1); % 2nd frame just after push
            indRm=logical(indRm_1.*indRm_2);
        elseif aimg.rmFrRF==3
            indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
            indRm_2=circshift(indRm_1,1); % 2nd frame just after push
            indRm_3=circshift(indRm_2,1); % 2nd frame just after push
            indRm=logical(indRm_1.*indRm_2.*indRm_3);    
        elseif aimg.rmFrRF==4
            indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
            indRm_2=circshift(indRm_1,1); % 2nd frame just after push
            indRm_3=circshift(indRm_2,1); % 2nd frame just after push
            indRm_4=circshift(indRm_3,1); % 2nd frame just after push
            indRm=logical(indRm_1.*indRm_2.*indRm_3.*indRm_4);   
         elseif aimg.rmFrRF==5
            indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
            indRm_2=circshift(indRm_1,1); % 2nd frame just after push
            indRm_3=circshift(indRm_2,1); % 2nd frame just after push
            indRm_4=circshift(indRm_3,1); % 2nd frame just after push
            indRm_5=circshift(indRm_4,1); % 2nd frame just after push
            indRm=logical(indRm_1.*indRm_2.*indRm_3.*indRm_4.*indRm_5);  
        else
            indRm=true(size(aimg.scatSum,3),1);
        end
        if aimg.rmFrRF>0
                aimg.track_time=aimg.track_time(indRm);
                aimg.scatSum=aimg.scatSum(:,:,indRm);
                if isfield(aimg,'track_timeCom')
                    aimg.track_timeCom=aimg.track_timeCom(indRm);
                end
                aimgIni.rmFrDisp=0;
        else
            if contains(filePath,'\PI')
                aimgIni.rmFrDisp=0;
            else
                aimgIni.rmFrDisp=1;
            end
        end

        if  flagCum==1
            
            aimg.interpFac=interpFac;
            aimg.parLoopNum=coreNumDisp;
            aimg.flagCum=flagCum;
            if flagDisp==1 
                [arfidata, RF, axial, lat, t] = ncorrWrapper_loupass(aimg);
                clc
                [dispFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial);
                  dispFilt=motion_filter(dispFilt,tNew,tNew(end));
                Peak=max(dispFilt,[],3);
                PeakV=max(diff(dispFilt,1,3),[],3);
                parSave(fullfile(filePath,saveName),'dataFilt',dispFilt,'tNew',tNew,...
                    'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                    'aimgIni',aimgIni,'Peak',Peak,'PeakV',PeakV)
            elseif flagDisp==3
%                 [arfidata, CCcoeff, RF, axial, lat, t] = ncorrWrapper_HARF_xCorr(aimg);
                if isfield(aimg,'track_timeCom')
                    ConfigParam.track_timeCom=aimg.track_timeCom;
                end
                if isfield(aimg,'track_time')
                
                    ConfigParam.track_time=aimg.track_time;
                end
                 tic
                 [arfidata,CCcoeff,axial,lat,t] = dispCalc_xCorr_parfor(aimg.scatSum,ConfigParam);
                 toc
                clc
                [dispFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial);
                 dispFilt=motion_filter(dispFilt,tNew,tNew(end));
                Peak=max(dispFilt,[],3);
                PeakV=max(diff(dispFilt,1,3),[],3);
                parSave(fullfile(filePath,saveName),'dataFilt',dispFilt,'tNew',tNew,...
                    'arfidata',arfidata,'axial',axial,'lat',lat,'t',t,...
                    'aimgIni',aimgIni,'Peak',Peak,'CCcoeff',CCcoeff,'PeakV',PeakV)
            elseif  flagDisp==4
                aimg.latKern=latKern;
                [arfidata, RF, axial, lat, t] = ncorrWrapper_BSS(aimg);
                [dispFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial);
                dispFilt=motion_filter(dispFilt,tNew,tNew(end));
                Peak=max(dispFilt,[],3);
                 PeakV=max(diff(dispFilt,1,3),[],3);
                parSave(fullfile(filePath,saveName),'dataFilt',dispFilt,'tNew',tNew,...
                    'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                    'aimgIni',aimgIni,'Peak',Peak,'PeakV',PeakV)
            end

        elseif  flagCum==0
            aimg.flagCum=flagCum;
            scatSum=aimg.scatSum; 
            aimg.interpFac=interpFac;
            aimg.parLoopNum=coreNumDisp;
            aimg.Lag=lag;
            if flagDisp==2
                aimg.lag=1;
                [arfidata, RF, axial, lat, t] = ncorrWrapper_loupass_iF(aimg);
                arfidata=arfidata(:,:,(lag+1):end);
                t=t((lag+1):end);
                aimgIni.numPreFrames= aimgIni.numPreFrames-lag;
                [dispFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial);
                PeakV=max(dispFilt,[],3);
                dispFilt1=cumsum(dispFilt,3);
                dispFilt1=motion_filter(dispFilt1,tNew,tNew(end));
                Peak=max(dispFilt1,[],3);
                parSave(fullfile(filePath,saveName),'dataFilt',dispFilt,'tNew',tNew,...
                    'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                    'aimgIni',aimgIni,'Peak',Peak,'PeakV',PeakV)
            elseif flagDisp==5
                aimg.latKern=latKern;
                aimg.lag=1;
                [arfidata, RF, axial, lat, t] = ncorrWrapper_BSS(aimg);
                arfidata=arfidata(:,:,(lag+1):end);
                t=t((lag+1):end);
                aimgIni.numPreFrames= aimgIni.numPreFrames-lag;
                [dispFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial);
                PeakV=max(dispFilt,[],3);
                dispFilt1=cumsum(dispFilt,3);
                dispFilt1=motion_filter(dispFilt1,tNew,tNew(end));
                Peak=max(dispFilt1,[],3);
                parSave(fullfile(filePath,saveName),'dataFilt',dispFilt,'tNew',tNew,...
                    'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                    'aimgIni',aimgIni,'Peak',Peak,'PeakV',PeakV)
            end
        end
        toc
        close all
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                close all
    
    end

end

 delete(gcp('nocreate'))
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
end

clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')

%% 
% arfidata=dispData{1,1};
% latLoc=[-5 0 5];
% for ii=1:numel(latLoc)
%     figure(ii);hold all
%      plot(t(1:end),squeeze(arfidata(knnsearch(axial,32),knnsearch(lat,latLoc(ii)),:)))
% %      plot(t(1:end-1),squeeze(data1(knnsearch(axial*1e3,32),knnsearch(lat*1e3,latLoc(ii)),:)))
% %     plot(squeeze(arfidata(knnsearch(axial,25),knnsearch(lat,latLoc(ii)),:)),'r')
% end
