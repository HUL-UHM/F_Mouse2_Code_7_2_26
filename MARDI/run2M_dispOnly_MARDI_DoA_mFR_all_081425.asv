clear all
% close all
clc
addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\Mouse2';

axL=[8 26];
flagDisp=2; % 1 = Loup, 2 = Loup inter Frame, 3= Xcorr, 4 = BSS, 5 = BSSiF
fSh='bf1p5rF0';
if flagDisp==1
    rSh='LouL';
    flagCum=1;
    rmF=0;
elseif flagDisp ==2
    rSh='LouIfL';
    lag=1;
    rmF=1;
    flagCum=0;
elseif flagDisp ==3
    rSh='xCorL';
    flagCum=1;
    rmF=0;
elseif flagDisp ==4
    rSh='Bss';
    flagCum=1;
    latKern=1;
    rmF=0;
elseif flagDisp ==5
    rSh='BssIf';
    lag=1;
    rmF=0;
    flagCum=0;
    latKern=1;
end
knLen =4;
interpFac=4;
coreNumDisp=2;
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
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),fSh(end-5:end),'Lg',num2str(lag),'rF',num2str(rmF)]);
    elseif contains(rSh,'Bss') && flagCum==1
         saveName=strrep(fileName,fSh,[rSh,num2str(knLen),num2str(latKern),fSh(end-5:end),'rF',num2str(rmF)]);
    elseif contains(rSh,'BssIf') && flagCum~=1
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),num2str(latKern),fSh(end-5:end),'Lg',num2str(lag),'rF',num2str(rmF),]);
    else
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen),fSh(end-5:end),'rF',num2str(rmF)]);
    end
    if ~exist(fullfile(filePath,[saveName,'.mat']),'file') && contains(filePath,'MARDI') && contains(filePath,'mm')
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
                    aimg.track_timeCom;
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
                parSave(fullfile(filePath,saveName),'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                'aimgIni',aimgIni)
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
     
                parSave(fullfile(filePath,saveName),'arfidata',arfidata,'axial',axial,'lat',lat,'t',t,...
                'aimgIni',aimgIni,'CCcoeff',CCcoeff)
            elseif  flagDisp==4
                aimg.latKern=latKern;
                [arfidata, RF, axial, lat, t] = ncorrWrapper_BSS(aimg);
                parSave(fullfile(filePath,saveName),'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                'aimgIni',aimgIni)
            end

        elseif  flagCum==0
            aimg.flagCum=flagCum;
            scatSum=aimg.scatSum;
            
            aimg.interpFac=interpFac;
            aimg.parLoopNum=coreNumDisp;
            aimg.lag=lag;
            if flagDisp==2
                [arfidata, RF, axial, lat, t] = ncorrWrapper_loupass_iF(aimg);
                arfidata=arfidata(:,:,(lag+1):end);
                t=t((lag+1):end);
                aimgIni.numPreFrames= aimgIni.numPreFrames-lag;
         
                parSave(fullfile(filePath,saveName),'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                'aimgIni',aimgIni)
            elseif flagDisp==5
                aimg.latKern=latKern;
                [arfidata, RF, axial, lat, t] = ncorrWrapper_BSS(aimg);
                arfidata=arfidata(:,:,(lag+1):end);
                t=t((lag+1):end);
                aimgIni.numPreFrames= aimgIni.numPreFrames-lag;
         
                parSave(fullfile(filePath,saveName),'arfidata',arfidata,'RF',RF,'axial',axial,'lat',lat,'t',t,...
                'aimgIni',aimgIni)
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
% run3M_p2pdOnly_MARDI_mFR_all_081425
% run4M_slope_MARDI_mFR_all_081425

