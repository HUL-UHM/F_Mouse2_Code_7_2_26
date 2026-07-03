clear all
% close all
clc
addpath(genpath('D:\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='D:\Zac\2025.12.1';
%ind=ind+1;main_path_all{ind,1}='D:\Zac\2025.11.20';
blat= 0.2985*((1:128)-mean(1:128));
factor.core=10; % number of core for parallel processing
%% 20kpa_5_0mmHg_P1
factor.axL=[5 30];
factor.interpolate=0;
factor.dx=0.1e-3;factor.spaDs=1;
factor.NT=4096;
factor.axAvg=4e-3; %  [m] axial average ROI in 
factor.tAvg=1; %[sample] axial average ROI in t
factor.tWind1=0.3;% tuckey window value
factor.tWind2=0.3;% tuckey window value

factor.tInterFac=4;
factor.NT=round(4096*factor.tInterFac/2);

Flag.dirFiltSong=1;
if Flag.dirFiltSong 
    factor.Dir_cutoffsVel=[7 0.2]; % in m/s
    factor.maxFreq=2.0e3;
    factor.Dir_order=2; % bandpass order
    factor.Dir_power=2; % bandpass power
end
% 2D TOF method
factor.axW=4e-3;factor.laW=factor.axW;
factor.windSize=factor.axW;
factor.windSizeCC=round(factor.windSize/factor.dx);
factor.stepSize=round(factor.axW/factor.dx);




nameSh='_F1p5LouL4';
rSh=['gV',nameSh(3:end)];
foldSh='pSWEI';
if contains(nameSh,'If')
    factor.calDiff=0; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.rmReverbFram=0;% number of frame to remove after push

else
    factor.calDiff=1; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.samDiff=4; % difference sample
    factor.rmReverbFram=1;% number of frame to remove after push

end

for mm=1:size(main_path_all,1)
%     if mm<4
%       factor.axL=[15 50];
%     else
%        factor.axL=[4 41];
%     end
    main_path=main_path_all{mm,1}
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;

    file_name_all=[];

    for kk=dummy_var:length(temp)

%         file_name=subdir(fullfile(main_path,temp(kk).name,['*',nameSh,'*F*_*nC*.mat']));
          file_name=subdir(fullfile(main_path,temp(kk).name,['*',nameSh,'*ang*Pos*.mat']));
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

    for bb=1:size(file_name_all,2)
          disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
        [filePath,fileName]=fileparts(file_name_all{1,bb});
        savePath=fullfile(filePath,['gV_RegCC_',num2str(factor.rmReverbFram),'']);
        if ~exist(savePath,'dir')
            mkdir(savePath)
        end
        saveName=strrep(fileName,nameSh,rSh);
        saveName=saveName(1:end);
        load(file_name_all{1,bb},'aimgIni');
        
            
        if  Flag.dirFiltSong
             saveNameT=[saveName,'_dirSong'];
        else
            saveNameT=saveName;
        end

        pushN=10*(aimgIni.swPushLoc)+(1:numel(aimgIni.swPushLoc));
        
        if  contains(filePath,foldSh)
            if numel(aimgIni.swPushLoc)>1
                saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushN(end)),'_']),'RL.mat'];
            else
                if pushN(1)<64
                    saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushN(1)),'_']),'LR.mat'];
                else
                     saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushN(1)),'_']),'RL.mat'];
                end
            end
        end
        if ~exist(fullfile(savePath,saveNameCh),'file') && contains(filePath,foldSh) && contains(filePath,'dispAll_F')
            tic
            load(file_name_all{1,bb},'dispData','lat','axial','t','aimgIni');

            for uu=1:size(dispData,1)%min(size(dispData,1),10)
                if uu==1
                    axialOld=axial;latOld=lat;tOld=t;
                else
                    axial=axialOld;lat=latOld;t=tOld;
                end

                pushL=aimgIni.swPushLoc(uu);


                if pushL<64
                     factor.Dir_angles=0; %directional angle
                     saveNameT1=[saveNameT,'LR'];
                     direc='lr'
                else
                     factor.Dir_angles=pi; %directional angle
                     saveNameT1=[saveNameT,'RL'];
                     direc='rl'
                end
               
                arfidata=dispData{uu,1};
                factor.pushLoc=blat(pushL)*1e-3;
                factor.latCheck=round(pushL*numel(lat)/128);
                aimg=aimgIni;           
                factor.SubPreFrTime=0;
                 aI=knnsearch(axial,factor.axL(1)):knnsearch(axial,factor.axL(2));
                if factor.calDiff
                    [data,axial,lat,t]=rmRevIsoPixel(arfidata(aI,:,:),axial(aI),lat,t,factor,aimgIni);
                else
                    [data,axial,lat,t]=rmRevIsoPixelIF(arfidata(aI,:,:),axial(aI),lat,t,factor,aimgIni);
                end
%                 dataOld=data;
%                 data=dataOld(:,:,knnsearch(t'*1e3,-0):end);
%                 t=t(knnsearch(t'*1e3,-0):end);
%                 dataF=data;
%                 data=motion_filter(data,t,t(end));
                % do filtering to remove some noise
                [b,a]=butter(2,2*[50 2000]*median(diff(t)),'bandpass'); 
                dataF=permute(filtfilt(b,a,permute(data,[3 1 2])),[2 3 1]);

                %  average
                dataF=movmean(dataF,round(factor.axAvg/factor.dx));% axial direction
                dataF=permute(movmean(permute(dataF,[3 1 2]),factor.tAvg),[2 3 1]); % 5 point along time
                if factor.calDiff
%                     dataF=diff(dataF,1,3); % particle velocity
%                     t=t(1:end-1);
                     factor.samDiff=4;
                    dataF1=circshift(dataF, factor.samDiff,3); % particle velocity
                    dataF=dataF-dataF1;
                    dataF(:,:,1:  factor.samDiff)=0;
%     t=t(1:end-1);
                end                    
                % Tukey window
                tWin=tukeywin(size(dataF,3),factor.tWind1);
                tWinMat=permute(repmat(tWin,[1 size(dataF(:,:,1))]),[2 3 1]);
                dataw=dataF.*tWinMat;

                fixedHiFreq=factor.maxFreq;
                factor.Dir_cutoffs=fixedHiFreq*factor.dx*0.5./factor.Dir_cutoffsVel;
                if Flag.dirFiltSong
                    dataw = df2d_Song_V4(dataw, factor.dx, factor.dx, factor.Dir_cutoffs,factor.Dir_order,factor.Dir_power,factor.Dir_angles);
                end
                % 2D TOF method
                [cs.TOF,r2Lval,csx.TOF,csy.TOF,ccx.TOF,ccy.TOF,cc.TOF,PeakL]=SWV_calc_2D_2(abs(dataw),axial,lat,t,factor.axW,factor.laW,factor.tInterFac,factor.windSize);
                
                
                % 2D CC
                prf=1/mean(diff(t));
                
                [cs.CC,cc.CC,csx.CC,csy.CC] = CUSEShearWaveSpeed2DCCAHParFor(dataw,prf,factor.dx,factor.dx,direc,factor.stepSize,factor.windSizeCC,factor.tInterFac);
 
                
                saveName1=strrep(saveNameT1,'Ac',['L',num2str(pushN(uu)),'_']);

               save(fullfile(savePath,saveName1),'dataw','lat','axial','t',...
                                    'cs','csx','csy','cc','ccx','ccy','factor')     
            

            end
            toc
%             end        
        end
    end
end
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')