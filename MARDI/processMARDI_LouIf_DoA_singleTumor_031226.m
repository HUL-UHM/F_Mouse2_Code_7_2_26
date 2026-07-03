clear all
close all
clc
addpath(genpath('D:\Code'))
drawMask=1;
saveFig=1;
yR_2=[1 3];
tranParam.translateMm(2)=.95; % 1 mm ish difference between our data and Verasonics Bmode
fileNameSh='LouIfL41p5rF0Lg1rF1Filt_25_';saveSh='Combine_LouIfL41p5rF0Lg1rF1Filt.mat';
yLV=[10 26]; xLV=[-4 12];
yTk=yLV(1):4:yLV(2);
xTk=xLV(1):4:xLV(2);

mainPath='F:\Mouse2\2026.02.23\38_62_Blue_2mm_D12\MARDI\Fun\T_76_Fs_2';
figName='untitled_mardi.fig';
stepIntv=2.0;
rotNum=7;
acqNum=zeros(rotNum,1);
kerN=9;
axNum=110;
axCheck=[15.5 21.5]; % background range
axRange = [10 24]; % image range with tumor


cd(mainPath)
tempDir=dir;
tempDir(~[tempDir.isdir]) = [];  %remove non-directories
stepNum=size(tempDir,1)-2;
latTemp=0:stepIntv:(stepNum-1)*stepIntv;

fTot=10;pushNum=2;
PeakCom=nan(axNum,pushNum,rotNum,stepNum);
PeakVCom=nan(axNum,pushNum,rotNum,stepNum);
DoA.PeakCom=nan(axNum,pushNum,stepNum);
DoA.PeakVCom=nan(axNum,pushNum,stepNum);


P2P_inteCom=nan(axNum,pushNum,rotNum,fTot,stepNum);
P2PCom=nan(axNum,pushNum,rotNum,fTot,stepNum);
DoA.P2P_inteCom=nan(axNum,pushNum,stepNum,fTot);
DoA.P2PCom=nan(axNum,pushNum,stepNum,fTot);

indFile=0;
for ii=3:size(tempDir,1)
    cd(fullfile(mainPath,tempDir(ii).name))
    tempDir1=dir;
    tempDir1(~[tempDir1.isdir]) = [];  %remove non-directories

    for jj=3:size(tempDir1,1)
        indT=strfind(tempDir1(jj).name,'_Ang');
        acqNum(jj-2,1)=str2double(tempDir1(jj).name(4:indT-1));
    end
    [acqNum,indS]=sort(acqNum);
    for jj=1:numel(acqNum)
        if acqNum(jj)<10
            fileName1=[fileNameSh,'0',num2str(acqNum(jj)),'.mat'];
            fileName2=['0',num2str(acqNum(jj))];

        else
            fileName1=[fileNameSh,num2str(acqNum(jj)),'.mat'];
            fileName2=num2str(acqNum(jj));
        end
        indT=strfind(tempDir1(indS(jj)+2).name,'_Ang');
        folName=['Acq',fileName2,tempDir1(indS(jj)+2).name(indT:end)];
        loadFile=fullfile(mainPath,tempDir(ii).name,folName,'Process_Fun','filtDisp2D',fileName1);
        load(loadFile,'axial','avg_p2p','avg_p2p_inte','dataF')
        axInd=knnsearch(axial,axRange(1)):knnsearch(axial,axRange(2));
        axial=axial(axInd);
        axIndCheck=knnsearch(axial,axCheck(1)):knnsearch(axial,axCheck(2));
        if contains(fileNameSh,'If')
            PeakV=max(dataF(axInd,:,1:26),[],3);
            Peak=max(cumsum(dataF(axInd,:,1:26),3),[],3);
        else
            Peak=max(dataF(axInd,:,1:26),[],3);
            PeakV=max(diff(dataF(axInd,:,1:26),1,3),[],3);
        end

        PeakCom(1:size(Peak,1),:,jj,ii-2)=movmean(Peak,kerN);
        PeakVCom(1:size(PeakV,1),:,jj,ii-2)=movmean(PeakV,kerN);

        P2P_inteCom(1:size(Peak,1),:,jj,:,ii-2)=movmean(squeeze(avg_p2p_inte(axInd,:,:,:,:)),kerN);
        P2PCom(1:size(PeakV,1),:,jj,:,ii-2)=movmean(squeeze(avg_p2p(axInd,:,:,:,:)),kerN);
    end
    loopDim=2;
    for fI=1:2*fTot+2
        if fI < 11
            tempVal=outlierRemove(squeeze(P2P_inteCom(:,:,:,fI,ii-2)),loopDim,axIndCheck);
            P2P_inteCom(:,:,:,fI,ii-2)=tempVal;
        elseif fI > 10 && fI < 21
            tempVal=outlierRemove(squeeze(P2PCom(:,:,:,fI-10,ii-2)),loopDim,axIndCheck);
            P2PCom(:,:,:,fI-10,ii-2)=tempVal;
        elseif fI==21
            tempVal=outlierRemove(squeeze(PeakCom(:,:,:,ii-2)),loopDim,axIndCheck);
            PeakCom(:,:,:,ii-2)=tempVal;
        elseif fI==22
            tempVal=outlierRemove(squeeze(PeakVCom(:,:,:,ii-2)),loopDim,axIndCheck);
            PeakVCom(:,:,:,ii-2)=tempVal;
        end
        minVal=squeeze(min(tempVal,[],3,"omitnan"));
        maxVal=squeeze(max(tempVal,[],3,"omitnan"));
        if fI < 11
            DoA.P2P_inteCom(1:size(maxVal,1),:,ii-2,fI)=movmean(maxVal./minVal,kerN);
        elseif fI > 10 && fI < 21
            DoA.P2PCom(1:size(maxVal,1),:,ii-2,fI-10)=movmean(maxVal./minVal,kerN);
        elseif fI==21
            DoA.PeakCom(1:size(maxVal,1),:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fI==22
            DoA.PeakVCom(1:size(maxVal,1),:,ii-2)=movmean(maxVal./minVal,kerN);
        end
    end
end
latDoA=min(latTemp):0.1:max(latTemp);
axDoA=min(axial):0.1:max(axial);
[xN,yN]=meshgrid(latDoA,axDoA);
[xO,yO]=meshgrid(latTemp,axial);
tempP2pI=zeros([size(xN) pushNum fTot]);
tempP2p=zeros([size(xN) pushNum fTot]);
tempP=zeros([size(xN) pushNum]);
tempPv=zeros([size(xN) pushNum]);
for mm=1:pushNum
     for fI=1:2*fTot+2
            if fI < 11
                tempData=squeeze(DoA.P2P_inteCom(:,mm,:,fI));        
            elseif fI > 10 && fI < 21
                tempData=squeeze(DoA.P2PCom(:,mm,:,fI-10)); 
            elseif fI==21
                tempData=squeeze(DoA.PeakCom(:,mm,:)); 
            elseif fI==22
                tempData=squeeze(DoA.PeakVCom(:,mm,:)); 
            end
            tempData(tempData>5)=nan;
            prc10th=prctile(tempData',10)';
            prc90th=prctile(tempData',90)';
            tempData(tempData<prc10th)=nan;
            tempData(tempData>prc90th)=nan;
            tempFil=inpaint_nans(tempData,4);
            if fI < 11 
                tempP2pI(:,:,mm,fI)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);      
            elseif fI > 10 && fI < 21
                tempP2p(:,:,mm,fI-10)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);      
            elseif fI==21
                tempP(:,:,mm)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);  
            elseif fI==22
               tempPv(:,:,mm)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);  
            end
     end
end
DoA.P2P_inteCom=tempP2pI;
DoA.P2PCom=tempP2p;
DoA.PeakCom=tempP;
DoA.PeakVCom=tempPv;
%% show the image overlaid on Bmode

figPath=fileparts(fileparts(fileparts(mainPath)));
savePath=fullfile(figPath,'LouIf_F');
if ~exist(savePath,'dir')
    mkdir(savePath)
end

uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
freqVal=100:100:1000;
for mm=1:pushNum+1
    for fI=1:2*fTot+2
        if mm<pushNum+1
            if fI < 11
                doaVal=DoA.P2P_inteCom(:,:,mm,fI);
                tName=[num2str(mm),'DoA:P2PI-',num2str(freqVal(fI))];
                tSh=[num2str(mm),'_P2PI_',num2str(freqVal(fI))];
            elseif fI > 10 && fI < 21
                doaVal=DoA.P2PCom(:,:,mm,fI-10);
                tName=[num2str(mm),'DoA:P2P-',num2str(freqVal(fI-10))];
                tSh=[num2str(mm),'_P2P_',num2str(freqVal(fI-10))];
            elseif fI==21
                doaVal=DoA.PeakCom(:,:,mm);
                tName=[num2str(mm),'DoA: Peak'];
                tSh=[num2str(mm),'_Peak_'];
            elseif fI==22
                doaVal=DoA.PeakVCom(:,:,mm);
                tName=[num2str(mm),'DoA: PeakV'];
                tSh=[num2str(mm),'_PeakV_'];
            end
        else
            if fI < 11
                doaVal=squeeze(mean(DoA.P2P_inteCom(:,:,:,fI),3));
                tName=['AvgDoA:P2PI-',num2str(freqVal(fI))];
                tSh=['Avg_P2PI_',num2str(freqVal(fI))];
            elseif fI > 10 && fI < 21
                doaVal=squeeze(mean(DoA.P2PCom(:,:,:,fI-10),3));
                tName=['AvgDoA:P2P-',num2str(freqVal(fI-10))];
                tSh=['Avg_P2P_',num2str(freqVal(fI-10))];
            elseif fI==21
                doaVal=squeeze(mean(DoA.PeakCom,3));
                tName='AvgDoA: Peak';
                tSh='Avg_Peak_';
            elseif fI==22
                doaVal=squeeze(mean(DoA.PeakVCom,3));
                tName='AvgDoA: PeakV';
                tSh='Avg_PeakV_';
            end
        end

    if fI==1 && mm==1 && drawMask==1
        figure;
        arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');hold all;axis image
        xlabel('Lat (mm)');ylabel('Axial (mm)');title(tName)
        box on;setfigparms; set(gca,'fontname','arial') % used to be 16
        title('Draw Mask')
        roiMask=roipoly;
        title('Draw Tumor Boundary')
        roiHMI =roipoly;
        close all
    end
 
    if drawMask
        h1=figure;
        arfiBmodeOverlay(latDoA,axDoA,doaVal.*double(roiMask),yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');
        hold all;axis image
        contour(latDoA,axDoA,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
    else
        h1=figure;
        arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');
    end

    xlabel('Lat (mm)');ylabel('Axial (mm)');title(tName)
    ylim(yLV);set(gca,'yTick',yTk)
    xlim(xLV);set(gca,'xTick',xTk)
    box on;setfigparms; set(gca,'fontname','arial') % used to be 16
    if saveFig
        export_fig(fullfile(savePath,[tSh,'Bmode_DoA']),'-png','-transparent','-r612','-painters','-rgb',h1)       
    end
    end

    if mm<pushNum+1
        close all
    end
end
save(fullfile(savePath,saveSh),'roiHMI','roiMask','axial','latDoA','axDoA','P2PCom','PeakVCom','P2P_inteCom','PeakCom','DoA')

%% shift DoA image to match our ROI
close all
clc
addpath(genpath('D:\Code'))
saveFig=1;
shitMm=5.6;
mainPath='D:\Mouse\2025.10.20_DoA\17_D10_Rot\ARFI\Fun\T_76_Fs_2';
figPath='D:\Mouse\2025.10.20\17_D10\MARDI\Fun\T_76_Fs_2\f100_n10_nC4_PRF10_PI\Fac_2_HP_4_6_5_nC_4_20-October-2025_20-17-43';
figName='35_HARF_x2f1Fac_2_HP_4_6_5_nC_4_Pos2_Ac1.fig';
uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));

load('D:\Mouse\2025.10.20_DoA\17_D10_Rot\ARFI\Fun\T_76_Fs_2\Combine_LouIfL41p5rF0Lg2rF1.mat')
latDoA=latDoA-shitMm;
[xx,yy]= meshgrid(latDoA,axDoA);
load(fullfile(figPath,'BmodeMat_Fun','ROI_Final'))
bkdROINum=3;
latR=lat;
axR=axial;
[xxR,yyR]=meshgrid(latR,axR);
roiMask=interp2(xxR,yyR,double(roiMask),xx,yy,'nearest')>0;
roiHMI=interp2(xxR,yyR,double(roiHMI),xx,yy,'nearest')>0;
roiBKD_rect=interp2(xxR,yyR,double(roiBKD_rect),xx,yy,'nearest')>0;
roiINC_rect=interp2(xxR,yyR,double(roiINC_rect),xx,yy,'nearest')>0;
yR_2=[0 5];
doaValAll=cell(6,1);
for fTot=1:6
    if fTot==1
        doaVal=DoA.Peak_1;tName='DoA: Peak_1';tSh='_Peak_1';
    elseif fTot==2
        doaVal=DoA.Peak_2;tName='DoA: Peak_2';tSh='_Peak_2';
    elseif fTot==3
        doaVal=DoA.Peak_Avg;tName='DoA: Peak_Avg';tSh='_Peak_Avg';
    elseif fTot==4
        doaVal=DoA.PeakV_1;tName='DoA: PeakV_1';tSh='_PeakV_1';
    elseif fTot==5
        doaVal=DoA.PeakV_2;tName='DoA: PeakV_2';tSh='_PeakV_2';
    elseif fTot==6
        doaVal=DoA.PeakV_Avg;tName='DoA: PeakV_Avg';tSh='_PeakV_Avg';
    end
    doaVal(doaVal>5)=nan;
    doaVal=inpaint_nans(doaVal,4);
    doaVal=nanmedfilt2(doaVal,[11 11]);
    doaValAll{fTot,1}=doaVal;
    h1=figure;
    arfiBmodeOverlay(latDoA,axDoA,doaVal.*double(roiMask),yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');hold all;axis image
    hold all;axis image
    contour(latDoA,axDoA,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
    xlabel('Lat (mm)');ylabel('Axial (mm)')
    box on;setfigparms; set(gca,'fontname','arial') % used to be 16 

    data_med=doaVal;
    data_med(~roiMask)=nan;
    medDoA.Abs(fTot,1)=round(nanmedian(data_med(roiHMI==1))*100)/100;
    meanDoA.Abs(fTot,1)=round(nanmean(data_med(roiHMI==1))*100)/100;
    iqrDoA.Abs(fTot,1)=round(iqr(data_med(roiHMI==1))*100)/100;
    stdDoA.Abs(fTot,1)=round(nanstd(data_med(roiHMI==1))*100)/100;
    
    if bkdROINum<3       
        bMed.rectAbs(fTot,1)=nanmedian(data_med(roiBKD_rect==bkdROINum));
        bVar.rectAbs(fTot,1)=nanvar(data_med(roiBKD_rect==bkdROINum));  
    else
       
        bMed.rectAbs(fTot,1)=nanmedian(data_med(roiBKD_rect==1))+nanmedian(data_med(roiBKD_rect==2));
        bVar.rectAbs(fTot,1)=nanvar(data_med(roiBKD_rect==1))+nanvar(data_med(roiBKD_rect==2));  
    end

    iMed.rectAbs(fTot,1)=nanmedian(data_med(roiINC_rect==1));
    iVar.rectAbs(fTot,1)=nanvar(data_med(roiINC_rect==1));

    if saveFig
        export_fig(fullfile(mainPath,['ShiftBmode_DoA',tSh]),'-png','-transparent','-r128','-painters','-rgb',h1)  
    end
end
        save(fullfile(mainPath,'ImQ_DoA'),'medDoA','iqrDoA','stdDoA','meanDoA','bMed','iMed','bVar','iVar','latDoA','axDoA','doaValAll','roiHMI');
