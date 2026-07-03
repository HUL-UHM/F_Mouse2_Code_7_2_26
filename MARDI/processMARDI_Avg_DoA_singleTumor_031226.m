clear all
close all
clc
addpath(genpath('D:\Code'))
fileNameSh='xCorL41p5rF0rF0Filt_25_';saveSh='Combine_xCorL41p5rF0rmF1.mat';

mainPath='F:\Mouse2\2026.03.09\37_102_Red_2p7mm_D26\MARDI\Fun\T_76_Fs_2';
figName='untitled_mardi.fig';
stepIntv=2.7;
rotNum=7;
acqNum=zeros(rotNum,1);
kerN=9;
axNum=110;
axCheck=[8.5 16]; % background range
axRange = [11 24]; % image range with tumor


cd(mainPath)
tempDir=dir;
tempDir(~[tempDir.isdir]) = [];  %remove non-directories
stepNum=size(tempDir,1)-2;
latTemp=0:stepIntv:(stepNum-1)*stepIntv;

fTot=10;
PeakCom=nan(axNum,rotNum,stepNum);
P2P_inteCom=nan(axNum,rotNum,fTot,stepNum);
PeakVCom=nan(axNum,rotNum,stepNum);
P2PCom=nan(axNum,rotNum,fTot,stepNum);

DoA.PeakCom=nan(axNum,stepNum);
DoA.P2P_inteCom=nan(axNum,stepNum,fTot);

DoA.PeakVCom=nan(axNum,stepNum);
DoA.P2PCom=nan(axNum,stepNum,fTot);

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
        Peak=max(dataF(axInd,:,1:26),[],3);
        PeakV=max(diff(dataF(axInd,:,1:26),1,3),[],3);

        PeakCom(1:size(Peak,1),jj,ii-2)=movmean(nanmedian(Peak,2),kerN);
        PeakVCom(1:size(PeakV,1),jj,ii-2)=movmean(nanmedian(PeakV,2),kerN);

        P2P_inteCom(1:size(Peak,1),jj,:,ii-2)=movmean(squeeze(nanmedian(avg_p2p_inte(axInd,:,:,:,:),2)),kerN);
        P2PCom(1:size(PeakV,1),jj,:,ii-2)=movmean(squeeze(nanmedian(avg_p2p(axInd,:,:,:,:),2)),kerN);
    end
    loopDim=3;
    for fI=1:2*fTot+2
        if fI < 11
            tempVal=outlierRemove(squeeze(P2P_inteCom(1:size(Peak,1),:,fI,ii-2)),loopDim,axIndCheck);
            P2P_inteCom(1:size(Peak,1),:,fI,ii-2)=tempVal;
        
        elseif fI > 11 && fI < 21
            tempVal=outlierRemove(squeeze(P2PCom(1:size(Peak,1),:,fI-10,ii-2)),loopDim,axIndCheck);
            P2PCom(1:size(Peak,1),:,fI-10,ii-2)=tempVal;
        elseif fI==21
            tempVal=outlierRemove(squeeze(PeakCom(:,:,ii-2)),loopDim,axIndCheck);
            PeakCom(1:size(Peak,1),:,ii-2)=tempVal;
        elseif fI==22
            tempVal=outlierRemove(squeeze( PeakVCom(:,:,ii-2)),loopDim,axIndCheck);
             PeakVCom(1:size(Peak,1),:,ii-2)=tempVal;
        end
        minVal=min(tempVal,[],2,"omitnan");
        maxVal=max(tempVal,[],2,"omitnan");
        if fI < 11
            DoA.P2P_inteCom(1:size(Peak,1),ii-2,fI)=movmean(maxVal./minVal,kerN);
        elseif fI > 11 && fI < 21
            DoA.P2PCom(1:size(Peak,1),ii-2,fI-10)=movmean(maxVal./minVal,kerN);
        elseif fI==21
            DoA.PeakCom(1:size(Peak,1),ii-2)=movmean(maxVal./minVal,kerN);
        elseif fI==22
            DoA.PeakVCom(1:size(Peak,1),ii-2)=movmean(maxVal./minVal,kerN);
        end
    end
end
 for fI=1:2*fTot+2
        if fI < 11
            tempData=DoA.P2P_inteCom(:,:,fI);        
        elseif fI > 11 && fI < 21
            tempData=DoA.P2PCom(:,:,fI-10); 
        elseif fI==21
            tempData=DoA.PeakCom; 
        elseif fI==22
            tempData=DoA.PeakVCom; 
        end
        tempData(tempData>5)=nan;
        prc10th=prctile(tempData',10)';
        prc90th=prctile(tempData',90)';
        tempData(tempData<prc10th)=nan;
        tempData(tempData>prc90th)=nan;
        if fI < 11
            DoA.P2P_inteCom(:,:,fI)=inpaint_nans(tempData,4);        
        elseif fI > 11 && fI < 21
            DoA.P2PCom(:,:,fI-10)=inpaint_nans(tempData,4);
        elseif fI==21
            DoA.PeakCom=inpaint_nans(tempData,4);
        elseif fI==22
           DoA.PeakVCom=inpaint_nans(tempData,4);
        end
 end

latDoA=min(latTemp):0.1:max(latTemp);
axDoA=min(axial):0.1:max(axial);
[xN,yN]=meshgrid(latDoA,axDoA);
[xO,yO]=meshgrid(latTemp,axial);
tempPI=zeros([size(xN) fTot]);
tempP=zeros([size(xN) fTot]);
for fI=1:2*fTot+2
        if fI < 11
            tempPI(:,:,fI)=nanmedfilt2(interp2(xO,yO,DoA.P2P_inteCom(:,:,fI),xN,yN,'cubic'),[11 11]);       
        elseif fI > 11 && fI < 21
            tempP(:,:,fI-10)=nanmedfilt2(interp2(xO,yO,DoA.P2PCom(:,:,fI-10),xN,yN,'cubic'),[11 11]);       
        elseif fI==21
            DoA.PeakCom=nanmedfilt2(interp2(xO,yO,DoA.PeakCom,xN,yN,'cubic'),[11 11]);
        elseif fI==22
            DoA.PeakVCom=nanmedfilt2(interp2(xO,yO,DoA.PeakVCom,xN,yN,'cubic'),[11 11]);
        end
end
DoA.P2P_inteCom=tempPI;
DoA.P2PCom=tempP;
%% show the image overlaid on Bmode
saveFig=1;
figPath=fileparts(fileparts(fileparts(mainPath)));
savePath=fullfile(figPath,'xCor_Avg');
if ~exist(savePath,'dir')
    mkdir(savePath)
end

uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
freqVal=100:100:1000;
for fI=1:2*fTot+2
        if fI < 11
            doaVal=DoA.P2P_inteCom(:,:,fI);tName=['DoA:P2PI-',num2str(freqVal(fI))];tSh=['_P2PI_',num2str(freqVal(fI))];
        elseif fI > 11 && fI < 21
            doaVal=DoA.P2PCom(:,:,fI-10);tName=['DoA:P2P-',num2str(freqVal(fI-10))];tSh=['_P2P_',num2str(freqVal(fI-10))];
        elseif fI==21
           doaVal=DoA.PeakCom;tName='DoA: Peak';tSh='_Peak_';
        elseif fI==22
            doaVal=DoA.PeakVCom;tName='DoA: PeakV';tSh='_PeakV_';
        end
tranParam.translateMm(2)=1; % 1 mm ish difference between our data and Verasonics Bmode
yR_2=[1 5];
h1=figure;
arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');hold all;axis image
xlabel('Lat (mm)');ylabel('Axial (mm)');title(tName)
box on;setfigparms; set(gca,'fontname','arial') % used to be 16
if saveFig
    export_fig(fullfile(savePath,['Bmode_DoA',tSh]),'-png','-transparent','-r128','-painters','-rgb',h1)       
end
end
save(fullfile(savePath,saveSh),'axial','latDoA','axDoA','P2PCom','PeakVCom','P2P_inteCom','PeakCom','DoA')

%% shift DoA image to match our ROI
close all
clc
addpath(genpath('D:\Code'))
saveFig=1;
shitMm=5.6;
mainPath='F:\tumor\tumor_rot_1p5_mm\ARFI\Fun\T_76_Fs_2\';
figPath='F:\tumor\tumor_rot_1p5_mm';
figName='untitled.fig';
uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
load('F:\tumor\tumor_rot_1p5_mm\ARFI\Fun\T_76_Fs_2\Combine_xCorL41p5rF0rmF1.mat')
latDoA=latDoA-shfitMm;
[xx,yy]= meshgrid(latDoA,axDoA);
load(fullfile(figPath,'BmodeMat_Fun','ROI_Final'))

latR=lat;
axR=axial;
[xxR,yyR]=meshgrid(latR,axR);
roiMed=interp2(xxR,yyR,double(roiMed),xx,yy,'nearest')>0;
roiLtCor=interp2(xxR,yyR,double(roiLtCor),xx,yy,'nearest')>0;
roiRtCor=interp2(xxR,yyR,double(roiRtCor),xx,yy,'nearest')>0;
roiUpCor=interp2(xxR,yyR,double(roiUpCor),xx,yy,'nearest')>0;
roiBtCor=interp2(xxR,yyR,double(roiBtCor),xx,yy,'nearest')>0;
roiWholeCor=interp2(xxR,yyR,double(roiWholeCor),xx,yy,'nearest')>0;
kidBouInterp=interp2(xxR,yyR,double(kidBouInterp),xx,yy,'nearest')>0;

yR_2=[0 4];
doaValAll=cell(6,1);
note ='6x6 : 6 rows = PD_1, PD_2, PD_avg, PV_1, PV_2, PV_avg,; 6 column = 6 rois; 1 = medulla, 2 = Left cortex, 3= right cortex, 4= upper cortex, 5 = bottom cortex, 6= whole cortex';
for fTot=1:6
    if fTot==1
        doaVal=DoA.PeakCom;tName='DoA: Peak_1';tSh='_Peak_1';
    elseif fTot==2
        doaVal=DoA.P2P_inteCom;tName='DoA: Peak_2';tSh='_Peak_2';
    elseif fTot==3
        doaVal=DoA.Peak_Avg;tName='DoA: Peak_Avg';tSh='_Peak_Avg';
    elseif fTot==4
        doaVal=DoA.PeakVCom;tName='DoA: PeakV_1';tSh='_PeakV_1';
    elseif fTot==5
        doaVal=DoA.P2PCom;tName='DoA: PeakV_2';tSh='_PeakV_2';
    elseif fTot==6
        doaVal=DoA.PeakV_Avg;tName='DoA: PeakV_Avg';tSh='_PeakV_Avg';
    end
    doaVal(doaVal>4)=nan;
    doaVal=inpaint_nans(doaVal,4);
    doaVal=nanmedfilt2(doaVal,[11 11]);
    doaValAll{fTot,1}=doaVal;
    h1=figure;
    arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');hold all;axis image
    hold all;axis image
    contour(latDoA,axDoA,roiMed,[1 1],'k--','LineWidth',3)
    xlabel('Lat (mm)');ylabel('Axial (mm)')
    box on;setfigparms; set(gca,'fontname','arial') % used to be 16 

    dataMed=doaVal;
    medP2P(fTot,1)=round(nanmedian(dataMed(roiMed==1))*100)/100;
    medP2P(fTot,2)=round(nanmedian(dataMed(roiLtCor==1))*100)/100;
    medP2P(fTot,3)=round(nanmedian(dataMed(roiRtCor==1))*100)/100;
    medP2P(fTot,4)=round(nanmedian(dataMed(roiUpCor==1))*100)/100;
    medP2P(fTot,5)=round(nanmedian(dataMed(roiBtCor==1))*100)/100;
    medP2P(fTot,6)=round(nanmedian(dataMed(roiWholeCor==1))*100)/100;

    meanP2P(fTot,1)=round(nanmean(dataMed(roiMed==1))*100)/100;
    meanP2P(fTot,2)=round(nanmean(dataMed(roiLtCor==1))*100)/100;
    meanP2P(fTot,3)=round(nanmean(dataMed(roiRtCor==1))*100)/100;
    meanP2P(fTot,4)=round(nanmean(dataMed(roiUpCor==1))*100)/100;
    meanP2P(fTot,5)=round(nanmean(dataMed(roiBtCor==1))*100)/100;
    meanP2P(fTot,6)=round(nanmean(dataMed(roiWholeCor==1))*100)/100;

    iqrP2P(fTot,1)=round(iqr(dataMed(roiMed==1))*100)/100;
    iqrP2P(fTot,2)=round(iqr(dataMed(roiLtCor==1))*100)/100;
    iqrP2P(fTot,3)=round(iqr(dataMed(roiRtCor==1))*100)/100;
    iqrP2P(fTot,4)=round(iqr(dataMed(roiUpCor==1))*100)/100;
    iqrP2P(fTot,5)=round(iqr(dataMed(roiBtCor==1))*100)/100;
    iqrP2P(fTot,6)=round(iqr(dataMed(roiWholeCor==1))*100)/100;

    stdP2P(fTot,1)=round(nanstd(dataMed(roiMed==1))*100)/100;
    stdP2P(fTot,2)=round(nanstd(dataMed(roiLtCor==1))*100)/100;
    stdP2P(fTot,3)=round(nanstd(dataMed(roiRtCor==1))*100)/100;
    stdP2P(fTot,4)=round(nanstd(dataMed(roiUpCor==1))*100)/100;
    stdP2P(fTot,5)=round(nanstd(dataMed(roiBtCor==1))*100)/100;
    meanP2P(fTot,6)=round(nanstd(dataMed(roiWholeCor==1))*100)/100;
  if saveFig
        export_fig(fullfile(mainPath,['ShiftBmode_DoA',tSh]),'-png','-transparent','-r128','-painters','-rgb',h1)  
  end
end
  save(fullfile(mainPath,'ImQ_DoA'),'medP2P','iqrP2P','stdP2P','meanP2P','note','latDoA','axDoA','doaValAll');
