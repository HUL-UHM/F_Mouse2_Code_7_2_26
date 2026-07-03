clear all
close all
clc
addpath(genpath('D:\Code'))
saveFig=1;
tranParam.translateMm(2)=1.0; % 1 mm ish difference between our data and Verasonics Bmode
yR_2=[1 3];
mainPath='F:\Mouse2\2026.02.17\38_125_Violet_1mm_D5\ARFI\Fun\T_76_Fs_2';
% fileNameSh='LouIfL41p5rF0Lg2rF1_30_';saveSh='Combine_LouIfL41p5rF0Lg2rF1.mat';folN='LouIf';
fileNameSh='xCorL41p5rF0rmF1_25_';saveSh='Combine_xCorL41p5rF0rmF1.mat';folN='xCor';
figName='untitled.fig';
yLV=[5 30]; xLV=[-7 20];
yTk=yLV(1):4:yLV(2);
xTk=xLV(1):4:xLV(2);

stepIntv=3;%
rotNum=8;%
pushNum=4;
acqNum=zeros(rotNum,1);
kerN=43;
axNum=100;
axCheck=[15.5 25]; % background range
axRange = [13 25]; % image range with tumor


cd(mainPath)
tempDir=dir;
tempDir(~[tempDir.isdir]) = [];  %remove non-directories
% tempDir(end)=[]; % remove the last translation as not all angle was collected
stepNum=size(tempDir,1)-2;
latTemp=0:stepIntv:(stepNum-1)*stepIntv;


PeakCom=nan(axNum,pushNum,rotNum,stepNum);
PeakVCom=nan(axNum,pushNum,rotNum,stepNum);
DoA.PeakCom=nan(axNum,pushNum,stepNum);
DoA.PeakVCom=nan(axNum,pushNum,stepNum);

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
        loadFile=fullfile(mainPath,tempDir(ii).name,folName,'Process_Fun',fileName1);
        load(loadFile,'axial','Peak','PeakV')
        latTot=1:size(Peak,2);
        axInd=knnsearch(axial,axRange(1)):knnsearch(axial,axRange(2));
        axial=axial(axInd);
        axIndCheck=knnsearch(axial,axCheck(1)):knnsearch(axial,axCheck(2));

        PeakCom(1:numel(axInd),latTot,jj,ii-2)=movmean(Peak(axInd,:),kerN);
        PeakVCom(1:numel(axInd),latTot,jj,ii-2)=movmean(PeakV(axInd,:),kerN);
    end
    loopDim=2;
    for fTot = 1:2
        if fTot==1
            tempVal=outlierRemove(squeeze(PeakCom(:,:,:,ii-2)),loopDim,axIndCheck);
            PeakCom(:,:,:,ii-2)=tempVal;
        elseif fTot==2
            tempVal=outlierRemove(squeeze(PeakVCom(:,:,:,ii-2)),loopDim,axIndCheck);
            PeakVCom(:,:,:,ii-2)=tempVal;
        end
        minVal=squeeze(min(tempVal,[],3,"omitnan"));
        maxVal=squeeze(max(tempVal,[],3,"omitnan"));
        if fTot==1
            DoA.PeakCom(1:size(maxVal,1),latTot,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==2
            DoA.PeakVCom(1:size(maxVal,1),latTot,ii-2)=movmean(maxVal./minVal,kerN);
        end
       
    end

end
pushNum=max(latTot);
latDoA=min(latTemp):0.1:max(latTemp);
axDoA=min(axial):0.1:max(axial);
axIndDoACheck=knnsearch(axDoA',axCheck(1)):knnsearch(axDoA',axCheck(2));
[xN,yN]=meshgrid(latDoA,axDoA);
[xO,yO]=meshgrid(latTemp,axial);
tempP2pI=zeros([size(xN) pushNum fTot]);
tempP2p=zeros([size(xN) pushNum fTot]);
tempP=zeros([size(xN) pushNum]);
tempPv=zeros([size(xN) pushNum]);
for mm=1:pushNum
     for fI=1:2
            if fI==1 
                tempData=squeeze(DoA.PeakCom(:,mm,:)); 
            elseif fI==2
                tempData=squeeze(DoA.PeakVCom(:,mm,:)); 
            end
            tempData(tempData>5)=nan;
            prc10th=prctile(tempData',10)';
            prc90th=prctile(tempData',90)';
            tempData(tempData<prc10th)=nan;
            tempData(tempData>prc90th)=nan;
            tempFil=inpaint_nans(tempData,4);
            if fI==1
                tempP(:,:,mm)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);  
            elseif fI==2
               tempPv(:,:,mm)=nanmedfilt2(interp2(xO,yO,tempFil,xN,yN,'cubic'),[11 11]);  
            end
     end
end
DoA.PeakCom=tempP;
DoA.PeakVCom=tempPv;
%% show the image overlaid on Bmode

figPath=fileparts(fileparts(fileparts(mainPath)));
savePath=fullfile(figPath,folN);
if ~exist(savePath,'dir')
    mkdir(savePath)
end
uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
for mm=1:pushNum+1
    for fI=1:2
        if mm<pushNum+1
            if fI==1
                doaVal=DoA.PeakCom(:,:,mm);
                tName=[num2str(mm),'DoA: Peak'];
                tSh=[num2str(mm),'_Peak_'];
            elseif fI==2
                doaVal=DoA.PeakVCom(:,:,mm);
                tName=[num2str(mm),'DoA: PeakV'];
                tSh=[num2str(mm),'_PeakV_'];
            end
        else
            if fI==1
         
                doaVal=squeeze(mean(DoA.PeakCom,3));
                 % doaVal=nanmedfilt2(squeeze(nanmean(outlierRemove(DoA.PeakCom,loopDim,axIndDoACheck),3)),[11 11]);
                tName='AvgDoA: Peak';
                tSh='Avg_Peak_';
            elseif fI==2
                doaVal=squeeze(mean(DoA.PeakVCom,3));
                % doaVal=nanmedfilt2(squeeze(nanmean(outlierRemove(DoA.PeakVCom,loopDim,axIndDoACheck),3)),[11 11]);
                tName='AvgDoA: PeakV';
                tSh='Avg_PeakV_';
            end
        end
        h1=figure;
        arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');
        xlabel('Lat (mm)');ylabel('Axial (mm)');title(tName)
        ylim(yLV);set(gca,'yTick',yTk)
        xlim(xLV);set(gca,'xTick',xTk)
        box on;setfigparms; set(gca,'fontname','arial') % used to be 16
        if saveFig
            export_fig(fullfile(savePath,[tSh,'Bmode_DoA']),'-png','-transparent','-r612','-painters','-rgb',h1)       
        end
    end
end
save(fullfile(savePath,saveSh),'axial','latDoA','axDoA','PeakVCom','PeakCom','DoA')

%% shift DoA image to match our ROI
close all
clc
roiPath='F:\Mouse2\2026.02.17\38_125_Violet_D5\MARDI\Fun\T_76_Fs_2\f100_n10_nC4_PRF10_PI\Fac_2_HP_4_6_5_nC_4_22-March-2026_11-10-04';
figName='25_HARF_x2f1Fac_2_HP_4_6_5_nC_4_Pos1_Ac1.fig';
shitMm=6.6;%
xLV=[-10 12];
xTk=xLV(1):4:xLV(2);

%latDoA=-latDoA+shitMm; %if untitled.fig is not available
latDoA=latDoA-shitMm;
[xx,yy]= meshgrid(latDoA,axDoA);
uiopen(fullfile(roiPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));

load(fullfile(roiPath,'BmodeMat_Fun','ROI_Final'))
latR=lat;
axR=axial;
[xxR,yyR]=meshgrid(latR,axR);
roiHMI=interp2(xxR,yyR,double(roiHMI),xx,yy,'nearest')>0;
roiBKD_rect=interp2(xxR,yyR,double(roiBKD_rect),xx,yy,'nearest')>0;
roiINC_rect=interp2(xxR,yyR,double(roiINC_rect),xx,yy,'nearest')>0;
roiMask=interp2(xxR,yyR,double(roiMask),xx,yy,'nearest')>0;
doaValAll=cell(6,1);
note ='pushNum x Feature x RoI';
pushNum=size(DoA.PeakCom,3);
for mm=1:pushNum+1
    for fI=1:2
        if mm<pushNum+1
            if fI==1
                doaVal=DoA.PeakCom(:,:,mm);
                tName=[num2str(mm),'DoA: Peak'];
                tSh=[num2str(mm),'_Peak_'];
            elseif fI==2
                doaVal=DoA.PeakVCom(:,:,mm);
                tName=[num2str(mm),'DoA: PeakV'];
                tSh=[num2str(mm),'_PeakV_'];
            end
        else
            if fI==1
                doaVal=squeeze(mean(DoA.PeakCom,3));
                tName='AvgDoA: Peak';
                tSh='Avg_Peak_';
            elseif fI==2
                doaVal=squeeze(mean(DoA.PeakVCom,3));
                tName='AvgDoA: PeakV';
                tSh='Avg_PeakV_';
            end
        end

        h1=figure;
        arfiBmodeOverlay(latDoA,axDoA,doaVal.*double(roiMask),yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');
        hold all;
        contour(latDoA,axDoA,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
        xlabel('Lat (mm)');ylabel('Axial (mm)');title(tName)
        ylim(yLV);set(gca,'yTick',yTk)
        xlim(xLV);set(gca,'xTick',xTk)
        box on;setfigparms; set(gca,'fontname','arial') % used to be 16
        if saveFig
            export_fig(fullfile(savePath,['Shift_',tSh,'Bmode_DoA']),'-png','-transparent','-r612','-painters','-rgb',h1)       
        end
    
        dataMed=doaVal;
        medP2P(mm,fI,1)=round(nanmedian(dataMed(roiHMI==1))*100)/100;
        medP2P(mm,fI,2)=round(nanmedian(dataMed(roiBKD_rect==1))*100)/100;
        medP2P(mm,fI,3)=round(nanmedian(dataMed(roiINC_rect==1))*100)/100;
    
        meanP2P(mm,fI,1)=round(nanmean(dataMed(roiHMI==1))*100)/100;
        meanP2P(mm,fI,2)=round(nanmean(dataMed(roiBKD_rect==1))*100)/100;
        meanP2P(mm,fI,3)=round(nanmean(dataMed(roiINC_rect==1))*100)/100;
    
        iqrP2P(mm,fI,1)=round(iqr(dataMed(roiHMI==1))*100)/100;
        iqrP2P(mm,fI,2)=round(iqr(dataMed(roiBKD_rect==1))*100)/100;
        iqrP2P(mm,fI,3)=round(iqr(dataMed(roiINC_rect==1))*100)/100;
    
        stdP2P(mm,fI,1)=round(nanstd(dataMed(roiHMI==1))*100)/100;
        stdP2P(mm,fI,2)=round(nanstd(dataMed(roiBKD_rect==1))*100)/100;
        stdP2P(mm,fI,3)=round(nanstd(dataMed(roiINC_rect==1))*100)/100;
    end

end
save(fullfile(savePath,'DoA_ImQ'),'medP2P','iqrP2P','stdP2P','meanP2P','note');
