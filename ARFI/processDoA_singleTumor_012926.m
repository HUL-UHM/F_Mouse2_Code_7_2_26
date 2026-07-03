clear all
close all
clc
addpath(genpath('D:\Code'))

mainPath='F:\tumor\tumor_rot_1p5_mm\ARFI\Fun\T_76_Fs_2';
fileNameSh='xCorL41p5rF0rmF1_35_';saveSh='Combine_xCorL41p5rF0rmF1.mat';
% fileNameSh='xCorL41p5rF0rmF1_35_';saveSh='Combine_xCorL41p5rF0rmF1.mat';
stepIntv=1.5;
rotNum=7;
acqNum=zeros(rotNum,1);
kerN=43;
axNum=988;


cd(mainPath)
tempDir=dir;
tempDir(~[tempDir.isdir]) = [];  %remove non-directories
stepNum=size(tempDir,1)-2;
latTemp=0:stepIntv:(stepNum-1)*stepIntv;

PeakCom_1=nan(axNum,rotNum,stepNum);
PeakCom_2=nan(axNum,rotNum,stepNum);
PeakVCom_1=nan(axNum,rotNum,stepNum);
PeakVCom_2=nan(axNum,rotNum,stepNum);
DoA.Peak_1=nan(axNum,stepNum);
DoA.Peak_2=nan(axNum,stepNum);
DoA.Peak_Avg=nan(axNum,stepNum);
DoA.PeakV_1=nan(axNum,stepNum);
DoA.PeakV_2=nan(axNum,stepNum);
DoA.PeakV_Avg=nan(axNum,stepNum);
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
        PeakCom_1(1:size(Peak,1),jj,ii-2)=movmean(Peak(:,1),kerN);
        PeakCom_2(1:size(Peak,1),jj,ii-2)=movmean(Peak(:,2),kerN);
        PeakVCom_1(1:size(PeakV,1),jj,ii-2)=movmean(PeakV(:,1),kerN);
        PeakVCom_2(1:size(PeakV,1),jj,ii-2)=movmean(PeakV(:,2),kerN);
    end
    for fTot = 1:6
        if fTot==1
             dataMed=PeakCom_1(:,:,ii-2);
        elseif fTot==2
             dataMed=PeakCom_2(:,:,ii-2);
        elseif fTot==4
             dataMed=PeakVCom_1(:,:,ii-2);
        elseif fTot==5
             dataMed=PeakVCom_2(:,:,ii-2);
        end
        if fTot~=3 || fTot~=6
            prc10th=prctile(dataMed',15)';
            prc90th=prctile(dataMed',85)';
            dataMed(dataMed<prc10th)=nan;
            dataMed(dataMed>prc90th)=nan;
        elseif fTot==3
            dataMed1=PeakCom_1(:,:,ii-2);
            prc10th=prctile(dataMed1',15)';
            prc90th=prctile(dataMed1',85)';
            dataMed1(dataMed1<prc10th)=nan;
            dataMed1(dataMed1>prc90th)=nan;
            
            dataMed2=PeakCom_2(:,:,ii-2);
            prc10th=prctile(dataMed2',15)';
            prc90th=prctile(dataMed2',85)';
            dataMed2(dataMed2<prc10th)=nan;
            dataMed2(dataMed2>prc90th)=nan;

            dataMed=0.5*(dataMed1+dataMed2);
        elseif fTot==6
            dataMed1=PeakVCom_1(:,:,ii-2);
            prc10th=prctile(dataMed1',15)';
            prc90th=prctile(dataMed1',85)';
            dataMed1(dataMed1<prc10th)=nan;
            dataMed1(dataMed1>prc90th)=nan;
            
            dataMed2=PeakVCom_2(:,:,ii-2);
            prc10th=prctile(dataMed2',15)';
            prc90th=prctile(dataMed2',85)';
            dataMed2(dataMed2<prc10th)=nan;
            dataMed2(dataMed2>prc90th)=nan;

            dataMed=0.5*(dataMed1+dataMed2);
        end
        minVal=min(dataMed,[],2,"omitnan");
        maxVal=max(dataMed,[],2,"omitnan");
        
        if fTot==1
            DoA.Peak_1(:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==2
             DoA.Peak_2(:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==3
             DoA.Peak_Avg(:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==4
             DoA.PeakV_1(:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==5
            DoA.PeakV_2(:,ii-2)=movmean(maxVal./minVal,kerN);
        elseif fTot==6
           DoA.PeakV_Avg(:,ii-2)=movmean(maxVal./minVal,kerN);
        end
    end

end
for fTot=1:6
    if fTot==1
       tempData=DoA.Peak_1;
    elseif fTot==2
         tempData=DoA.Peak_1;
    elseif fTot==3
        tempData=DoA.Peak_Avg;
    elseif fTot==4
         tempData=DoA.PeakV_1;
    elseif fTot==5
        tempData=DoA.PeakV_2;
    elseif fTot==6
       tempData=DoA.PeakV_Avg;
    end
    prc10th=prctile(tempData',15)';
    prc90th=prctile(tempData',85)';
    tempData(tempData<prc10th)=nan;
    tempData(tempData>prc90th)=nan;
    if fTot==1
        DoA.Peak_1=inpaint_nans(tempData,4);
    elseif fTot==2
        DoA.Peak_2=inpaint_nans(tempData,4);
    elseif fTot==3
        DoA.Peak_Avg=inpaint_nans(tempData,4);
    elseif fTot==4
        DoA.PeakV_1=inpaint_nans(tempData,4);
    elseif fTot==5
        DoA.PeakV_2=inpaint_nans(tempData,4);
    elseif fTot==6
        DoA.PeakV_Avg=inpaint_nans(tempData,4);
    end
end

latDoA=min(latTemp):0.1:max(latTemp);
axDoA=min(axial):0.1:max(axial);
[xN,yN]=meshgrid(latDoA,axDoA);
[xO,yO]=meshgrid(latTemp,axial);


DoA.Peak_1=nanmedfilt2(interp2(xO,yO,DoA.Peak_1,xN,yN,'cubic'),[11 11]);
DoA.Peak_2=nanmedfilt2(interp2(xO,yO,DoA.Peak_2,xN,yN,'cubic'),[11 11]);
DoA.Peak_Avg=nanmedfilt2(interp2(xO,yO,DoA.Peak_Avg,xN,yN,'cubic'),[11 11]);
DoA.PeakV_1=nanmedfilt2(interp2(xO,yO,DoA.PeakV_1,xN,yN,'cubic'),[11 11]);
DoA.PeakV_2=nanmedfilt2(interp2(xO,yO,DoA.PeakV_2,xN,yN,'cubic'),[11 11]);
DoA.PeakV_Avg=nanmedfilt2(interp2(xO,yO,DoA.PeakV_Avg,xN,yN,'cubic'),[11 11]);

save(fullfile(mainPath,saveSh),'axial','latDoA','axDoA','PeakVCom_2','PeakVCom_1','PeakCom_2','PeakCom_1','DoA')
%% show the image overlaid on Bmode
saveFig=1;
figPath='F:\tumor\tumor_rot_1p5_mm';
figName='untitled.fig';
uiopen(fullfile(figPath,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
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
tranParam.translateMm(2)=1; % 1 mm ish difference between our data and Verasonics Bmode
yR_2=[0.8 1.5];
h1=figure;
arfiBmodeOverlay(latDoA,axDoA,doaVal,yR_2,blatVer,baxialVer+tranParam.translateMm(2),bmodeVer,0.3,' ');hold all;axis image
xlabel('Lat (mm)');ylabel('Axial (mm)')
box on;setfigparms; set(gca,'fontname','arial') % used to be 16
if saveFig
    export_fig(fullfile(mainPath,['Bmode_DoA',tSh]),'-png','-transparent','-r128','-painters','-rgb',h1)       
end
end
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
