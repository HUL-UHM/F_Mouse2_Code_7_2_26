clear all
close all
clc
addpath(genpath('D:\Code'))

mainPath='D:\Mouse\2025.10.20_DoA\17_D10_Rot\ARFI\Fun\T_76_Fs_2';
fileNameSh='LouIfL41p5rF0Lg2rF1_35_';saveSh='Combine_LouIfL41p5rF0Lg2rF1.mat';
stepIntv=0.75;
rotNum=8;
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
             p2p_norm=PeakCom_1(:,:,ii-2);
        elseif fTot==2
             p2p_norm=PeakCom_2(:,:,ii-2);
        elseif fTot==3
             p2p_norm=0.5*(PeakCom_1(:,:,ii-2)+PeakCom_1(:,:,ii-2));
        elseif fTot==4
             p2p_norm=PeakVCom_1(:,:,ii-2);
        elseif fTot==5
             p2p_norm=PeakVCom_2(:,:,ii-2);
        elseif fTot==6
            p2p_norm=0.5*(PeakVCom_1(:,:,ii-2)+PeakVCom_1(:,:,ii-2));
        end

        prc10th=prctile(p2p_norm(:),10);
        prc90th=prctile(p2p_norm(:),90);
        p2p_norm(p2p_norm<prc10th)=nan;
        p2p_norm(p2p_norm>prc90th)=nan;
        minVal=min(p2p_norm,[],2,"omitnan");
        maxVal=max(p2p_norm,[],2,"omitnan");
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

save(fullfile(mainPath,saveSh),'latDoA','axDoA','PeakVCom_2','PeakVCom_1','PeakCom_2','PeakCom_1','DoA')
%% show the image overlaid on Bmode
saveFig=0;
figPath='D:\Mouse\2025.10.20_DoA\17_D10_Rot\ARFI\Fun\T_76_Fs_2';
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
yR_2=[0 6];
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
