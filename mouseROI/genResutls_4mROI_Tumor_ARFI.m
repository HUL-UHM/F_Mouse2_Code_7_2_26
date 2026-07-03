clear all
close all
clc
addpath(genpath('D:\Code\functions'))
saveFig=1;
%% 
mpTR='F:\Mouse2\2026.02.16\38_61_Red_D4\MARDI\Fun\T_76_Fs_2\f100_n10_nC4_PRF10_PI\';
main_pathR=fullfile(mpTR,'Fac_2_HP_4_6_5_nC_4_16-February-2026_18-50-04');

mpT='F:\Mouse2\2026.02.16\38_61_Red_D4\ARFI\Fun\T_76_Fs_2\';
main_path=fullfile(mpT,'400_16-February-2026_18-47-38');

figName='25_ARF_x2_cycle_400_Pos1_Ac1.fig';
yLV=[5 25]; xLV=[-10 10];

bkdROINum=1; % 1=left only, 2= right only, 3=both;

yTk=yLV(1):5:yLV(2);
xTk=xLV(1):5:xLV(2);

axL=[16 20.5]; % change
latRg=[-10.5 -8.3]; % change

% hF1='LouIfL41p5rF0Lg2rF1_35_01.mat';
% hF2='LouIfL41p5rF0Lg2rF1_35_02.mat';
hF1='xCorL41p5rF0rmF1_25_01.mat';
hF2='xCorL41p5rF0rmF1_25_01.mat';

roiName='ROI_Final';
LogThreshold=-60;
addLogThresh=10;
saveImPath=fullfile(main_path,'Process_Fun','Results_xCor_R');
if ~exist(saveImPath,'dir')
    mkdir(saveImPath)
end

load(fullfile(main_pathR,'BmodeMat_Fun',roiName))
latR=lat;
axR=axial;
[xxR,yyR]=meshgrid(latR,axR);
tumorSizeMm2=sum(roiHMI(:)).*0.1*0.1;
%% load Verasonics Bmode with ROI
uiopen(fullfile(main_path,figName),1);
set(gcf,'position',[1   236   472   730])
[x,y,bmodeVer] = getimage(gcf);
baxialVer=linspace(y(1),y(2),size(bmodeVer,1));
blatVer=linspace(x(1),x(2),size(bmodeVer,2));
h1=figure;
imagesc(blatVer,baxialVer+tranParam.translateMm(2),bmodeVer);colormap gray;
hold all;axis image
contour(latR,axR,roiMask,[1 1],'m','LineWidth',3)
contour(latR,axR,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
contour(latR,axR,roiINC_rect.*double(roiMask),[1 1],'r','LineWidth',3)
if bkdROINum<3
    contour(latR,axR,double(roiBKD_rect==bkdROINum).*double(roiMask),[1 1],'b','LineWidth',3) 
else
    contour(latR,axR,double(roiBKD_rect).*double(roiMask),[1 1],'b','LineWidth',3)
end
xlabel('Lat (mm)');ylabel('Axial (mm)')
ylim(yLV);set(gca,'yTick',yTk)
xlim(xLV);set(gca,'xTick',xTk)
box on;setfigparms; set(gca,'fontname','arial') % used to be 16
if saveFig
    export_fig(fullfile(saveImPath,'Bmode_Ver'),'-png','-transparent','-r128','-painters','-rgb',h1)
end
%% Das bmode
bName=strrep(figName,'ARF','BmodeF1p5');
bName=strrep(bName,'.fig','.mat');
bmodePath=fullfile(main_path,'BmodeMat_Fun');
load(fullfile(bmodePath,bName))
bmodeLog= 20*log10(bmode/max(bmode(:)));
bmodeLog(bmodeLog<LogThreshold)=LogThreshold;
h1=figure;imagesc(blat,baxial,bmodeLog,[LogThreshold+addLogThresh 0]);colormap gray
hold all;axis image
contour(latR,axR,roiMask,[1 1],'m','LineWidth',3)
contour(latR,axR,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
contour(latR,axR,roiINC_rect.*double(roiMask),[1 1],'r','LineWidth',3)
if bkdROINum<3
    contour(latR,axR,double(roiBKD_rect==bkdROINum).*double(roiMask),[1 1],'b','LineWidth',3) 
else
    contour(latR,axR,double(roiBKD_rect).*double(roiMask),[1 1],'b','LineWidth',3)
end
xlabel('Lat (mm)');ylabel('Axial (mm)')
ylim(yLV);set(gca,'yTick',yTk)
xlim(xLV);set(gca,'xTick',xTk)
box on;setfigparms; set(gca,'fontname','arial') % used to be 16
if saveFig
    export_fig(fullfile(saveImPath,'Bmode_Das'),'-png','-transparent','-r128','-painters','-rgb',h1)
end
%% load HMI data
load(fullfile(main_path,'Process_Fun',hF1),'Peak','PeakV','lat','axial')
[xx,yy]= meshgrid(lat,axial);
PeakI1=interp2(xx,yy,Peak,xxR,yyR,'spline');
PeakVI1=interp2(xx,yy,PeakV,xxR,yyR,'spline');
load(fullfile(main_path,'Process_Fun',hF2),'Peak','PeakV','lat','axial')
[xx,yy]= meshgrid(lat,axial);
PeakI2=interp2(xx,yy,Peak,xxR,yyR,'spline');
PeakVI2=interp2(xx,yy,PeakV,xxR,yyR,'spline');

axInd=knnsearch(axR,axL(1)):knnsearch(axR,axL(2));
kerN=[round(1.1/median(diff(axR))),round(1.1./median(diff(latR)))];
normLat=knnsearch(latR,latRg(1)):knnsearch(latR,latRg(2));
yR_2=[0.0 4]; 
%%
note ='2x1 : 2 rows = PV, PD';

medP2P.Norm=zeros(2,1); medP2P.Abs=zeros(2,1);
meanP2P.Norm=zeros(2,1);meanP2P.Abs=zeros(2,1);
iqrP2P.Norm=zeros(2,1);iqrP2P.Abs=zeros(2,1);
stdP2P.Norm=zeros(2,1);stdP2P.Abs=zeros(2,1);

for fTot=1:2
    if fTot==1
         data_med1=nanmedfilt2(PeakVI1,kerN);
         data_med2=nanmedfilt2(PeakVI2,kerN);
         tName='PV   ';
    elseif fTot==2
         data_med1=nanmedfilt2(PeakI1,kerN);
         data_med2=nanmedfilt2(PeakI2,kerN);
         tName='PD   ';
    end   
    normAL=movmean(median(data_med1(axInd,normLat),2),kerN(1));
    [fitresult] = disp_ax_gauss2_fit(axR(axInd),normAL./max(normAL),startPoint,1);
    normAL1=feval(fitresult,axR);
    normMat=repmat(normAL1,1,size(data_med1,2));
    p2p_norm1=nanmedfilt2(data_med1./normMat,kerN);  
    normAL=movmean(median(data_med2(axInd,normLat),2),kerN(1));
    [fitresult] = disp_ax_gauss2_fit(axR(axInd),normAL./max(normAL),startPoint,1);
    normAL1=feval(fitresult,axR);
    normMat=repmat(normAL1,1,size(data_med1,2));
    p2p_norm2=nanmedfilt2(data_med2./normMat,kerN);  

    p2p_norm=mean(cat(3,p2p_norm2,p2p_norm1),3);
    data_med=mean(cat(3,data_med2,data_med1),3);        

    p2p_norm(~roiMask)=nan;
    data_med(~roiMask)=nan;
    
    medP2P.Norm(fTot,1)=round(nanmedian(p2p_norm(roiHMI==1))*100)/100;
    meanP2P.Norm(fTot,1)=round(nanmean(p2p_norm(roiHMI==1))*100)/100;
    iqrP2P.Norm(fTot,1)=round(iqr(p2p_norm(roiHMI==1))*100)/100;
    stdP2P.Norm(fTot,1)=round(nanstd(p2p_norm(roiHMI==1))*100)/100;
    
    medP2P.Abs(fTot,1)=round(nanmedian(data_med(roiHMI==1))*100)/100;
    meanP2P.Abs(fTot,1)=round(nanmean(data_med(roiHMI==1))*100)/100;
    iqrP2P.Abs(fTot,1)=round(iqr(data_med(roiHMI==1))*100)/100;
    stdP2P.Abs(fTot,1)=round(nanstd(data_med(roiHMI==1))*100)/100;
    
    if bkdROINum<3
        bMed.rect(fTot,1)=nanmedian(p2p_norm(roiBKD_rect==bkdROINum));
        bVar.rect(fTot,1)=nanvar(p2p_norm(roiBKD_rect==bkdROINum));  
        
        bMed.rectAbs(fTot,1)=nanmedian(data_med(roiBKD_rect==bkdROINum));
        bVar.rectAbs(fTot,1)=nanvar(data_med(roiBKD_rect==bkdROINum));  
    else
        bMed.rect(fTot,1)=nanmedian(p2p_norm(roiBKD_rect==1))+nanmedian(p2p_norm(roiBKD_rect==2));
        bVar.rect(fTot,1)=nanvar(p2p_norm(roiBKD_rect==1))+nanvar(p2p_norm(roiBKD_rect==2)); 
        
        bMed.rectAbs(fTot,1)=nanmedian(data_med(roiBKD_rect==1))+nanmedian(data_med(roiBKD_rect==2));
        bVar.rectAbs(fTot,1)=nanvar(data_med(roiBKD_rect==1))+nanvar(data_med(roiBKD_rect==2));  
    end

    iMed.rect(fTot,1)=nanmedian(p2p_norm(roiINC_rect==1));
    iVar.rect(fTot,1)=nanvar(p2p_norm(roiINC_rect==1));
    
    iMed.rectAbs(fTot,1)=nanmedian(data_med(roiINC_rect==1));
    iVar.rectAbs(fTot,1)=nanvar(data_med(roiINC_rect==1));

    imQ.Con_rect(fTot,1)=abs(squeeze(iMed.rect(fTot,1))-squeeze(bMed.rect(fTot,1)))...
                            ./squeeze(bMed.rect(fTot,1));
    imQ.CNR_rect(fTot,1)=abs(squeeze(iMed.rect(fTot,1))-squeeze(bMed.rect(fTot,1)))...
                            ./sqrt(squeeze(bVar.rect(fTot,1))+squeeze(iVar.rect(fTot,1)));

    medVal=nanmedian(p2p_norm(:));
    madVal=mad(p2p_norm(:));
    lowRan=max(0,medVal-2.5*madVal);
    hiRan=medVal+2.5*madVal;
    yR_2=[lowRan,hiRan];
  
    h1=figure;imagesc(lat,axial,(p2p_norm).*double(roiMask),yR_2);
    hold all;axis image;colorbar
    contour(latR,axR,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
    xlabel('Lat (mm)');ylabel('Axial (mm)')
    ylim(yLV);set(gca,'yTick',yTk)
    xlim([min(lat) max(lat)])
    colormap jet;title(tName)
    box on;setfigparms
    % colorbar('off')
    set(gca,'fontname','arial') % used to be 16
    if saveFig
        export_fig(fullfile(saveImPath,['P2PD_',tName(1:end-3)]),'-png','-transparent','-r128','-painters','-rgb',h1)       
    end
    h1=figure;
    arfiBmodeOverlay(lat,axial,(p2p_norm).*double(roiMask),yR_2,blat,baxial,bmodeLog,0.3,' ');
    hold all;axis image
    contour(latR,axR,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
    xlabel('Lat (mm)');ylabel('Axial (mm)')
    ylim(yLV);set(gca,'yTick',yTk)
    xlim(xLV);set(gca,'xTick',xTk)
    box on;setfigparms; set(gca,'fontname','arial') % used to be 16
    colormap jet;title(tName)
    % colorbar('off')
    if saveFig
        export_fig(fullfile(saveImPath,['Bmode_P2PD_',tName(1:end-3)]),'-png','-transparent','-r128','-painters','-rgb',h1)       
    end


    medVal=nanmedian(data_med(:));
    madVal=mad(data_med(:));
    lowRan=max(0,medVal-2.5*madVal);
    hiRan=medVal+2.5*madVal;
    yR_3=[lowRan,hiRan];

   
    h11=figure;
    arfiBmodeOverlay(lat,axial,data_med.*double(roiMask),yR_3,blat,baxial,bmodeLog,0.4,' ');
    hold all;axis image
    contour(latR,axR,roiHMI.*double(roiMask),[1 1],'k--','LineWidth',3)
    xlabel('Lat (mm)');ylabel('Axial (mm)')
    ylim(yLV);set(gca,'yTick',yTk)
    xlim(xLV);set(gca,'xTick',xTk)
    box on;setfigparms; set(gca,'fontname','arial') % used to be 16
    colormap jet;title(tName)
    if saveFig
        export_fig(fullfile(saveImPath,['Bmode_AbsP2PD_',tName(1:end-3)]),'-png','-transparent','-r128','-painters','-rgb',h11)       
    end

end
% figure(h1)
% colorbar 
% export_fig(fullfile(saveImPath,['MF_colorbar']),'-png','-transparent','-r256','-painters','-rgb',h1)

vStat=regionprops(roiHMI,'EquivDiameter','Area');
diaVal=vStat.EquivDiameter.*mean(diff(lat));
areaVal=vStat.Area.*mean(diff(lat).^2);
save(fullfile(saveImPath,'ImQ_P2PD'),'medP2P','iqrP2P','stdP2P','meanP2P','bMed','iMed','bVar','iVar','imQ','vStat','diaVal','areaVal','latR','axR','roiHMI');