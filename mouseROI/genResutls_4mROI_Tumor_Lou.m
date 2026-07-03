clear all
close all
clc
saveFig=1;
%% 
mpT='F:\Mouse2\2026.03.22\37_120_Violet_D39\MARDI\Fun\T_76_Fs_2\f100_n10_nC4_PRF10_PI';
main_path=fullfile(mpT,'Fac_2_HP_4_6_5_nC_4_22-March-2026_10-03-27');
figName='25_HARF_x2f1Fac_2_HP_4_6_5_nC_4_Pos1_Ac1.fig';
yLV=[5 25]; xLV=[-10 10];

bkdROINum=1; % 1=left only, 2= right only, 3=both;

yTk=yLV(1):5:yLV(2);
xTk=xLV(1):5:xLV(2);

axL=[14 18]; % change
latRg=[-6 -4]; % change

hF1='LouIfL41p5rF0Lg1rF0Filt_35_01.mat';
hF2='LouIfL41p5rF0Lg1rF0Filt_35_02.mat';
hFSlope1='LouIfL41p5rF0Lg1rF0FiltLinFit_35_01.mat';
hFSlope2='LouIfL41p5rF0Lg1rF0FiltLinFit_35_02.mat';
roiName='ROI_Final1';
LogThreshold=-60;
addLogThresh=10;
saveImPath=fullfile(main_path,'Process_Fun','filtDisp2D','Results_Lou');
if ~exist(saveImPath,'dir')
    mkdir(saveImPath)
end

load(fullfile(main_path,'BmodeMat_Fun',roiName))
latR=lat;
axR=axial;
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
bName=strrep(figName,'HARF','BmodeF1p5');
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
load(fullfile(main_path,'Process_Fun','filtDisp2D',hF1),'avg_p2p_inte','dataF')
avgData1=squeeze(avg_p2p_inte(:,:,1,:));dataF1=dataF;
load(fullfile(main_path,'Process_Fun','filtDisp2D',hF2),'avg_p2p_inte','lat','axial','factor','dataF')
avgData2=squeeze(avg_p2p_inte(:,:,1,:));
dataF2=dataF; clear dataF
axInd=knnsearch(axial,axL(1)):knnsearch(axial,axL(2));
kerN=[round(1.1/median(diff(axial))),round(1.1./median(diff(lat)))];
normLat=knnsearch(lat,latRg(1)):knnsearch(lat,latRg(2));
yR_2=[0.0 2]; 
%%
note ='13x1 : 13 rows = 10 frequency, PV, PD, and inverse slope';

medP2P.Norm=zeros(size(avgData1,3)+3,1); medP2P.Abs=zeros(size(avgData1,3)+3,1); 
meanP2P.Norm=zeros(size(avgData1,3)+3,1); meanP2P.Abs=zeros(size(avgData1,3)+3,1); 
iqrP2P.Norm=zeros(size(avgData1,3)+3,1);iqrP2P.Abs=zeros(size(avgData1,3)+3,1);
stdP2P.Norm=zeros(size(avgData1,3)+3,1);stdP2P.Abs=zeros(size(avgData1,3),1);
harmFreq=100:100:1e3;
for fTot=1:size(avgData1,3)+3
    if fTot<11
        data_med1=nanmedfilt2(avgData1(:,:,fTot),kerN);    
        data_med2=nanmedfilt2(avgData2(:,:,fTot),kerN);  
        tName=[num2str(factor.harmFreq(fTot)),' Hz'];
    elseif fTot==11
         data_med1=nanmedfilt2(max(dataF1(:,:,1:26),[],3),kerN);
         data_med2=nanmedfilt2(max(dataF2(:,:,1:26),[],3),kerN);
         tName='PV   ';
    elseif fTot==12
         data_med1=nanmedfilt2(max(cumsum(dataF1(:,:,1:26),3),[],3),kerN);
         data_med2=nanmedfilt2(max(cumsum(dataF2(:,:,1:26),3),[],3),kerN);
         tName='PD   ';
    elseif fTot==13
        load(fullfile(main_path,'Process_Fun','filtDisp2D',hFSlope1),'Param','GOF')
        invSlope=squeeze(1./Param.Lin1norm(:,:,1));
        R2=GOF.Lin1R2norm;
        invSlope(R2<0.7)=nan;
        invSlope=inpaint_nans(invSlope,4);
        data_med1=nanmedfilt2(invSlope,kerN);

        load(fullfile(main_path,'Process_Fun','filtDisp2D',hFSlope2),'Param','GOF')
        invSlope=squeeze(1./Param.Lin1norm(:,:,1));
        R2=GOF.Lin1R2norm;
        invSlope(R2<0.7)=nan;
        invSlope=inpaint_nans(invSlope,4);
        data_med2=nanmedfilt2(invSlope,kerN);
        tName='invSlope   ';
    end   
    normAL=movmean(median(data_med1(axInd,normLat),2),kerN(1));
    [fitresult] = disp_ax_gauss2_fit(axial(axInd),normAL./max(normAL),startPoint,1);
    normAL1=feval(fitresult,axial);
    normMat=repmat(normAL1,1,size(data_med1,2));
    p2p_norm1=nanmedfilt2(data_med1./normMat,kerN);  
    normAL=movmean(median(data_med2(axInd,normLat),2),kerN(1));
    [fitresult] = disp_ax_gauss2_fit(axial(axInd),normAL./max(normAL),startPoint,1);
    normAL1=feval(fitresult,axial);
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

    % data_med=p2p_norm.*double(roiMask);
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
