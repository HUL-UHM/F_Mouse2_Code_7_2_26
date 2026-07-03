clear all
close all
clc
addpath(genpath('D:\Code'))
% addpath(genpath('E:\SWEI\code\PaperFigCode'))
blatTemp=0.3*((1:128)-mean(1:128)); %[mm]
saveFig=1;

sTh='D:\Zac\2025.12.1\20kpa_5_2p61mmHg_P1\pSWEI\dispAll_F';
mPath=fullfile(sTh,'gV_RegCC_1');

saveImPath=fullfile(mPath,'results');

if ~exist(saveImPath,'dir')
    mkdir(saveImPath);
end

xlimVal=[-10 10];xTickVal=-10:5:10;
ylimVal=[10 30];axMed=[15 25];latMed=[-5 5];
yTickVal=[10:5:30];
% yTickName={'10';'20';'25';'30';'35'};

Sh1='';
posN='1';facN='2';acN='1';sSh=[Sh1,posN];
%% load group velocity
fN{1,1}=['35gV1p5LouL4_c800_e17_ang3_Pos',posN,'_L351_2_1_1_dirSongLR.mat'];
fN{2,1}=['35gV1p5LouL4_c800_e17_ang3_Pos',posN,'_L931_2_1_2_dirSongRL.mat'];
% fN{3,1}=['55gVLouIfL5F5p0_xx_c1400_ang3_Pos',posN,'_L834_dirSongRL.mat'];
% fN{4,1}=['55gVLouIfL5F5p0_xx_c1400_ang3_Pos',posN,'_L1044_dirSongRL.mat'];

rPV=[25 50];
kerNMm=[3.5 3.5];
threshMm=2.5;
swRL=0.1; swRH=7;
numCluster=2;downSamFac=1;flag.avg=1; % 1= mean; 2 =median;

fNum=size(fN,1);
load(fullfile(mPath,fN{1,1}),'cs','lat','axial')
if median(diff(lat))< 0.1
    lat=lat*1e3;
    axial=axial*1e3;
end

kerN=round(kerNMm./median(diff(lat)));
dx=median(diff(lat)); % mm
% tempLat=(xlimVal(1)-3.5):dx:(xlimVal(2)+3.5);
tempLat=blatTemp(1):dx:blatTemp(end);

axL=size(cs.CC,1);latL=numel(tempLat);
gVelAll_cc=nan([axL latL fNum]);

threshVal=round(threshMm/median(diff(lat)));
for ii=1:fNum
    load(fullfile(mPath,fN{ii,1}),'lat','axial','cs','csx','factor')
    if median(diff(lat))< 0.5
        lat=lat*1e3;
        axial= axial*1e3;
    end
    if contains(fN{ii,1},'RL')
            latInd=knnsearch(lat,min(lat)+4):knnsearch(lat,factor.pushLoc*1e3)-threshVal;
%         latInd=knnsearch(lat,xlimVal(1)-3.5):knnsearch(lat,factor.pushLoc*1e3)-threshVal;
    else
        latInd=knnsearch(lat,factor.pushLoc*1e3)+threshVal:knnsearch(lat,max(lat)-4); 
%         latInd=knnsearch(lat,factor.pushLoc*1e3)+threshVal:knnsearch(lat,xlimVal(2)+3.5); 
    end
    latSel=lat(latInd);
    latIndSel=knnsearch(tempLat',latSel);
    tempD=nanmedfilt2(cs.CC,kerN);
    gVelAll_cc(:,latIndSel,ii)=tempD(:,latInd);
end
gVelAll_cc(gVelAll_cc>swRH)=nan;
gVelAll_cc(gVelAll_cc<swRL)=nan;
gvAt=gVelAll_cc;

if numCluster>0%   
    [temgvgvP_cc,gVelAll_cc]=sumSWVimage(gVelAll_cc,numCluster,flag.avg);
else
    gVelAll_cc=squeeze(nanmean(gVelAll_cc,3));
end
%%
maskMed=zeros(size(gVelAll_cc(:,:,1)));
% axInd=knnsearch(axial,ylimVal(1)+1.5):knnsearch(axial,ylimVal(2)-1.5);
% latInd=knnsearch(tempLat',xlimVal(1)+1.5):knnsearch(tempLat',xlimVal(2)-1.5);
axInd=knnsearch(axial,axMed(1)):knnsearch(axial,axMed(2));
latInd=knnsearch(tempLat',latMed(1)):knnsearch(tempLat',latMed(2));
maskMed(axInd,latInd)=1;
%%
tempD=3*nanmedfilt2(gVelAll_cc.^2,kerN); 
h1=figure;
imagesc(tempLat,axial,tempD,rPV);
% hold all;contour(tempLat,axial,maskMed,[1 1],'k','lineWidth',3)
axis image
ylim(ylimVal)
xlim(xlimVal)
set(gca,'xTick',[xTickVal])
set(gca,'yTick',[yTickVal])
xlabel('Lat (mm)')
ylabel('Axial (mm)')
colormap jet
box on
setfigparms
set(gca,'fontname','arial') % used to be 16
setfigparms
if saveFig
    export_fig(fullfile(saveImPath,['gvP_',sSh]),'-png','-transparent','-r256','-painters','-rgb',h1)
end
medSWV=round(nanmedian(tempD(maskMed==1))*100)/100;
iqrSWV=round(iqr(tempD(maskMed==1))*100)/100;
meanSWV=round(nanmedian(tempD(maskMed==1))*100)/100;
stdSWV=round(nanstd(tempD(maskMed==1))*100)/100;
covSWV=stdSWV./meanSWV
save(fullfile(saveImPath,['gvP_ImQ',sSh]),'medSWV','iqrSWV','covSWV','stdSWV','meanSWV','gVelAll_cc');
disp('done')
