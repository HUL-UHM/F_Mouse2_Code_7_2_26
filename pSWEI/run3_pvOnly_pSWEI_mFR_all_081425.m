clear all
% close all
clc
addpath(genpath('D:\Lab_Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='D:\Zac\2025.11.20';
blat= 0.2985*((1:128)-mean(1:128));
factor.core=10; % number of core for parallel processing
%%
factor.axL=[5 30];
factor.interpolate=0;
factor.dx=0.1e-3;factor.spaDs=1;
factor.axAvg=1e-3; %  [m] axial average ROI in 
factor.tAvg=0.2e-3; %[ms] axial average ROI in t
factor.tWind1=0.3;% tuckey window value
factor.tWind2=0.3;% tuckey window value

factor.windSize=4.0e-3; % window size in mm;
factor.sTsize=4.0e-3;
factor.downSamplePhaseCalc=5; %[sample] downsample factor for speeding up phase velocity calc                             calcculatin
% if Flag.dirFilt
%     factor.dirAngles=0; % [rad] angle in direct fitlering  
% end
Flag.dirFilt_selFreq=1; % flag for direction filter
if Flag.dirFilt_selFreq
    factor.Dir_order=2; % bandpass order
    factor.Dir_power=2; % bandpass power
    factor.Dir_angles=0; % directional angle
    Flag.bpFilter=0; % flag for bandpass filter
    factor.Dir_cutoffs=[0.2 8]; % in m/s
    factor.sTsize=4.0e-3;
    factor.searchFreqRg=[-40 40]; % range for frequency band
end

Flag.dirFilt=0; % flag for direction filter
Flag.linearReg=0;
Flag.LinSong=0;Flag.LinSelF=~Flag.LinSong;
Flag.dirFiltSong=0;
Flag.dirFiltAdaptive=0; % flag for adaptive filter with FFT
factor.interpolate=0;
% Flag.bpFilter=0; % flag for bandpass filter

%%
factor.tInterFac=2;
factor.NT=round(3072*factor.tInterFac);
flagFilt=1;

nameSh='_F1p5LouL4';
rSh=['pV',nameSh(3:end)];
foldSh='pSWEI';
if contains(nameSh,'If')
    factor.calDiff=0; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.rmReverbFram=0;% number of frame to remove after push
else
    factor.calDiff=1; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.samDiff=4; % difference sample
  factor.rmReverbFram=1;% number of frame to remove after push
end

factor.findFreq=1;

for mm=1:size(main_path_all,1)
    main_path=main_path_all{mm,1}
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;
    file_name_all=[];
    for kk=dummy_var:length(temp)

        file_name=subdir(fullfile(main_path,temp(kk).name,['*',nameSh,'*_c1*ang*Pos*.mat']));
%           file_name=subdir(fullfile(main_path,temp(kk).name,['*',nameSh,'*F4*.mat']));
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
        if contains(filePath,'f150_10') || contains(filePath,'f150_n10')
            factor.harmFreq=150:150:1500;% initial oscillation frequency
        elseif contains(filePath,'f75_14') || contains(filePath,'f75_n14')
            factor.harmFreq=75:75:1200;% initial oscillation frequency
        elseif contains(filePath,'f100_n10') || contains(filePath,'f100_10')
            factor.harmFreq=100:100:1000;% initial oscillation frequency
        else
            factor.harmFreq=100:100:1000;% initial oscillation frequency
        end
        if contains(fileName,'_nAng_1')
            factor.rmReverbFram=3;% number of frame to remove after push
        end
        savePath=fullfile(filePath,['pV_FFTselFrM_',num2str(factor.rmReverbFram),'_2']);
        if ~exist(savePath,'dir')
            mkdir(savePath)
        end
        saveName=strrep(fileName,nameSh,rSh);
%         saveName=saveName(1:end-2);
        load(file_name_all{1,bb},'aimgIni');
        
            
        if Flag.dirFilt_selFreq
            saveNameT=[saveName,'_SelF'];
        elseif Flag.linearReg==1
            if Flag.LinSong
                saveNameT=[saveName,'_RegSong'];
            else
                saveNameT=[saveName,'_RegSelF'];
            end
        elseif Flag.dirFiltAdaptive
            saveNameT=[saveName,'_Adap'];
        elseif Flag.dirFiltSong
             saveNameT=[saveName,'_dirSong'];
        else
            saveNameT=saveName;
            flagFilt=0;
        end
        
        if numel(aimgIni.swPushLoc)>1 && flagFilt 
            pushL=aimgIni.swPushLoc(end);
            if pushL>64
                saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushL)]),'Pi.mat'];
            else
                saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushL)]),'0.mat'];
            end
        else
            pushL=aimgIni.swPushLoc(1);
            if pushL>64
                saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushL)]),'Pi.mat'];
            else
                saveNameCh=[strrep(saveNameT,'Ac',['L',num2str(pushL)]),'0.mat'];
            end
        end
%         if ~exist(fullfile(savePath,[saveNameT1,'Pi.mat']),'file') && contains(filePath,foldSh) && contains(fileName,'nAng_3')
        if ~exist(fullfile(savePath,saveNameCh),'file')  %&& contains(fileName,'Pos1')%&& contains(filePath,'f100_n10_14') && contains(fileName,'Pos2')
            %if contains(filePath,'pwHSWEI')  
                tic
                load(file_name_all{1,bb},'dispData','lat','axial','t','aimgIni');
                
                for uu=1:size(dispData,1)
                    if uu==1
                        axialOld=axial;latOld=lat;tOld=t;
                    else
                        axial=axialOld;lat=latOld;t=tOld;
                    end
                    
                    pushL=aimgIni.swPushLoc(uu);                   
                    arfidata=dispData{uu,1};
                    factor.pushLoc=blat(pushL)*1e-3;
                    factor.latCheck=round(pushL*numel(lat)/128);
                    aimg=aimgIni;           
                    factor.SubPreFrTime=0;
                    aI=knnsearch(axial,factor.axL(1)):knnsearch(axial,factor.axL(2));
                    if uu==1
                        aimgIni.PRFKHz=aimgIni.PRFKHz*factor.tInterFac;
                    end
                    if factor.calDiff
                        [data,axial,lat,t]=rmRevIsoPixel(arfidata(aI,:,:),axial(aI),lat,t,factor,aimgIni);
                    else
                        [data,axial,lat,t]=rmRevIsoPixelIF(arfidata(aI,:,:),axial(aI),lat,t,factor,aimgIni);
                    end  
                    [b,a]=butter(2,2*[50 2500]*median(diff(t)),'bandpass'); 
                    data=permute(filtfilt(b,a,permute(data,[3 1 2])),[2 3 1]);
                    
                    FD=aimgIni.focalDepthmm*1e-3;
                    
                    if Flag.dirFilt~=0 ||  Flag.dirFilt_selFreq==1
                        if pushL<64
                             factor.Dir_angles=0; %directional angle
                             saveNameT1=[saveNameT,'0'];
                        else
                             factor.Dir_angles=pi; %directional angle
                             saveNameT1=[saveNameT,'Pi'];
                        end
                    else

                         saveNameT1=saveNameT;

                    end 
                    if factor.tAvg<1
                        factor.tAvg =round(factor.tAvg/median(diff(t)));
                    else 
                        factor.tAvg =round(0.3e-3/median(diff(t)));
                    end
                    saveName1=strrep(saveNameT1,'Ac',['L',num2str(pushL)]);
                    tic

%                     [phaseVel,selFreq,fftAmp,mu,eta,dim,factor]=calc_pV_selF_DirBP_2DFFTwAvg(data,lat,axial,t,FD,factor);
%                         else   
                            [phaseVel,fftAmp,mu,eta,dim,factor]=calc_phaseV_DirBP_2DFFTwAvg(data,lat,axial,t,FD,factor);
%                         end
                    toc
                    axial=dim.axial*1e3;
                    lat=dim.lat*1e3;
                    t=dim.t*1e3;
                    dS=factor.spaDs;aI=1:dS:numel(axial);lI=1:dS:numel(lat);
                    axial=axial(aI);lat=lat(lI);dim.axial= axial*1e-3;dim.lat=lat*1e-3;
                    mu.Fmin=mu.Fmin(aI,lI);mu.CF=mu.CF(aI,lI);
                    eta.Fmin=eta.Fmin(aI,lI);eta.CF=eta.CF(aI,lI);
                    fftAmp=fftAmp(aI,lI,:,:);phaseVel=phaseVel(aI,lI,:,:);
                    savefast(fullfile(savePath,saveName1),'phaseVel','fftAmp','mu','eta','dim','factor','Flag','lat','axial')
                    
                end
                toc
%             end        
        end
    end
end
%%
close all
YRp=[0.5 4];
xlimVal=[-10 10];
ylimVal=[20 40];
yTickVal=ylimVal(1):5:ylimVal(2);
xTickVal=xlimVal(1):5:xlimVal(2);
kerN=round([2.1e-3 2.1e-3]./factor.dx);
if mod(kerN,2)==0
    kerN=kerN+1;
end
YRyoung=[5 45];
YRv=[0 5];

for ll=1:size(phaseVel,3)
    figure('position',[40 335 560 420]);
    imagesc(lat,axial,nanmedfilt2(squeeze(phaseVel(:,:,ll,1)),kerN),YRp);
    axis image
    ylim(ylimVal)
    xlim(xlimVal)
    set(gca,'xTick',[xTickVal])
    set(gca,'yTick',[yTickVal])
    colormap jet
    title(['LR: ', num2str(factor.harmFreq(ll))])
    box on
    setfigparms
end

disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')