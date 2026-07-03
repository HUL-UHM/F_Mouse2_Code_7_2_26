clear all
% close all
clc
addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\Mouse2';

axL=[8 26];

% common parameter
factor.NT=4096;
factor.tAvg=3; %[sample] axial average ROI in t
factor.axAvg=1; %[mm] axial average ROI in t
factor.order=8;% filter order
factor.BW=1.0;% passband factor for filtering
factor.core=2; % number of core for parallel processing
factor.dx=0.1;
% for 1D code
factor.axInd=[15 20 25 32 35]; % get a sense of pass band
factor.filterName=2;% 1= butter, 2= IIR; 3 = FIR
factor.commonPB=0; % 1 = PB will generate from the latINd and axInd
factor.SubPreFrTime=0;

fSh='LouIfL41p5rF0Lg1rF1';
rSh=[fSh(1:end),'Filt'];
factor.flagInterp2D=0; % don't 2-D interpolate if lat dim is sparsely sampled
factor.flagInterp1D=~factor.flagInterp2D; % interpolate in axial dimension
if factor.flagInterp1D
   % for 2D code
    factor.windSize=[0.5 0]; % window size in mm [axial lat];
    factor.downSamplePhaseCalc=[5 1]; %[sample] downsample factor for speeding up displacement calc                             calcculatin
else
       % for 2D code
    factor.windSize=[0.7 0.7]; % window size in mm [axial lat];
    factor.downSamplePhaseCalc=[5 5]; %[sample] downsample factor for speeding up displacement calc                             calcculatin

end

if contains(fSh,'If')
    factor.calDiff=0; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.rmReverbFram=0;% number of frame to remove after push
else
    factor.calDiff=1; % 1 = xCorDisp, 0 =LoupassDisplacement
    factor.rmReverbFram=1;% number of frame to remove after push
    factor.diffSam=2; % lag for calculating the difference
    factor.CCthresh=0.95;
end

for mm=1:size(main_path_all,1)
    main_path=main_path_all{mm,1};
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;

    file_name_all=[];

    for kk=dummy_var:length(temp)

        file_name=subdir(fullfile(main_path,temp(kk).name,[fSh,'_','*.mat']));
        if ~isempty(file_name)

            file_name_cell=struct2cell(file_name);
            if kk==dummy_var && done_already==0;
                 done_already=1;
                file_name_all=file_name_cell(1,1:size(file_name_cell,2));
            else
                file_name_all=[file_name_all,file_name_cell(1,1:size(file_name_cell,2))];
            end
        else
           dummy_var=kk+1;
        end

    end
    %%
% parpool(factor.core)
    for bb=1:size(file_name_all,2)
         disp(['!!!! start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
        [filePath,fileName]=fileparts(file_name_all{1,bb});     
        [~,folName]=fileparts(fileparts(filePath));  
        saveName=strrep(fileName,fSh,rSh);
        savePath=fullfile(filePath,'filtDisp2D');
        if ~exist(fullfile(savePath,[saveName,'.mat']),'file') && contains(filePath,'MARDI')  && contains(filePath,'mm')
            if ~exist(savePath,'dir')
                mkdir(savePath)
            end
            if contains(fSh,'xC')
               load(file_name_all{1,bb},'arfidata','lat','axial','t','aimgIni','CCcoeff');
               factor.CC=CCcoeff;
            else
               load(file_name_all{1,bb},'arfidata','lat','axial','t','aimgIni');
            end
           tLen=1:numel(t);
            if median(diff(lat))<0
                lat=lat(end:-1:1,1);
                arfidata=arfidata(:,end:-1:1,:);
            end
             if numel(lat)>4  
                factor.latInd=lat(1:2:round((lat(end-2)-lat(1))/median(diff(lat)))); % get a sense of pass band
             else
                factor.latInd=1:numel(lat);
             end% get a sense of pass band
                indHARF=strfind(folName,'_nC');
            factor.NCycle=str2double(folName(indHARF+4));%: oscillation cycle number
            if isnan(factor.NCycle) || factor.NCycle==1
                factor.NCycle=str2double(folName(indHARF+4:indHARF+5));%: oscillation cycle number
            end
            if isnan(factor.NCycle) 
                factor.NCycle=4;%: oscillation cycle number
            end
            if contains(filePath,'f100_10') || contains(filePath,'f100_n10') 
                factor.harmFreq=100:100:1e3;% initial oscillation frequency
                factor.fundFreq=100;% fundamental oscillation frequency
            elseif contains(filePath,'f100_15') || contains(filePath,'f100_n15') 
                factor.harmFreq=100:100:1.5e3;% initial oscillation frequency
                factor.fundFreq=100;% fundamental oscillation frequency
            elseif contains(filePath,'f150_10') || contains(filePath,'f150_n10') 
                factor.harmFreq=150:150:1.5e3;% initial oscillation frequency
                factor.fundFreq=150;% fundamental oscillation frequency
            elseif contains(filePath,'f150_7') || contains(filePath,'f150_n7') 
                factor.harmFreq=150:150:1.05e3;% initial oscillation frequency
                factor.fundFreq=150;% fundamental oscillation frequency
            elseif contains(filePath,'f100_1') || contains(filePath,'f100_n1') 
                factor.harmFreq=100;% initial oscillation frequency
                factor.fundFreq=100;% fundamental oscillation frequency
            elseif contains(filePath,'f200_5') || contains(filePath,'f200_n5') 
                factor.harmFreq=200:200:1000;% initial oscillation frequency
                factor.fundFreq=200;% fundamental oscillation frequency
            elseif contains(filePath,'f300_5') || contains(filePath,'f300_n5') 
                factor.harmFreq=300:300:1500;% initial oscillation frequency
                factor.fundFreq=300;% fundamental oscillation frequency
            else
                factor.harmFreq=100:100:1000;% initial oscillation frequency
                factor.fundFreq=100;% fundamental oscillation frequency
            end

            if median(diff(lat))<0
                lat=lat(end:-1:1);
                arfidata=arfidata(:,end:-1:1,:);
            end
            ind=knnsearch(axial,axL(1)):knnsearch(axial,axL(2));
            axial=axial(ind);
            arfidata=arfidata(ind,:,:);
            if isfield(factor,'CC')
                factor.CC=factor.CC(ind,:,:);
            end
            if ~isnan(factor.NCycle) && ~isnan(factor.fundFreq)
                tic
                [avg_p2p,avg_p2p_inte,dataF,dim,passParams,dataFilt]=calc_P2P_disp_2D_3I(arfidata(:,:,tLen),axial,lat,t(tLen),aimgIni,factor);
                toc
                axial=dim.axial*1e3;
                if iscolumn(dim.lat)                    
                    lat=dim.lat*1e3;
                else
                    lat=dim.lat'*1e3;
                end
                t=dim.t*1e3;
                if contains(fSh,'xC')
                factor=rmfield(factor,'CC');
                end
                 Peak=max(dataF(:,:,1:25),[],3);
                savefast(fullfile(savePath,saveName),'Peak','dataFilt','avg_p2p','avg_p2p_inte','passParams','dataF','factor','dim','axial','lat','t')
            end
        end
       
    end
    
%  delete(gcp('nocreate'))
        disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])

end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
run4M_slope_MARDI_mFR_all_081425
