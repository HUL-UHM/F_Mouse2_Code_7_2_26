clear all
% close all
clc
addpath(genpath('D:\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='D:\Zac\2025.12.1';
axL=[5 30];
 
flagDisp=1; % 1 = Loup, 2 = Loup inter Frame, 3= Xcorr
fSh='dasF1p5';
if flagDisp ==1
    rSh=[fSh(4:7),'LouL'];lag=1;dwS=1;
elseif flagDisp ==2
    rSh=[fSh(1:2),'LouIf'];
    lag=2;rmFr=3;
elseif flagDisp ==3
    rSh=[fSh,'CorL'];dwS=1;
end
    
knLenAll = [4];
interpFac=4;
coreNumDisp=10;

for oo=1:numel(knLenAll)
    knLen=knLenAll(oo);
for mm=1:size(main_path_all,1)

main_path=main_path_all{mm,1}
cd(main_path)
temp=dir;
temp(~[temp.isdir]) = [];  %remove non-directories
dummy_var=3;
done_already=0;

file_name_all=[];

for kk=dummy_var:length(temp)

    file_name=subdir(fullfile(main_path,temp(kk).name,['*_',fSh,'_*Ac*.mat']));
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
% parpool(coreNumDisp)
% parfor bb=1:size(file_name_all,2)
for bb=1:size(file_name_all,2)
    disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
    [filePath,fileName]=fileparts(file_name_all{1,bb});
    savePath=fullfile(filePath,'dispAll_F');
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    if contains(rSh,'If')
        saveName=strrep(fileName,fSh,[rSh,num2str(lag),'Rm',num2str(rmFr),'L',num2str(knLen)]);
    else
        saveName=strrep(fileName,fSh,[rSh,num2str(knLen)]);
    end
        % &&
        %(contains(filePath,'phan3_inc1_Dia6')... % && contains(fileName,'Pos2') 
            %
        tempData=load(file_name_all{1,bb});
        scatSum_all=tempData.scatSum;
        aimg=tempData.aimgIni;
        aimgIni=tempData.aimgIni;      
        aimg.kernelLengthWavs = knLen;       
        aimg.maxDisplacement =300e-6; 
%         if exist('axL','var')
        aimg.axL=axL;
        aimgIni.axL=axL;
        
%         dispData=cell(aimg.SWIFrames,1);
        dispData=cell(1,1);
        RF=cell(1,1);
        CC=cell(1,1);
        pushLoc=aimgIni.swPushLoc;
        aimgIni.sweiFrames=1;
        for ii=1:aimg.SWIFrames
            saveName1=[saveName,'_',num2str(ii)];
            if ~exist(fullfile(savePath,[saveName1,'.mat']),'file')   %&&  contains(fileName,'Pos1') %&& contains(filePath,'f100_n10_14_nC5_PRF') % %&&  contains(fileName,'Pos2')&& contains(fileName,'Ac4_43') %&&  % 

            scatSum= scatSum_all{ii,1};
%             if exist('axL','var')
                aimg.axL=axL;
                aimgIni.axL=axL;
%             end
            if flagDisp==1
                aimg.interpFac=interpFac;
                aimg.parLoopNum=coreNumDisp;
                
                scatSum=scatSum(:,1:dwS:end,:);
                aimg1=aimg;
                aimg1.ARF_loc=aimg.ARF_loc(1:dwS:end);  %aimgIni.ARF_loc=aimgIni.ARF_loc(1:dwS:end);  
                aimg1.numRays=numel(aimg1.ARF_loc);     %aimgIni.numRays=numel(aimg.ARF_loc);
                aimg1.scatSum=reshape(permute(scatSum,[1 3 2]),[size(scatSum,1),size(scatSum,2)*size(scatSum,3)]);     

                [dispData{1,1},~, axial, lat, t] = ncorrWrapper_loupass(aimg1);
                aimgIni.swPushLoc=pushLoc(ii);
                savefast(fullfile(savePath,saveName1),'dispData','axial','lat','t','aimgIni')
                disp(['!!!! Frame: ', mat2str(ii),'/',mat2str(aimg.SWIFrames),' !!!!!'])
            elseif flagDisp==2
                aimg.rmFrame=rmFr;
                aimg.lag=lag;
                aimgN=aimg;
                PRF_interval=median(diff(aimg.track_timeCom));
                diffT=round(diff(aimg.track_timeCom)*100)/100;
                if aimg.rmFrame==1
                    indRm=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
                elseif aimg.rmFrame==2
                    indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
                    indRm_2=circshift(indRm_1,1); % 2nd frame just after push
                    indRm=logical(indRm_1.*indRm_2);
                elseif aimg.rmFrame==3
                    indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
                    indRm_2=circshift(indRm_1,1); % 2nd frame just after push
                    indRm_3=circshift(indRm_2,1); % 2nd frame just after push
                    indRm=logical(indRm_1.*indRm_2.*indRm_3);    
                elseif aimg.rmFrame==4
                    indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
                    indRm_2=circshift(indRm_1,1); % 2nd frame just after push
                    indRm_3=circshift(indRm_2,1); % 2nd frame just after push
                    indRm_4=circshift(indRm_3,1); % 2nd frame just after push
                    indRm=logical(indRm_1.*indRm_2.*indRm_3.*indRm_4);   
                 elseif aimg.rmFrame==5
                    indRm_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
                    indRm_2=circshift(indRm_1,1); % 2nd frame just after push
                    indRm_3=circshift(indRm_2,1); % 2nd frame just after push
                    indRm_4=circshift(indRm_3,1); % 2nd frame just after push
                    indRm_5=circshift(indRm_4,1); % 2nd frame just after push
                    indRm=logical(indRm_1.*indRm_2.*indRm_3.*indRm_4.*indRm_5);  
                else
                    indRm=true(size(scatSum,4),1);
                end
%                 if aimg.rmFrame>0
%                     t_1=aimg.track_timeCom(indRm);
%                     scatSum_1=permute(scatSum(:,:,:,indRm),[4 1 2 3]);
%                     scatSum=permute(interp1(t_1,scatSum_1,aimg.track_timeCom,'spline'),[2 3 4 1]);
%                     clear scatSum_1
%                 end
                if aimgN.rmFrame>0
                        aimgN.track_time=aimg.track_time(indRm);
                        aimgN.track_timeCom=aimg.track_timeCom(indRm);
                        scatSum=scatSum(:,:,indRm);
                        aimgN.numTrackLines=sum(indRm);
                end
                aimgN.scatSum=reshape(permute(scatSum,[1 3 2]),[size(scatSum,1),size(scatSum,2)*size(scatSum,3)]);    
                aimgN.interpFac=interpFac;
                aimgN.parLoopNum=coreNumDisp;
                 [dispData{1,1}, ~, axial, lat, t] = ncorrWrapper_loupass_iF(aimgN);
                  aimgIni.swPushLoc=pushLoc(ii);
                savefast(fullfile(savePath,saveName1),'dispData','axial','lat','t','aimgIni')
                disp(['!!!! Frame: ', mat2str(ii),'/',mat2str(aimg.SWIFrames),' !!!!!'])
            elseif flagDisp==3
                scatSum=scatSum(:,1:dwS:end,:);
                aimg1=aimg;
                aimg1.ARF_loc=aimg.ARF_loc(1:dwS:end);  %aimgIni.ARF_loc=aimgIni.ARF_loc(1:dwS:end);  
                aimg1.numRays=numel(aimg1.ARF_loc);     %aimgIni.numRays=numel(aimg.ARF_loc);
                aimg1.scatSum=reshape(permute(scatSum,[1 3 2]),[size(scatSum,1),size(scatSum,2)*size(scatSum,3)]);     
                [dispData{1,1},CC{1,1},~,axial, lat, t] = ncorrWrapper_HARF_xCorr(aimg1);
                aimgIni.swPushLoc=pushLoc(ii);
                parSave(fullfile(savePath,saveName1),'dispData',dispData,'axial',axial,'lat',lat,'t',t,...
                            'aimgIni',aimgIni,'CC',CC)
                disp(['!!!! Frame: ', mat2str(ii),'/',mat2str(aimg.SWIFrames),' !!!!!'])
            end

            end
        end
%              parSave(fullfile(filePath,saveName),'dispData',dispData,'axial',axial,'lat',lat,'t',t,...
%             'aimgIni',aimgIni,'CC',CC)
    
end
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
  delete(gcp('nocreate'))
end
end
clc

disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
