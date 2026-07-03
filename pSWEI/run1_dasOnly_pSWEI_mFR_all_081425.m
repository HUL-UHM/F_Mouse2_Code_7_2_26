clear all
% close all
clc
addpath(genpath('D:\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='E:\Walid\2025.12.16';

fSh='pRF';
rSh='das';   
tFnumAll=[1.5];fNAll=['F1p5'];

interpFac=4;
coreNumDas=10;


for pp=1:numel(tFnumAll)
    tFnum=tFnumAll(pp);
    fN=fNAll(pp,:);

for mm=1:size(main_path_all,1)
%     if mm==2
%         buffNum=3;
%     else
        buffNum=2;
%     end
main_path=main_path_all{mm,1}
cd(main_path)
temp=dir;
temp(~[temp.isdir]) = [];  %remove non-directories
dummy_var=3;
done_already=0;

file_name_all=[];

for kk=dummy_var:length(temp)

    file_name=subdir(fullfile(main_path,temp(kk).name,['*_',fSh,'*Pos*_Ac*.mat']));
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
%%
% parpool(coreNumDas)
for bb=size(file_name_all,2):-1:1
    disp(['!!!! Start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
    [filePath,fileName]=fileparts(file_name_all{1,bb});
    if contains(rSh,'If')
        saveName=strrep(fileName,fSh,[rSh,num2str(lag),'Rm',num2str(rmFr),fN]);
    else
        saveName=strrep(fileName,fSh,[rSh,fN]);
    end
    if ~exist(fullfile(filePath,[saveName,'.mat']),'file') %&& (contains(filePath,'f100_n10_16_nC5_PRF15')...
            %|| contains(filePath,'f100_n10_17_nC5_PRF15')|| contains(filePath,'f100_n10_14_nC5_PRF15'))%&&  %&&  % 
        % &&
        %(contains(filePath,'phan3_inc1_Dia6')... % && contains(fileName,'Pos2') 
            %
        clear Receive 
        load(file_name_all{1,bb})

        tempPath=filePath;

        indTemp=strfind(fileName,'_Pos');
        tempName=strrep(fileName,['_Pos',fileName(indTemp+4)],'_Pos2');
        tempName=strrep(tempName,tempName(1:2),'50');
        if contains(tempName,'_Ac2')
             tempName=strrep(tempName,'Ac2','Ac1');
        end
        if ~exist('Receive','var')
            load(fullfile(tempPath,tempName),'Receive')
        end


         aimg.maxDisplacement =400e-6; 
        tic
        aimg.Receive=Receive;
        aimg.figuresOn = 0;
        aimg.ImgDataP = ImgDataP;
        % clearvars -EXCEPT aimg RcvData Receive Event
        if isfield(aimg,'nchan')
            nchan=aimg.nchan;
        else
            nchan=128;
        end
        aimg.trackFNum=tFnum;
        Receive=aimg.Receive;
        ARFlocs=aimg.ARF_loc;
        numpushlocs=length(ARFlocs);
        numRays = aimg.numRays;
        
        % make nAng for displacement tracking
        aimg.nAng=aimg.nAngSW;
        aimg.angles=aimg.anglesSW;
        nAng=aimg.nAng;
%         lengthseq=aimg.numTrackLines*nAng;

         lengthseq=numel(aimg.track_time);
            

        aimgIni=aimg;

        axLen=length(Receive(1).startSample:Receive(1).endSample);
        len_3rd=length(ARFlocs)*length(aimg.track_time);
        
        tic
        ind=find( cellfun(@(x)isequal(x,buffNum),{Receive.bufnum}),1,'first'); 
        dispData=cell(aimg.SWIFrames,1);
        RF=cell(aimg.SWIFrames,1);
        CC=cell(aimg.SWIFrames,1);
        scatSum=cell(aimg.SWIFrames,1);
        for ii=1:aimg.SWIFrames
           
            temp=unSWEI_data(Receive(ind).startSample:Receive(ind+lengthseq-1).endSample,1:nchan,ii);
            temp=permute(reshape(temp,[axLen,nAng,length(aimg.track_time)/nAng,nchan]),[1 4 2 3]);
            aimg.data=single(temp);        
            % chagne later if needed  %%%%%%%%%%%%%%%%%%%%%%%%%%
    %         aimg.frequencyMHz=6.097;
            aimg.bpf=1;  aimg.bpfVal=[0 0];    
            scatSum{ii,1}=squeeze(nanmean(DAS_SWEI_1202(aimg,coreNumDas),3));
            
%         if flagDisp==3
%         else
%             savefast(fullfile(filePath,saveName),'dispData','RF','axial','lat','t','aimgIni')
%         end
         toc
         close all
         ii
%         clear aimg
        end
        savefast(fullfile(filePath,saveName),'scatSum','aimgIni')

    end
    
end
  delete(gcp('nocreate'))
    disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!'])
    main_path_all{mm,1}
end
end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
run2_dispOnly_pSWEI_mFR_all_081425
run3_gvOnly_pSWEI_mFR_all_081425
dasBmode_pSWEI_mFR_all_081425