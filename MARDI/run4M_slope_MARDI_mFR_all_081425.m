clear all
% close all
clc
addpath(genpath('D:\Code\functions'))
addpath(genpath('F:\Mouse2\Code'))
ind =0;
ind=ind+1;main_path_all{ind,1}='F:\mouse3';
parLoopNum =8;
Ds=2;
windS=[1.5 1.5];
fSh='CorL41p5rF0rF0Filt_';
rSh='CorL41p5rF0rF0FiltLinFit_';
nCycle=1;


for mm=1:size(main_path_all,1)
    main_path=main_path_all{mm,1};
    cd(main_path)
    temp=dir;
    temp(~[temp.isdir]) = [];  %remove non-directories
    dummy_var=3;
    done_already=0;

    file_name_all=[];

    for kk=dummy_var:length(temp)

        file_name=subdir(fullfile(main_path,temp(kk).name,['*',fSh,'*.mat']));
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

    for bb=1:size(file_name_all,2)
         disp(['!!!! start : ', mat2str(bb),'/',mat2str(size(file_name_all,2)),' !!!!!'])
        [filePath,fileName]=fileparts(file_name_all{1,bb});     
        saveName=strrep(fileName,fSh,rSh);
        savePath=filePath;
        if ~exist(savePath,'dir')
            mkdir(savePath)
        end
        comN=filePath;
        if ~exist(fullfile(savePath,[saveName,'.mat']),'file') 
            % if  contains(filePath,'HP_4_6_5')
                TCName='D:\Code\ProcessData_Rat\MARDI\TC_HMI_L115\f100_n10_nC4_PRF10_PI\f1Fac_2_HP_4_6_5_nC_4.mat';
            % end
            load(TCName, 'fftAmpFreq')
            norm1=fftAmpFreq.^1./sum(fftAmpFreq.^1);
            norm1=norm1(1)./norm1;
            
            load(file_name_all{1,bb},'avg_p2p_inte','lat','axial','factor','passParams');
            factor.windS=ceil(0.5*windS./median(diff(lat)));
            wSax=factor.windS(1);
            wSlt=factor.windS(2);
            factor.downSamplePhaseCalc=Ds;
            factor.core=parLoopNum;
%             freqAll=factor.harmFreq;

            DispVal=squeeze(avg_p2p_inte(:,:,nCycle,:));
            fVal=squeeze(passParams.centerFreq(:,:,nCycle,:));
            [M,N,T]=size(DispVal);spaTotI=M*N; 
            downSamp=factor.downSamplePhaseCalc;
            indM1=1:downSamp:M;mTot=length(indM1);
            indN1=1:downSamp:N;nTot=length(indN1);
            spaTot=mTot*nTot;
%             disp1=reshape(DispVal(indM1,indN1,:),spaTot,T);
%             fVal1=reshape(fVal(indM1,indN1,:),spaTot,T);
            disp1=reshape(DispVal,spaTotI,T);
            fVal1=reshape(fVal,spaTotI,T);
            
            Lin1ori=nan(spaTot,2);Lin1norm=nan(spaTot,2);
            Lin1R2ori=nan(spaTot,1);Lin1R2norm=nan(spaTot,1);
            Lin1RMSEori=nan(spaTot,1);Lin1RMSEnorm=nan(spaTot,1);
            
            Lin2ori=nan(spaTot,3);Lin2norm=nan(spaTot,3);
            Lin2R2ori=nan(spaTot,1);Lin2R2norm=nan(spaTot,1);
            Lin2RMSEori=nan(spaTot,1);Lin2RMSEnorm=nan(spaTot,1);
                         parpool(factor.core);
            % D = parallel.pool.Constant(disp1);
            D.Value=disp1;
            for ii=1:spaTot
                [rIt,cIt]=ind2sub([mTot,nTot],ii);
                indM=1:downSamp:M;
                indN=1:downSamp:N;
                rI=indM(rIt);
                cI=indN(cIt);
                axLen=max(1,rI-wSax):min(M,rI+wSax);  
                latLen=max(1,cI-wSlt):min(N,cI+wSlt);
                [lI,aI]=meshgrid(latLen,axLen);
                iiAll=sub2ind([M,N],aI(:),lI(:));
                
                freqAll=fVal1(iiAll,:);
                dd=abs(D.Value(iiAll,:));
                ddN=repmat(dd(:,1),1,size(dd,2))./dd;
                ddN=ddN(:);
                if sum(isnan(dd))<5
                    p=polyfitn(freqAll(:),log10(ddN),1);   
                    Lin1R2ori(ii,1)=p.AdjustedR2;Lin1RMSEori(ii,1)=p.RMSE;
                    Lin1ori(ii,:)=p.Coefficients';
                    
                    p=polyfitn(freqAll(:),log10(ddN),2);   
                    Lin2R2ori(ii,1)=p.AdjustedR2;Lin2RMSEori(ii,1)=p.RMSE;
                    Lin2ori(ii,:)=p.Coefficients';
                end
                dd=abs(D.Value(iiAll,:)).*repmat(norm1,numel(iiAll),1);
                ddN=repmat(dd(:,1),1,size(dd,2))./dd;
                ddN=ddN(:);
                if sum(isnan(dd))<5
                    p=polyfitn(freqAll(:),log10(ddN),1);   
                    Lin1R2norm(ii,1)=p.AdjustedR2;Lin1RMSEnorm(ii,1)=p.RMSE;
                    Lin1norm(ii,:)=p.Coefficients';
                    
                    p=polyfitn(freqAll(:),log10(ddN),2);   
                    Lin2R2norm(ii,1)=p.AdjustedR2;Lin2RMSEnorm(ii,1)=p.RMSE;
                    Lin2norm(ii,:)=p.Coefficients';
                end
            end
            delete(gcp('nocreate'));
            
            [nGdI,mGdI]=meshgrid(1:N,1:M);
            [nGd,mGd]=meshgrid(indN1,indM1);
            
            GOF.Lin1R2ori=interp2(nGd,mGd,reshape(Lin1R2ori,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin1RMSEori=interp2(nGd,mGd,reshape(Lin1RMSEori,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin1R2norm=interp2(nGd,mGd,reshape(Lin1R2norm,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin1RMSEnorm=interp2(nGd,mGd,reshape(Lin1RMSEnorm,mTot,nTot),nGdI,mGdI,'cubic');
            
            GOF.Lin2R2ori=interp2(nGd,mGd,reshape(Lin2R2ori,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin2RMSEori=interp2(nGd,mGd,reshape(Lin2RMSEori,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin2R2norm=interp2(nGd,mGd,reshape(Lin2R2norm,mTot,nTot),nGdI,mGdI,'cubic');
            GOF.Lin2RMSEnorm=interp2(nGd,mGd,reshape(Lin2RMSEnorm,mTot,nTot),nGdI,mGdI,'cubic');
            
            Param.Lin1ori=nan(M,N,2);Param.Lin1norm=nan(M,N,2);
            Param.Lin2ori=nan(M,N,3);Param.Lin2norm=nan(M,N,3);
            for jj=1:3  
                if jj<3
                    Param.Lin1ori(:,:,jj)=interp2(nGd,mGd,reshape(Lin1ori(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                    Param.Lin1norm(:,:,jj)=interp2(nGd,mGd,reshape(Lin1norm(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                    Param.Lin2ori(:,:,jj)=interp2(nGd,mGd,reshape(Lin2ori(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                    Param.Lin2norm(:,:,jj)=interp2(nGd,mGd,reshape(Lin2norm(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                else
                    Param.Lin2ori(:,:,jj)=interp2(nGd,mGd,reshape(Lin2ori(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                    Param.Lin2norm(:,:,jj)=interp2(nGd,mGd,reshape(Lin2norm(:,jj),mTot,nTot),nGdI,mGdI,'cubic');
                end
            end          
            
           savefast(fullfile(savePath,saveName),'GOF','Param','factor','axial','lat');
        end
       
    end

        disp(['!!!! Done : ', mat2str(mm),'/',mat2str(size(main_path_all,1)),' !!!!!']);

end
clc
disp('!!!!!!!!  YOU CAN TURN OFF IF YOU SEE THIS MESSAGE !!!!!!!!!!!!!!!')
