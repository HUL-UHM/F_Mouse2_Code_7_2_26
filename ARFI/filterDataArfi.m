function [arfidataFilt,tNew]=filterDataArfi(t,arfidata,aimgIni,axial)
FD=aimgIni.focalDepthmm;
PRF_interval=median(diff(t));
diffT=round(diff(t)*100)/100;
if  aimgIni.rmFrDisp==1
    ind=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
elseif  aimgIni.rmFrDisp==2
    ind_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
    ind_2=circshift(ind_1,1); % 2nd frame just after push
    ind=logical(ind_1.*ind_2);
elseif  aimgIni.rmFrDisp==3
    ind_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
    ind_2=circshift(ind_1,1); % 2nd frame just after push
    ind_3=circshift(ind_2,1); % 2nd frame just after push
    ind=logical(ind_1.*ind_2.*ind_3);
elseif  aimgIni.rmFrDisp==4
    ind_1=[true ~(diffT~=(round((100*PRF_interval)))/100)]'; % find frame just after push
    ind_2=circshift(ind_1,1); % 2nd frame just after push
    ind_3=circshift(ind_2,1); % 2nd frame just after push
    ind_4=circshift(ind_3,1); % 2nd frame just after push
    ind=logical(ind_1.*ind_2.*ind_3.*ind_4);
else
    ind=1:length(t);
end
t_1=t(ind);
arfidata_1=arfidata(:,:,ind);
[b,a]=butter(2,2*[100 2500]*median(diff(t_1*1e-3)),'bandpass'); 
arfidata_1=permute(filtfilt(b,a,permute(arfidata_1,[3 1 2])),[2 3 1]);
dx=median(diff(axial));
avgSam=round(1/dx);
arfidata_1=movmean(arfidata_1,avgSam);% axial direction
%% noise estimation
tTh=round(aimgIni.numPreFrames-3);
meanNoise=median(arfidata_1(:,:,1:aimgIni.numPreFrames),3);
meanNoiseRep=repmat(meanNoise,[1 1 size(arfidata_1,3)]);
arfidata_1=arfidata_1-meanNoiseRep;
% arfidata_1=abs(arfidata_1);

t_1=t_1-t_1(tTh);
% t_1=t_1(1:end);
% arfidata_1=arfidata_1(:,:,tTh:end);

t_interp_1=t_1(1):1/aimgIni.PRFKHz:t_1(end);

%%
FD_ind=knnsearch(axial,FD);
lat_id_1=2;
dispThresh=10;
disp_temp=(squeeze(arfidata_1(FD_ind,lat_id_1,:)));

ind_1=[true;~(abs(diff(disp_temp))>dispThresh)]; % find interframe displacement > 1 micron
if all(ind_1)
%             arfidata_1=interp1(t_1,permute(arfidata_1,[3 1 2]),t_interp_1,'spline');
    try arfidata_1=interp1(t_1,permute(arfidata_1,[3 1 2]),t_interp_1,'spline');
    catch ME        
        for ii=1:size(arfidata_1,1)
            for jj=1:size(arfidata_1,2)
                if all(isnan(squeeze(arfidata_1(ii,jj,:))))==1
                 arfidata_1(ii,jj,[1 size(arfidata_1,3)])=[0 1];
                end
            end
        end
        arfidata_1=interp1(t_1,permute(arfidata_1,[3 1 2]),t_interp_1,'spline');
    end
else
    t_1=t_1(ind_1);
    arfidata_1=arfidata(:,:,ind_1);
    try arfidata_1=interp1(t_1,permute(arfidata_1,[3 1 2]),t_interp_1,'spline');
    catch ME

        for ii=1:size(arfidata_1,1)
            for jj=1:size(arfidata_1,2)
                if all(isnan(squeeze(arfidata_1(ii,jj,:))))==1
                 arfidata_1(ii,jj,[1 round(size(arfidata_1,3)/2) size(arfidata_1,3)])=[0 0.5 1];
                end
            end
        end
        arfidata_1=interp1(t_1,permute(arfidata_1,[3 1 2]),t_interp_1,'spline');
    end
end
%%
arfidataFilt =permute(arfidata_1,[2 3 1]);
tNew=t_interp_1;
