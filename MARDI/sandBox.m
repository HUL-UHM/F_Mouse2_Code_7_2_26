clear all
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq36_Ang000_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_36.mat')
Id =1;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq37_Ang015_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_37.mat')
Id=2;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq38_Ang030_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_38.mat')
Id=3;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq39_Ang060_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_39.mat')
Id=4;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq40_Ang090_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_40.mat')
Id=5;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq41_Ang105_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_41.mat')
Id=6;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
load('F:\Mouse2\2026.03.03\39_67_Blue_2p3mm_D19\MARDI\Fun\T_76_Fs_2\Step06\Acq42_Ang-15_MARDI_03-March-2026_11-55-12\Process_Fun\filtDisp2D\xCorL41p5rF0rF0Filt_25_42.mat')
Id=7;
dd=nanmedian(avg_p2p_inte(:,:,1,1),2);
dd100(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,5),2);
dd500(:,Id)=dd;
dd=nanmedian(avg_p2p_inte(:,:,1,10),2);
dd1000(:,Id)=dd;
dd100=dd100(knnsearch(axial,10):knnsearch(axial,20),:);
dd500=dd500(knnsearch(axial,10):knnsearch(axial,20),:);
dd1000=dd1000(knnsearch(axial,10):knnsearch(axial,20),:);
%%
X = dd100;
Xc = X- mean(X,2);   % remove mean trend
[coeff,score,latent] = pca(Xc');
figure
scatter(score(:,1),score(:,2),100,'filled')
text(score(:,1),score(:,2),{'d1','d2','d3','d4','d5','d6','d7'})
xlabel('PC1')
ylabel('PC2')
title('Curve distribution in PCA space')
k = 2;   % number of PCs
X_recon = score(:,1:k)*coeff(:,1:k)';
X_recon = X_recon';
err = sqrt(mean((X - X_recon).^2));
threshold = median(err) + 0.5*std(err);

reject = err > threshold
%%
X = dd1000;
trend = median(X,2);
err = sqrt(mean((X - trend).^2));
mad_err = median(abs(err - median(err)));
threshold = median(err) + 1*mad_err;
reject = err > threshold;
figure('Position',[34   436   560   420]);plot(X)
figure('Position',[597   434   560   420]);plot(X(:,reject))
figure('Position',[1166         438         560         420]);plot(X(:,~reject))