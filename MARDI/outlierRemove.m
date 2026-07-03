function xFilter=outlierRemove(X,loopDim,axIndCheck)
% assume X can be maximum 3D
% loopDim: Dimesion index for looping
    xFilter = X;
    for ii=1:size(X,loopDim)
        if loopDim==2
           Xsh=squeeze(X(axIndCheck,ii,:));
        elseif loopDim==3
           Xsh=squeeze(X(axIndCheck,:,ii));
        end
        trend = median(Xsh,2);
        err = sqrt(mean((Xsh - trend).^2));
        mad_err = median(abs(err - median(err)));
        threshold = median(err) + 1*mad_err;
        reject = err > threshold;
        % Xsh(:,reject)=nan;
        if loopDim==2
           xFilter(:,ii,reject)=nan;
        elseif loopDim==3
           xFilter(:,reject,ii)=nan;
        end
    end
end