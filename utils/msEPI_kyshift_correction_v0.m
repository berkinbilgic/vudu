function [Kcorrected] = msEPI_kyshift_correction_v0(K, prot, kyshift, phs_sign)

% correction of kyshift value due to rotation of principle gradient axis

if nargin < 4
    phs_sign = 1;
end


    % When the value for these coordinate is very small, e.g.
    % sPosition_dSag = 1.343452e-9 then the read_meas_dat would not
    % recognize it and will leave the array empty so fix it here
    if isempty(prot.sSliceArray(1).sPosition_dSag) && isempty(prot.sSliceArray(2).sPosition_dSag)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dSag = 0;
            prot.sSliceArray(count).sNormal_dSag = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dCor) && isempty(prot.sSliceArray(2).sPosition_dCor)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dCor = 0;
            prot.sSliceArray(count).sNormal_dCor = 0;
        end
    end
    if isempty(prot.sSliceArray(1).sPosition_dTra) && isempty(prot.sSliceArray(2).sPosition_dTra)
        for count = 1:length(prot.sSliceArray(1:end))
            prot.sSliceArray(count).sPosition_dTra = 0;
            prot.sSliceArray(count).sNormal_dTra = 0;
        end
    end
    
    NormalVec = [prot.sSliceArray(1).sNormal_dSag, prot.sSliceArray(1).sNormal_dCor, prot.sSliceArray(1).sNormal_dTra].';   
    Pos(:,1) = [prot.sSliceArray(1:end).sPosition_dSag].';
    Pos(:,2) = [prot.sSliceArray(1:end).sPosition_dCor].';
    Pos(:,3) = [prot.sSliceArray(1:end).sPosition_dTra].';
    SlicePos = Pos*NormalVec;
    yPos=sqrt((Pos(:,1).^2+Pos(:,2).^2+Pos(:,3).^2)-SlicePos.^2);
    delta_ky=1/prot.sSliceArray(1).dPhaseFOV;   
    global_phase=exp(-1i*2*pi*(kyshift.*delta_ky).*yPos(1) * phs_sign);
    Kcorrected=K.*global_phase;
    
%     for ii=1:size(yPos,1)      
%         Kcorrected(:,:,:,ii)  = Kcorrected(:,:,:,ii).*exp(-1i*global_phase(ii));
%     end
    
end
