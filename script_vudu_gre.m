%--------------------------------------------------------------------------
%% load undersampled EPI data at R=3, AP encoding in shot1, PA in shot2
%% (after ghost correction, and removal of readout oversampling)
%% 2021 07 01
%--------------------------------------------------------------------------

addpath utils/

load data/Img_epi_ap.mat
load data/Img_epi_pa.mat
load data/receive.mat       % espirit coil sensitivity estimated from reference flash data   

mosaic( rsos(Img_epi_ap(:,:,:,:),3), 4, 6, 1, 'ap @ R=3', [0,1e-3], 90 )
mosaic( rsos(Img_epi_pa(:,:,:,:),3), 4, 6, 2, 'pa @ R=3', [0,1e-3], 90 )


% signal norm ratio:
norm2(Img_epi_pa) / norm2(Img_epi_ap)   


%--------------------------------------------------------------------------
%% forward model for subsampled data
%-------------------------------------------------------------------------- 

[N(1), N(2), num_chan, num_slice] = size(Img_epi_ap);

Ry = 3;             % accl factor
del_ky = [0,1];     % delta_ky shift between shots

A = fftc(eye(N(2)),1);          

A_ap = A(1+del_ky(1):Ry:end,:);
A_pa = A(1+del_ky(2):Ry:end,:);


% (x,y,chan,ap/pa,slc,dwi)
img_use = cat(4, permute(Img_epi_ap, [1,2,3,6,4,5]), permute(Img_epi_pa, [1,2,3,6,4,5]));


% (x, ky, chan) subsampled data
sgnl_ap = zeross([N ./ [1,Ry], num_chan, num_slice]);
sgnl_pa = zeross([N ./ [1,Ry], num_chan, num_slice]);

tic
for nslc = 1:num_slice
    for xn = 1:N(1)
        m_ap = sq( img_use(xn,:,:,1,nslc) );
        m_pa = sq( img_use(xn,:,:,2,nslc) );

        % signal in x,ky,chan
        sgnl_ap(xn, :, :, nslc) = A_ap  * sq(m_ap);
        sgnl_pa(xn, :, :, nslc) = A_pa  * sq(m_pa);    
    end
end
toc

  
%--------------------------------------------------------------------------
%% Sense for AP and PA separately without B0 model 
%--------------------------------------------------------------------------

sens = permute(receive, [1,2,4,3]);

lambda_tik = 1e-6;      % tikhonov regularization parameter

img_sense = zeross([N, 2, num_slice]);

tic
% recon middle slice
for nslc = 16%1:num_slice
    disp(['slice: ', num2str(nslc), ' / ', num2str(num_slice)])
    
    for xn = 1:N(1)
        AC_ap = zeross([num_chan * N(2) / Ry, N(2)]);
        AC_pa = zeross([num_chan * N(2) / Ry, N(2)]);

        for c = 1:num_chan       
            AC_ap(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :) = A_ap * diag( sens(xn,:,nslc,c) );
            AC_pa(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :) = A_pa * diag( sens(xn,:,nslc,c) );
        end

        [U_ap, S_ap, V_ap] = svd(AC_ap, 'econ');
        [U_pa, S_pa, V_pa] = svd(AC_pa, 'econ');    
    
        Einv_ap = V_ap * diag(diag(S_ap) ./ (diag(S_ap).^2 + lambda_tik)) * U_ap';
        Einv_pa = V_pa * diag(diag(S_pa) ./ (diag(S_pa).^2 + lambda_tik)) * U_pa';

        rhs_ap = sgnl_ap(xn, :, :, nslc);   
        rhs_ap = rhs_ap(:);
        img_sense(xn,:,1, nslc) = Einv_ap * rhs_ap;

        rhs_pa = sgnl_pa(xn, :, :, nslc);   
        rhs_pa = rhs_pa(:);
        img_sense(xn,:,2, nslc) = Einv_pa * rhs_pa;
    end
end
toc

img_sense = apodize(img_sense);     % perform k-space apodization

nslc = 16;  % slices to display

mosaic(img_sense(:,:,1,nslc),1,1,11, 'sense ap', [0,1.5e-3],90)
mosaic(img_sense(:,:,2,nslc),1,1,12, 'sense pa', [0,1.5e-3],90)

mosaic(fft2call(img_sense(:,:,1,nslc)).^.5,1,1,1, 'ap k-space', [0,5e-2],0), setGcf(.5)
mosaic(fft2call(img_sense(:,:,2,nslc)).^.5,1,1,2, 'pa k-space', [0,5e-2],0), setGcf(.5)



%--------------------------------------------------------------------------
%% Vudu joint recon with Hankel low-rank constraint
%--------------------------------------------------------------------------

load data/img_fieldmap.mat      % field map estimated using Topup

scale_fleet = norm2(img_sense(:,:,2,:)) / norm2(img_sense(:,:,1,:));    % scale factor between the two shots


esp = 3.3e-4;   % echo spacing in ms

t_axis_ap = [0:Ry:N(2)-1] * esp;
t_axis_pa = t_axis_ap(end:-1:1);

winSize = [1,1] * 11;
step_size = 0.5;

num_shot = 2;
use_scaling = 1;

vc = 0;
lambda_msl = (1+vc);

num_iter = 100;
tol = 0.1;

keep = 1:floor(lambda_msl*prod(winSize));

img_vudu = zeross([N,2,num_slice]);


tic
for nslc = 16%1:num_slice
    disp(['slice: ', num2str(nslc)])
    
    % create and store encoding matrices
    AWC_ap = zeross([N(2)*num_chan/Ry, N(2), N(1)]);
    AWC_pa = zeross([N(2)*num_chan/Ry, N(2), N(1)]);

    AWC_apH = zeross([N(2), N(2)*num_chan/Ry, N(1)]);
    AWC_paH = zeross([N(2), N(2)*num_chan/Ry, N(1)]);

    AWC_apN = zeros(size(AWC_ap, 2), size(AWC_ap, 2), size(AWC_ap, 3));
    AWC_paN = zeros(size(AWC_pa, 2), size(AWC_pa, 2), size(AWC_pa, 3));

    AWC_apHrhs =  zeros(size(AWC_apH, 1), size(AWC_apH, 3));
    AWC_paHrhs =  zeros(size(AWC_paH, 1), size(AWC_paH, 3));

    for xn = 1:N(1)
        b0 = img_fieldmap(xn,:,nslc) * 2 * pi;

        W_ap = exp(1i * t_axis_ap.' * b0);
        AW_ap = A_ap .* W_ap;      

        W_pa = exp(1i * t_axis_pa.' * b0);
        AW_pa = A_pa .* W_pa;

        for c = 1:num_chan       
            AWC_ap(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :, xn) = AW_ap * diag( sens(xn,:,nslc,c) );
            AWC_pa(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :, xn) = AW_pa * diag( sens(xn,:,nslc,c) );
        end

        AWC_apH(:,:,xn) = AWC_ap(:,:,xn)';
        AWC_paH(:,:,xn) = AWC_pa(:,:,xn)';

        AWC_apN(:, :, xn) = AWC_apH(:,:,xn) * AWC_ap(:,:,xn);
        AWC_paN(:, :, xn) = AWC_paH(:,:,xn) * AWC_pa(:,:,xn);

        rhs_ap = sgnl_ap(xn, :, :, nslc);
        AWC_apHrhs(:,xn) = AWC_apH(:,:,xn) * rhs_ap(:);
        
        rhs_pa = sgnl_pa(xn, :, :, nslc);
        
        if use_scaling
            rhs_pa = rhs_pa / scale_fleet;
        end
        
        AWC_paHrhs(:,xn) = AWC_paH(:,:,xn) * rhs_pa(:);
    end    
    
    % initialize with sense solution: 
    img_sense_ap = sq(img_sense(:,:,1,:,:));
    img_sense_pa = sq(img_sense(:,:,2,:,:));

    im_rec = permute(cat(3, img_sense_ap(:,:,nslc), img_sense_pa(:,:,nslc)), [2,3,1]);

    for iter = 1:num_iter
        im_prev = im_rec;

        for xn = 1:N(1)
            im_rec(:,1,xn) = im_rec(:,1,xn) - step_size * ( AWC_apN(:,:,xn) * im_rec(:,1,xn) - AWC_apHrhs(:, xn) );
            im_rec(:,2,xn) = im_rec(:,2,xn) - step_size * ( AWC_paN(:,:,xn) * im_rec(:,2,xn) - AWC_paHrhs(:, xn) );
        end

        im_rec = permute(im_rec, [3,1,2]);

        if vc
            im_rec = cat(3,im_rec, conj(im_rec));
        end

        A = Im2row( fft2call(im_rec), winSize );

        [U, S, V] = svd(A, 'econ');

        A = U(:,keep) * S(keep,keep) * V(:,keep)';

        k_pocs = Row2im(A, [N, num_shot*(vc+1)], winSize);

        im_rec_disp = ifft2call(k_pocs);

        im_rec_disp = im_rec_disp(:,:,1:num_shot);

        im_rec = permute(im_rec_disp, [2,3,1]);  

        update = rmse(im_prev,im_rec);
        
        if ~mod(iter,10)
            mosaic(im_rec_disp, 1, 2, 100+vc, ['iter: ', num2str(iter), '  update: ', num2str(rmse(im_prev,im_rec))], genCaxis(im_rec), 90)
        end
        
        if update < tol
            break
        end
    end

    img_vudu(:,:,:,nslc) = permute(im_rec, [3,1,2]);
end
toc

img_vudu = apodize(img_vudu);     % perform k-space apodization


mosaic(mean(abs(img_vudu(:,:,1,nslc,1)),3), 1, 1, 21, 'vudu ap', [0,1.5e-3], 90)
mosaic(mean(abs(img_vudu(:,:,2,nslc,1)),3), 1, 1, 22, 'vudu pa', [0,1.5e-3], 90)



