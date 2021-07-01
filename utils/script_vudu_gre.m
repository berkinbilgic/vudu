%--------------------------------------------------------------------------
%% load data
% dwi:     standard buda R3
%          fleet buda R3
%
% gre:     fleet buda R3
%
% looks like fleet msepi acquisition with ap-ap phase encoding overwrote
% the eco index, so only the second shot was recorded. need to fix this in the sequence 
% -> don't use m_bPEReverse flag to determine the eco to write
%--------------------------------------------------------------------------

addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/fleet_buda/mrir_toolbox
addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/LIBRARY/read_meas_dat__20140924112147
addpath /autofs/cluster/kawin/berkin/Matlab_Code_New/fleet_buda


data_path = '/autofs/space/marduk_001/users/berkin/2020_12_11_bay1_symri_buda_fleet_invivo/';
    
load([data_path, 'receive'], 'receive')
    

%--------------------------------------------------------------------------
% load buda fleet data
%--------------------------------------------------------------------------

% file_name = 'meas_MID00437_FID52932_a_ep2d_buda_fleet_dwi_R3_del10_b1k.dat';	
file_name = 'meas_MID00447_FID52942_a_ep2d_buda_fleet_gre_R3_del10.dat';	


dat = mapVBVD([data_path, file_name]);

meas.prot = read_meas_prot([data_path, file_name]);



%--------------------------------------------------------------------------
% select dwi direction
%--------------------------------------------------------------------------

shot_start = 1;     %b=0
% shot_start = 3;     %b=1k
num_shot = 1;           

select_shot = shot_start:shot_start+num_shot-1      
            
meas.data = dat{end}.image(:,:,:,:,:,:,:,:,select_shot,:,:,:,:,:,:,:);                     % data
meas.data = permute(meas.data, [1 3 2 10 8 7 9 11 4 5 6]);
  

% parse navigator
meas.data_phascor1d = dat{end}.phasecor(:,:,:,:,:,:,:,:,select_shot,:,:,:,:,:,:,:);        % navigator for data
meas.data_phascor1d = permute(meas.data_phascor1d, [1 3 2 10 8 7 9 11 4 5 6]);

sz = size(meas.data_phascor1d);
sz(2) = 3; 
sz(11) = 1;

test = zeros(sz);
test(:,1,:,:,:,:,:,1,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,1,:,:,1);
test(:,2,:,:,:,:,:,2,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,2,:,:,1);
test(:,3,:,:,:,:,:,1,:,:) = meas.data_phascor1d(:,:,:,:,:,:,:,1,:,:,2);
meas.data_phascor1d = test;


% deinterleave
meas.data = mrir_image_slice_deinterleave(meas.data); 
meas.data_phascor1d = mrir_image_slice_deinterleave(meas.data_phascor1d);
 

%--------------------------------------------------------------------------
% correct for delta-ky shift
%--------------------------------------------------------------------------


meas.data = meas.data(:,2:end,:,:,:,:,:,:,:,:);
meas.data(:,:,:,:,2,:,:,:,:,:) = circshift( meas.data(:,:,:,:,2,:,:,:,:,:), [0,-1] );


temp = sum(sq(meas.data), 5);

mosaic( rsos(temp(1:2:end,:,:,1,16),3).^.25,1,1,1), setGcf(.5)
mosaic( rsos(temp(1:2:end,:,:,2,16),3).^.25,1,1,2), setGcf(.5)



%--------------------------------------------------------------------------
% ghost correction
%--------------------------------------------------------------------------

startLine = 1;

AccY = meas.prot.lAccelFactPE;              % PAT factor

num_dwi = 1;%s(meas.data, 7) / 2               % number of volumes

nco_polarity = [-1, +1];                    % use -1 for ap, +1 for pa (on prisma)

del_ky = [0, 0];   % there is actually [0,1] shift, taken care of later

disp(['Pat Factor: ', num2str(AccY)])
    
size_data = s( meas.data(:,:,:,:,:,:,1:num_shot,:,:,:) );   % (dim5) -> ap/pa    (dim8) -> even/odd


Img_epi_ap = zeross( [size_data(1:3), size_data(10), num_dwi] );
Img_epi_pa = zeross( [size_data(1:3), size_data(10), num_dwi] );


tic
for ndwi = 1:num_dwi
    % load EPI-data 
    EPI_data_raw = meas.data(:,:,:,:,:,:,1+num_shot*(ndwi-1):num_shot*ndwi,:,:,:);

    % load EPI-nav
    EPI_nav = meas.data_phascor1d(:,:,:,:,:,:,1+num_shot*(ndwi-1):num_shot*ndwi,:,:,:);
    

    % Correct the global phase due to NCO
    EPI_data = zeross(size(EPI_data_raw));

    for ii = 1:2
       EPI_data(:,:,:,1,ii,1,1,:,1,:) = msEPI_kyshift_correction_v0(EPI_data_raw(:,:,:,1,ii,1,1,:,1,:), meas.prot, del_ky(ii), nco_polarity(ii)) ;
    end


    s_data = size_data;
    s_data(8) = 1;

    EPI_data_comb = zeross(size_data);      % even/odd separate
    EPI_data_cor_kspace = zeross(s_data);   % even/odd combined
    
    % shift the ky lines
    for sh = 1:2
        disp(['shot: ', num2str(sh)])
        
        ky_idx = 1 + del_ky(sh) : AccY : s(EPI_data, 2)  + del_ky(sh);
        ky_idx(ky_idx >= size_data(2)) = [] 

        EPI_data_comb(:, ky_idx,:,:,sh,:,:,:,:,:) = EPI_data(:,startLine:AccY:startLine+AccY*(length(ky_idx)-1),:,:,sh,:,:,:,:,:);

        % (ghost correct), grid, apodize EPI-data
        EPI_data_cor_kspace(:,:,:,:,sh,:,:,:,:,:) = ghost_correct_pat_ref_v1_STD_bb(meas.prot, EPI_data_comb(:,:,:,:,sh,:,:,:,:,:), EPI_nav, AccY);
    end


    Img_epi_ap(:,:,:,:,ndwi) = ifft2call(EPI_data_cor_kspace(:,:,:,1,1,1,1,1,1,:));
    Img_epi_pa(:,:,:,:,ndwi) = ifft2call(EPI_data_cor_kspace(:,:,:,1,2,1,1,1,1,:));
    
    
    mosaic(sq(rsos(EPI_data_cor_kspace(1:4:end,:,:,1,:,1,1,1,1,end/2), 3)), 2, 1, 11, 'acquired', [0,.25e-3]), setGcf(.5)

    mosaic( rsos(Img_epi_ap(:,:,:,:,ndwi),3), 4, 6, 1, 'ap', [0,.25e-3], 90 )
    mosaic( rsos(Img_epi_pa(:,:,:,:,ndwi),3), 4, 6, 2, 'pa', [0,.25e-3], 90 )
end
toc



%--------------------------------------------------------------------------
%% load undersampled EPI data at R=3, AP encoding in shot1, PA in shot2
%% (after ghost correction, and removal of readout oversampling)
%--------------------------------------------------------------------------

load data/Img_epi_ap
load data/Img_epi_pa

mosaic( rsos(Img_epi_ap(:,:,:,:,ndwi),3), 4, 6, 1, 'ap', [0,1e-3], 90 )
mosaic( rsos(Img_epi_pa(:,:,:,:,ndwi),3), 4, 6, 2, 'pa', [0,1e-3], 90 )


% signal norm ratio:
norm2(Img_epi_pa) / norm2(Img_epi_ap)   


%--------------------------------------------------------------------------
%% forward model for subsampled data
%--------------------------------------------------------------------------

mosaic( rsos(fft2call(Img_epi_ap(:,:,:,:,ndwi)),3), 1, 1, 1, 'ap', [0,.25e-3], 0 ), setGcf(.5)
mosaic( rsos(fft2call(Img_epi_pa(:,:,:,:,ndwi)),3), 1, 1, 2, 'pa', [0,.25e-3], 0 ), setGcf(.5)
 

[N(1), N(2), num_chan, num_slice] = size(Img_epi_ap);

Ry = 3;
del_ky = [0,1];     % somehow delta_ky shift did not take effect

A = fftc(eye(N(2)),1);          

A_ap = A(1+del_ky(1):Ry:end,:);
A_pa = A(1+del_ky(2):Ry:end,:);


% (x,y,chan,ap/pa,slc,dwi)
img_patref_use = cat(4, permute(Img_epi_ap, [1,2,3,6,4,5]), permute(Img_epi_pa, [1,2,3,6,4,5]));


% (x, ky, chan) subsampled data
sgnl_ap = zeross([N ./ [1,Ry], num_chan, num_slice, num_dwi]);
sgnl_pa = zeross([N ./ [1,Ry], num_chan, num_slice, num_dwi]);

tic
for ndwi = 1:num_dwi
    for nslc = 1:num_slice
        for xn = 1:N(1)
            m_ap = sq( img_patref_use(xn,:,:,1,nslc,ndwi) );
            m_pa = sq( img_patref_use(xn,:,:,2,nslc,ndwi) );

            % signal in x,ky,chan
            sgnl_ap(xn, :, :, nslc, ndwi) = A_ap  * sq(m_ap);
            sgnl_pa(xn, :, :, nslc, ndwi) = A_pa  * sq(m_pa);    
        end
    end
end
toc


  
%--------------------------------------------------------------------------
%% Sense for AP and PA separately without B0 model 
%--------------------------------------------------------------------------

load([data_path, 'receive'])

sens = permute(receive, [1,2,4,3]);


lambda_tik = 1e-6;

img_sense = zeross([N, 2, num_slice, num_dwi]);

tic
for nslc = 1:num_slice
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

        for ndwi = 1:num_dwi
            rhs_ap = sgnl_ap(xn, :, :, nslc, ndwi);   
            rhs_ap = rhs_ap(:);
            img_sense(xn,:,1, nslc, ndwi) = Einv_ap * rhs_ap;

            rhs_pa = sgnl_pa(xn, :, :, nslc, ndwi);   
            rhs_pa = rhs_pa(:);
            img_sense(xn,:,2, nslc, ndwi) = Einv_pa * rhs_pa;
        end
    end
end
toc

img_sense = apodize(img_sense);


mosaic(img_sense(:,:,1,nslc),1,1,11, '', [0,.25e-3],90)
mosaic(img_sense(:,:,2,nslc),1,1,12, '', [0,.25e-3],90)

mosaic(fft2c(img_sense(:,:,1,nslc)).^.5,1,1,1, '', [0,5e-2],0), setGcf(.5)
mosaic(fft2c(img_sense(:,:,2,nslc)).^.5,1,1,2, '', [0,5e-2],0), setGcf(.5)


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------


fleet_path = '/autofs/space/marduk_001/users/berkin/2020_12_11_bay1_symri_buda_fleet_invivo/fleet_buda_gre/';

% save([fleet_path, 'img_sense_MID447'], 'img_sense')
load([fleet_path, 'img_sense_MID447'])


mosaic(img_sense(:,:,1,1:4:end),2,4,11, '', [0,1.5e-3],90)
mosaic(img_sense(:,:,2,1:4:end),2,4,12, '', [0,1.5e-3],90)


scale_fleet = norm2(img_sense(:,:,2,:)) / norm2(img_sense(:,:,1,:))


%--------------------------------------------------------------------------
%% save as nifti
%--------------------------------------------------------------------------

% apply scaling factor
img_recon_pad = permute(img_sense, [1,2,4,3]);

img_recon_pad(:,:,:,2) = img_recon_pad(:,:,:,2) / scale_fleet;

num_slc_pad = ceil(num_slice/4)*4

img_recon_pad = padarray(img_recon_pad, [0,0, num_slc_pad - num_slice]/2);

img_recon_pad = img_recon_pad(:,2:end-1,:,:);   % divisible by 4

       
ii_dif = shot_start;
ii_rf = 1;


voxel_size = [meas.prot.dReadoutFOV, meas.prot.dPhaseFOV, meas.prot.dThickness] ./ [meas.prot.lBaseResolution, meas.prot.lPhaseEncodingLines, meas.prot.sSliceArray_lSize]

 
fpB0Topup = [fleet_path, 'img_ap_pa_',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'];

% genNii(single( abs(img_recon_pad) ), voxel_size, fpB0Topup)      
        

dt = load_untouch_nii([fleet_path, 'img_ap_pa_', num2str(ii_dif), 'rf', num2str(ii_rf), '.nii']);


imagesc3d2( dt.img(:,:,:,1), s(dt.img)/2, 31, [0,0,0]+90, [-0,1.25e-3])
imagesc3d2( dt.img(:,:,:,2), s(dt.img)/2, 32, [0,0,0]+90, [-0,1.25e-3])


%--------------------------------------------------------------------------
%% topup
%--------------------------------------------------------------------------

addpath /cluster/kawin/berkin/Matlab_Code_New/BUDA_noddi

config_path = '/autofs/cluster/kawin/Congyu/CODE/gSlider_BUDA/Recon_Data/';

save_path = fleet_path;


% create topup text file 
fpAcqp = [save_path, 'acq_param.txt'];

esp = meas.prot.iEffectiveEpiEchoSpacing * 1e-6

readout_duration = (meas.prot.lPhaseEncodingLines - 1) * esp 

acqp0 = [0 -1 0 readout_duration];
acqp1 = [0 1 0 readout_duration];
acq_par = [acqp0; acqp1];
               
        
write_acq_params(fpAcqp, acq_par)

% edit(fpAcqp)


tic
    fpB0Topup = [save_path, 'img_ap_pa_',num2str(ii_dif),'rf',num2str(ii_rf),'.nii'];

    fpOut = [save_path, 'out_ap_pa_',num2str(ii_dif),'rf',num2str(ii_rf)];
    fpIout = [save_path, 'epi_unwarp_',num2str(ii_dif),'rf',num2str(ii_rf)];
    fpField = [save_path, 'fieldmap_',num2str(ii_dif),'rf',num2str(ii_rf)];

    cmd = ['topup --imain=' fpB0Topup ' --config=', config_path, 'b02b0_topup.cnf --datain=' fpAcqp ' --out=' fpOut ' --iout=' fpIout ' --fout=' fpField]

    [status, result] = system(cmd, '-echo');   
toc



%--------------------------------------------------------------------------
%% load topup results
%--------------------------------------------------------------------------

ii_dif = shot_start;
ii_rf = 1;

save_path = fleet_path;


fpIout = [save_path, 'epi_unwarp_',num2str(ii_dif),'rf',num2str(ii_rf)];
fpField = [save_path, 'fieldmap_',num2str(ii_dif),'rf',num2str(ii_rf)];

% system(['gunzip ', save_path, '*.gz'])

dt1 = load_untouch_nii([fpField,'.nii']);       % field map
dt2 = load_untouch_nii([fpIout,'.nii']);        % distortion corrected volumes

img_fieldmap = padarray(dt1.img, [0,1]);


mosaic(img_fieldmap(:,:,1:2:end),2,5,1,'',2*[-50,50], 90)        
mosaic(dt2.img(:,:,2:3:end,1),2,5,2,'',[-0,1.25e-3], 90)        
mosaic(dt2.img(:,:,2:3:end,2),2,5,3,'',[-0,1.25e-3], 90)        


mosaic(img_recon_pad(:,:,2:3:end,1),2,5,12,'',[-0,1.25e-3], 90)        
mosaic(img_recon_pad(:,:,2:3:end,2),2,5,13,'',[-0,1.25e-3], 90)        
        

norm2(dt2.img(:,:,:,2))/norm2(dt2.img(:,:,:,1))


%--------------------------------------------------------------------------
%% apply topup on sense recon -> estimate shot-to-shot phase difference
%--------------------------------------------------------------------------

img_sense_topup = [];
    
file_path = fleet_path;

cd(file_path)

num_dwi = 1;
 
for ndwi = 1:num_dwi
    disp(['dwi: ', num2str(ndwi)])
    
    system(['rm ', file_path, 'my_blipup_re_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_blipup_im_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_blipdown_re_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_blipdown_im_', num2str(shot_start),'.nii'])
    
    system(['rm ', file_path, 'my_good_blipup_re_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_good_blipup_im_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_good_blipdown_re_', num2str(shot_start),'.nii'])
    system(['rm ', file_path, 'my_good_blipdown_im_', num2str(shot_start),'.nii'])
    

    genNii( real(img_recon_pad(:,:,:,1)), voxel_size, [file_path, 'my_blipup_re_', num2str(shot_start),'.nii'])
    genNii( imag(img_recon_pad(:,:,:,1)), voxel_size, [file_path, 'my_blipup_im_', num2str(shot_start),'.nii'])

    genNii( real(img_recon_pad(:,:,:,2)), voxel_size, [file_path, 'my_blipdown_re_', num2str(shot_start),'.nii'])
    genNii( imag(img_recon_pad(:,:,:,2)), voxel_size, [file_path, 'my_blipdown_im_', num2str(shot_start),'.nii'])



    cmd = ['applytopup --imain=my_blipup_re_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ' --method=jac --out=my_good_blipup_re_', num2str(shot_start), ' --interp=spline'];
    [status, result] = system(cmd, '-echo');  

    
    cmd = ['applytopup --imain=my_blipup_im_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=1 --topup=', fpOut, ' --method=jac --out=my_good_blipup_im_', num2str(shot_start), ' --interp=spline'];
    [status, result] = system(cmd, '-echo');  



    cmd = ['applytopup --imain=my_blipdown_re_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=jac --out=my_good_blipdown_re_', num2str(shot_start), ' --interp=spline'];
    [status, result] = system(cmd, '-echo');  
    

    cmd = ['applytopup --imain=my_blipdown_im_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=jac --out=my_good_blipdown_im_', num2str(shot_start), ' --interp=spline'];
    [status, result] = system(cmd, '-echo');  


    
%     cmd = ['applytopup --imain=my_blipdown_re_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=lsr --out=my_good_blipdown_re2_', num2str(shot_start), ' --interp=spline'];
%     [status, result] = system(cmd, '-echo');  
%     
% 
%     cmd = ['applytopup --imain=my_blipdown_im_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=2 --topup=', fpOut, ' --method=lsr --out=my_good_blipdown_im2_', num2str(shot_start), ' --interp=spline'];
%     [status, result] = system(cmd, '-echo');  

    
    
    system(['gunzip ', file_path, '*.gz'])


    dt_re = load_untouch_nii([file_path, 'my_good_blipup_re_', num2str(shot_start), '.nii']);
    dt_im = load_untouch_nii([file_path, 'my_good_blipup_im_', num2str(shot_start), '.nii']);

    img_1 = dt_re.img + 1i * dt_im.img;



    dt_re = load_untouch_nii([file_path, 'my_good_blipdown_re_', num2str(shot_start), '.nii']);
    dt_im = load_untouch_nii([file_path, 'my_good_blipdown_im_', num2str(shot_start), '.nii']);

    img_2 = dt_re.img + 1i * dt_im.img;


    img_sense_topup(:,:,:,:,ndwi) = permute(cat(4, img_1, img_2), [1,2,4,3]);
end


mosaic(img_recon_pad(:,:,2:5:end,1),2,3,2,'',[-0,1.25e-3], 90)        
mosaic(angle(img_recon_pad(:,:,2:5:end,1)),2,3,12,'',[-pi,pi], 90)        

mosaic(img_recon_pad(:,:,2:5:end,2),2,3,3,'',[-0,1.25e-3], 90)        
mosaic(angle(img_recon_pad(:,:,2:5:end,2)),2,3,13,'',[-pi,pi], 90)        
          

mosaic(img_sense_topup(:,:,1,2:5:end),2,3,12,'',[-0,1.25e-3], 90)      
mosaic(img_sense_topup(:,:,2,2:5:end),2,3,13,'',[-0,1.25e-3], 90)        
              

norm2(img_sense_topup(:,:,2,:)) / norm2(img_sense_topup(:,:,1,:))
norm2(img_recon_pad(:,:,:,2)) / norm2(img_recon_pad(:,:,:,1))


% pad by zeroes back to 150 mtx size
img_sense_topup = padarray(img_sense_topup, [0,1]);

 
% save([fleet_b1k_path, 'img_sense_topup_shot3_MID437'], 'img_sense_topup')


%--------------------------------------------------------------------------
%% try to get a single combined image from real_blipup and real_blipdown
%% and another combined image from imag_blipup and imag_blipdown
%--------------------------------------------------------------------------

file_path = fleet_path;

genNii( real(img_recon_pad(:,:,:,1)), voxel_size, [file_path, 'my_blipup_re_', num2str(shot_start),'.nii'])
genNii( imag(img_recon_pad(:,:,:,1)), voxel_size, [file_path, 'my_blipup_im_', num2str(shot_start),'.nii'])

genNii( real(img_recon_pad(:,:,:,2)), voxel_size, [file_path, 'my_blipdown_re_', num2str(shot_start),'.nii'])
genNii( imag(img_recon_pad(:,:,:,2)), voxel_size, [file_path, 'my_blipdown_im_', num2str(shot_start),'.nii'])

% applytopup --imain=my_blipup,my_blipdown --datain=my_parameters --inindex=1,2 --topup=my_field --out=my_good_images


fpAcqp = [file_path, 'acq_param.txt'];
fpOut = [file_path, 'out_ap_pa_',num2str(ii_dif),'rf',num2str(ii_rf)];

cmd = ['applytopup --imain=', file_path, 'my_blipup_re_', num2str(shot_start), ',', file_path, 'my_blipdown_re_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=1,2 --topup=', fpOut, ' --out=', file_path, 'my_good_images_re'];
[status, result] = system(cmd, '-echo');  

cmd = ['applytopup --imain=', file_path, 'my_blipup_im_', num2str(shot_start), ',', file_path, 'my_blipdown_im_', num2str(shot_start), ' --datain=', fpAcqp, ' --inindex=1,2 --topup=', fpOut, ' --out=', file_path, 'my_good_images_im'];
[status, result] = system(cmd, '-echo');  


   
system(['gunzip ', file_path, '*.gz'])


dt_re = load_untouch_nii([file_path, 'my_good_images_re.nii']);
dt_im = load_untouch_nii([file_path, 'my_good_images_im.nii']);

img_res = dt_re.img + 1i * dt_im.img;

mosaic(dt_re.img(:,:,2:5:end),2,3,3,'',[-0,1.25e-3], 90)     
mosaic(dt_im.img(:,:,2:5:end),2,3,4,'',[-0,1.25e-3], 90)     



%--------------------------------------------------------------------------
%% hybrid sense with B0 model 
%--------------------------------------------------------------------------

t_axis_ap = [0:Ry:N(2)-1] * esp;

t_axis_pa = t_axis_ap(end:-1:1);

lambda_tik = 1e-6;

rec = zeross([N, num_slice, num_dwi]);
img_hybrid_sense = zeross([N, num_slice, num_dwi]);


tic
for ndwi = 1%:num_dwi
    t1 = img_sense_topup(:,:,1,:,ndwi);
    t2 = img_sense_topup(:,:,2,:,ndwi);
    
    for nslc = 1:num_slice
        disp(['dwi: ', num2str(ndwi), ' / ', num2str(num_dwi), '   slc: ', num2str(nslc), ' / ', num2str(num_slice)])

        phs_diff = exp(-1i* angle(t1(:,:,:,nslc) ./ (eps+t2(:,:,:,nslc))));

        mosaic(angle(phs_diff), 1, 1, 100+ndwi, '', [-pi,pi]/1,90), 

                    
        for xn = 1:N(1)
            b0 = img_fieldmap(xn,:,nslc) * 2 * pi;

            W_ap = exp(1i * t_axis_ap.' * b0);
            AW_ap = A_ap .* W_ap;      

            W_pa = exp(1i * t_axis_pa.' * b0);
            AW_pa = A_pa .* W_pa;      

            AWC_ap = zeross([num_chan*N(2)/Ry, N(2)]);
            AWC_pa = zeross([num_chan*N(2)/Ry, N(2)]);

            for c = 1:num_chan       
                AWC_ap(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :) = AW_ap * diag( sens(xn,:,nslc,c) );

                AWC_pa(1 + (c-1)*N(2)/Ry : c*N(2)/Ry, :) = AW_pa * diag( sens(xn,:,nslc,c) .* phs_diff(xn,:) );
            end

            E = cat(1, AWC_ap, AWC_pa);

            [U,S,V] = svd(E, 'econ');

            E_inv = V * diag(diag(S) ./ (diag(S).^2 + lambda_tik)) * U';


            rhs_ap = sgnl_ap(xn, :, :, nslc, ndwi);   
            rhs_ap = rhs_ap(:);  

            rhs_pa = sgnl_pa(xn, :, :, nslc, ndwi);   
            rhs_pa = rhs_pa(:);  

            rhs = cat(1, rhs_ap, rhs_pa);
    
            rec(xn,:,nslc,ndwi) = E_inv * rhs;
        end
                
    end

    img_hybrid_sense(:,:,:,ndwi) = apodize(rec(:,:,:,ndwi));
end
toc


% save([fleet_b1k_path, 'img_hybrid_sense_shot3_MID437'], 'img_hybrid_sense')

mosaic(img_hybrid_sense(:,:,nslc),1,1,1,'',[0,.25e-3], 90)
mosaic(img_sense_topup(:,:,1,nslc),1,1,2,'',[0,.25e-3], 90)
mosaic(img_sense_topup(:,:,2,nslc),1,1,3,'',[0,.25e-3], 90)

mosaic(img_sense(:,:,1,nslc),1,1,4,'',[0,1e-3], 90)
mosaic(img_sense(:,:,2,nslc),1,1,5,'',[0,1e-3], 90)

mosaic(fft2c(img_hybrid_sense(:,:,nslc)),1,1,11,'',[0,8e-4], 0), setGcf(.6)
 

%--------------------------------------------------------------------------
%%  
%--------------------------------------------------------------------------

imagesc3d2( img_hybrid_sense, s(img_hybrid_sense)/2, 31, [0,0,0]+90, [-0,3e-4])
imagesc3d2( sq(img_sense(:,:,1,:)), s(img_hybrid_sense)/2, 32, [0,0,0]+90, [-0,3e-4])
imagesc3d2( sq(img_sense(:,:,2,:)), s(img_hybrid_sense)/2, 33, [0,0,0]+90, [-0,3e-4])

 

%--------------------------------------------------------------------------
%% fast Buda: gradient descent with B0 model 
%--------------------------------------------------------------------------

addpath(genpath('/autofs/cluster/kawin/Congyu/CODE/msEPI/NEATR_Siemens/'))

load([data_path, 'receive'])

sens = permute(receive, [1,2,4,3]);


% use standard buda fieldmap: improves recon significantly
load('/autofs/space/marduk_001/users/berkin/2020_12_11_bay1_symri_buda_fleet_invivo/standard_buda/img_fieldmap_MID435.mat')


esp = meas.prot.iEffectiveEpiEchoSpacing * 1e-6

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

img_fbuda = zeross([N,2,num_slice,num_dwi]);


tic
for nslc = 1:num_slice
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

        rhs_ap = sgnl_ap(xn, :, :, nslc, ndwi);
        AWC_apHrhs(:,xn) = AWC_apH(:,:,xn) * rhs_ap(:);
        
        rhs_pa = sgnl_pa(xn, :, :, nslc, ndwi);
        
        if use_scaling
            rhs_pa = rhs_pa / scale_fleet;
        end
        
        AWC_paHrhs(:,xn) = AWC_paH(:,:,xn) * rhs_pa(:);
    end    
    

    for ndwi = 1:num_dwi
        % initialize with sense solution: 
        img_sense_ap = sq(img_sense(:,:,1,:,:));
        img_sense_pa = sq(img_sense(:,:,2,:,:));

        im_rec = permute(cat(3, img_sense_ap(:,:,nslc,ndwi), img_sense_pa(:,:,nslc,ndwi)), [2,3,1]);

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

            mosaic(im_rec_disp, 1, 2, 100+vc, ['iter: ', num2str(iter), '  update: ', num2str(rmse(im_prev,im_rec))], genCaxis(im_rec), 90)

            if update < tol
                break
            end
        end

        img_fbuda(:,:,:,nslc,ndwi) = permute(im_rec, [3,1,2]);
    end
end
toc

mosaic(mean(abs(img_fbuda(:,:,1,nslc,1)),3), 1, 1, 16, 'fast buda', [0,1.5e-3], 90)
mosaic(mean(abs(img_fbuda(:,:,2,nslc,1)),3), 1, 1, 17, 'fast buda', [0,1.5e-3], 90)


% save([file_path, 'img_buda_MID447_using_standard_buda_fieldmap'], 'img_fbuda')


%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

mosaic(mean(abs(img_sense(:,:,1,nslc,1)),3), 1, 1, 11, 'sense', [0,1.5e-3], 90)
mosaic(mean(abs(img_sense(:,:,2,nslc,1)),3), 1, 1, 12, 'sense', [0,1.5e-3], 90)


img_sense_topup = padarray(dt2.img, [0,1]);

mosaic(mean(abs(img_sense_topup(:,:,nslc,1)),3), 1, 1, 21, 'sense', [0,1.5e-3], 90)
mosaic(mean(abs(img_sense_topup(:,:,nslc,2)),3), 1, 1, 22, 'sense', [0,1.5e-3], 90)



%--------------------------------------------------------------------------
%% 
%--------------------------------------------------------------------------

% load([file_path, 'img_buda_MID447']);     
load([file_path, 'img_buda_MID447_using_standard_buda_fieldmap']);     
 

Img_real = RealDiffusion_lowRes(permute(img_fbuda, [1,2,4,3]), 0, 0, 0);

% save([file_path, 'Img_real_MID44_using_standard_buda_fieldmap7'], 'Img_real')


imagesc3d2( mean(Img_real,4), s(Img_real)/2, 1, [0,0,0]+90, [-0,13e-4])

imagesc3d2( Img_real(:,:,:,1), s(Img_real)/2, 2, [0,0,0]+90, [-0,13e-4])
imagesc3d2( Img_real(:,:,:,2), s(Img_real)/2, 3, [0,0,0]+90, [-0,13e-4])

%--------------------------------------------------------------------------
%%
%--------------------------------------------------------------------------

mosaic(Img_real(:,:,1:3:end,1), 2, 4, 11, 'sense', [0,1.5e-3], 90)
mosaic(Img_real(:,:,1:3:end,2), 2, 4, 12, 'sense', [0,1.5e-3], 90)

mosaic(mean(Img_real(:,:,1:3:end,:),4), 2, 4, 10, 'sense', [0,1.5e-3], 90)


