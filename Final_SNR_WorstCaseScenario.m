% SNR computation 
% Rewritten by Alessandra, Stella Caporale, 25/10/2024
clearvars
clc
close all

% Directory with subject data
addpath(genpath('/home/stella/Documents/MATLAB/'));  % Add NIfTI toolbox
outpath='/storage/ekaterina/Charmed/Results_R1_SANDI_FR_Final'; % Set outputpath% Set outputpath

initpath = '/media/nas_rete/GLOVE_STUDY/DDC/derivatives';
cd(initpath);
folders = dir('pil*');
N_subj = length(folders);  % Number of subjects
N_diff_vol = 266;

cc_signal_arr = zeros(N_diff_vol,N_subj);
gm_signal_arr = zeros(N_diff_vol,N_subj);
noise_arr = zeros(N_diff_vol,N_subj);
SNR_arr_cc = zeros(N_diff_vol,N_subj);
SNR_arr_gm = zeros(N_diff_vol,N_subj);

for i = 1 : N_subj

    dwifolder = strcat(initpath,'/',folders(i).name,'/resting/dwi');
    cd(dwifolder);
    dwiname = strcat(folders(i).name,'_resting_dwi.nii.gz');
    dwidata = load_untouch_nii(dwiname);
    dwi = double(dwidata.img);
    xdim = size(dwi,1);
    ydim = size(dwi,2);

    maskname = strcat(dwifolder,'/SANDI_input/',folders(i).name,'_resting_data_mask.nii.gz');
    binmask = load_untouch_nii(maskname);
    binmask = double(binmask.img);
    se = strel('sphere',7);
    binmask2 = imerode(binmask,se);
    erobinmask = binmask2;

    bvalname = strcat(dwifolder,'/SANDI_input/','data_dwi.bval');
    bvals = load(bvalname);
    bvals = bvals';

    bvecname = strcat(dwifolder,'/SANDI_input/',folders(i).name,'_resting_data_dwi.bvec');
    bvecs = load(bvecname);
    bvecs = bvecs';

    vol_label = (1:1:length(bvals))';
    grad_table = cat(2,bvals,bvecs,vol_label);

    sorted_grad_table = sortrows(grad_table,1);

    FA = load_untouch_nii(strcat(dwifolder,'/DTI/',folders(i).name,'_resting_dti_FA.nii.gz'));
    FA = double(FA.img);
    V1 = load_untouch_nii(strcat(dwifolder,'/DTI/',folders(i).name,'_resting_dti_V1.nii.gz'));
    V1 = double(V1.img);

    %% masking of corpus callosum
    FA_thr = 0.75;
    FA_mask = FA;
    FA_mask(FA_mask>=FA_thr & FA_mask<1)=1;   % FA_mask<1 this eliminates spourious FA=1.2 due to eddy currents
    FA_mask(FA_mask~=1)=0;

    V1_thr = 0.9;
    V1_mask = squeeze(V1(:,:,:,1));
    V1_mask(V1_mask>=V1_thr)=1;
    V1_mask(V1_mask<=-V1_thr)=1;
    V1_mask(V1_mask~=1)=0;

    cc_mask = FA_mask.*V1_mask.*erobinmask;

    figure(1)
    for ss = 1 : 10
        subplot(2,5,ss);
        imshow(cc_mask(:,:,14+ss),[]);
    end
    figure(2)
    for ss = 1 : 5
        subplot(5,1,ss);
        tmp = squeeze(cc_mask(xdim/2+ss,:,:));
        imshow(rot90(tmp),[]);
    end

    %% masking of GM

    FA_thr1 = 0.15; FA_thr2 = 0.38;
    FA_mask = FA;
    FA_mask(FA_mask>=FA_thr1 & FA_mask<=FA_thr2)=1;
    FA_mask(FA_mask~=1)=0;

    gm_mask = FA_mask.*erobinmask;

    figure(1)
    for ss = 1 : 10
        subplot(2,5,ss);
        imshow(gm_mask(:,:,14+ss),[]);
    end
    figure(2)
    for ss = 1 : 5
        subplot(5,1,ss);
        tmp = squeeze(gm_mask(xdim/2+ss,:,:));
        imshow(rot90(tmp),[]);
    end

    %% bkg mask
    bkg_mask = zeros(size(FA));
    roi_size = xdim/5.5;
    bkg_mask(1:roi_size,1:roi_size,:) = 1;
    bkg_mask(end-(roi_size-1):end,end-(roi_size-1):end,:) = 1;
    bkg_mask(1:roi_size,end-(roi_size-1):end,:) = 1;
    bkg_mask(end-(roi_size-1):end,1:roi_size,:) = 1;

    dwi_bkg = dwi.*bkg_mask;

    figure(3)
    subplot(1,2,1)
    imshow(bkg_mask(:,:,round(size(FA,3)/2)));
    subplot(1,2,2)
    imshow(dwi_bkg(:,:,round(size(FA,3)/2)),[])

    %% SNR in CC

    for bb = 1 : length(bvals)
        vol = sorted_grad_table(bb,5);
        dwivol = squeeze(dwi(:,:,:,vol));
        tmp_signal = cc_mask.*dwivol;
        cc_signal_arr(bb,i) = mean(tmp_signal(tmp_signal~=0));

        tmp_signal = gm_mask.*dwivol;
        gm_signal_arr(bb,i) = mean(tmp_signal(tmp_signal~=0));

        tmp_bkg = bkg_mask.*dwivol;
        noise_arr(bb,i) = std(tmp_bkg(tmp_bkg~=0));
    end

    bval = sorted_grad_table(:,1);
    SNR_arr_cc(:,i) = cc_signal_arr(:,i)./noise_arr(:,i);
    SNR_arr_gm(:,i) = gm_signal_arr(:,i)./noise_arr(:,i);
    
    h = figure;
    plot(bval,SNR_arr_cc(:,i),'ro');
    ylabel('SNR');
    xlabel('b-value (s/mm^2)');
    xlim([-200 6200]);
    h.Color = [1 1 1];
    h.Position = [1 1 1000 700];
    set(gca,'FontSize',12,'FontWeight','Bold');
    hold on
    plot(bval,SNR_arr_gm(:,i),'bo');
    grid on
    legend('CC','GM');
    figpath = '/storage/ekaterina/Charmed/Results_R1_SANDI_FR_Final/';
    saveas(h,strcat(figpath,'SNR_cc_gm_',folders(i).name),'fig');
    % from the plot, I can pinpoint the min SNR for each b-value
    % and check that it corresponds to the L-R direction.
    % with this command I get an index N --> I go and check the Nth
    % line in the sorted_grad_table, and that should correspond to the
    % [1 0 0] direction
    % The max SNR for each b-value should correspond to the [0 1 0] or [0 0 1]
    % direction
    %[m1,n1]=find(abs(SNR_arr_cc-36.209)<0.001); % min for b=500 s/mm2 [1 0 0]
    %[m2,n2]=find(abs(SNR_arr_cc-71.011)<0.001); % max for b=500 s/mm2 [0 1 0]
    close all
end
%% SNR vs b-value - average across subjects 
SNR_across_cc = mean(SNR_arr_cc,2);
SNR_across_gm = mean(SNR_arr_gm,2);

h = figure;
plot(bval,SNR_across_cc,'ro');
ylabel('SNR');
xlabel('b-value (s/mm^2)');
xlim([-200 6200]);
h.Color = [1 1 1];
h.Position = [1 1 1000 700];
set(gca,'FontSize',12,'FontWeight','Bold');
hold on
plot(bval,SNR_across_gm,'bo');
grid on
legend('CC','GM');
figpath = '/storage/ekaterina/Charmed/Results_R1_SANDI_FR_Final/';
saveas(h,strcat(figpath,'SNR_cc_gm_across_subjs'),'fig');

cd(figpath);

%% SNR vs b-value - average across directions
N_bvals = 7;
bval_arr = [0 200 500 1200 2400 4000 6000];
SNR_acrossDir_cc = zeros(N_subj,N_bvals);
SNR_acrossDir_gm = zeros(N_subj,N_bvals);

idx_0 = find(bval==0);
idx_200 = find(bval==200);
idx_500 = find(bval==500);
idx_1200 = find(bval==1195 | bval==1200 | bval==1205);
idx_2400 = find(bval==2395 | bval==2400 | bval==2405);
idx_4000 = find(bval==3995 | bval==4000 | bval==4005);
idx_6000 = find(bval==5990 | bval==5995 | bval==6000 | bval==6005 | bval==6010);

% GM
tmp = SNR_arr_gm(idx_0,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,1) = tmp_mean;

tmp = SNR_arr_gm(idx_200,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,2) = tmp_mean;

tmp = SNR_arr_gm(idx_500,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,3) = tmp_mean;

tmp = SNR_arr_gm(idx_1200,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,4) = tmp_mean;

tmp = SNR_arr_gm(idx_2400,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,5) = tmp_mean;

tmp = SNR_arr_gm(idx_4000,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,6) = tmp_mean;

tmp = SNR_arr_gm(idx_6000,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_gm(:,7) = tmp_mean;

% WM
tmp = SNR_arr_cc(idx_0,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,1) = tmp_mean;

tmp = SNR_arr_cc(idx_200,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,2) = tmp_mean;

tmp = SNR_arr_cc(idx_500,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,3) = tmp_mean;

tmp = SNR_arr_cc(idx_1200,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,4) = tmp_mean;

tmp = SNR_arr_cc(idx_2400,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,5) = tmp_mean;

tmp = SNR_arr_cc(idx_4000,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,6) = tmp_mean;

tmp = SNR_arr_cc(idx_6000,:);
tmp_mean = mean(tmp,1);
SNR_acrossDir_cc(:,7) = tmp_mean;

h = figure;
h.Color = [1 1 1];
h.Position = [1 1 1000 700];

subplot(2,2,1);
plot(bval_arr,SNR_acrossDir_gm,'b*');
grid on
ylabel('SNR');
xlabel('b-value (s/mm^2)');
xlim([-200 6200]);
ylim([0 150]);
legend('GM SNR');
set(gca,'FontSize',12,'FontWeight','Bold');

subplot(2,2,2);
set(gca,'FontSize',12,'FontWeight','Bold');
plot(bval,SNR_across_gm,'bo');
grid on
ylabel('SNR');
xlabel('b-value (s/mm^2)');
xlim([-200 6200]);
ylim([0 150]);
legend('GM SNR');
set(gca,'FontSize',12,'FontWeight','Bold');

subplot(2,2,3);
plot(bval_arr,SNR_acrossDir_cc,'r*');
grid on
ylabel('SNR');
xlabel('b-value (s/mm^2)');
xlim([-200 6200]);
ylim([0 150]);
legend('CC SNR');
set(gca,'FontSize',12,'FontWeight','Bold');

subplot(2,2,4);
plot(bval,SNR_across_cc,'ro');
grid on
ylabel('SNR');
xlabel('b-value (s/mm^2)');
xlim([-200 6200]);
ylim([0 150]);
legend('CC SNR');
set(gca,'FontSize',12,'FontWeight','Bold');

cd(figpath);

figpath = '/storage/ekaterina/Charmed/Results_R1_SANDI_FR_Final/';
saveas(h,strcat(figpath,'SNR_cc_gm_across_Subjs_Dirs'),'fig');

close all