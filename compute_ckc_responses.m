function [ckc,pval] = compute_ckc_responses(EEG, Fs, target_freq, ChNames, eloc)

% A script for detection of corticokinematic coherence (CKC) responses from CSD signals (EEG).
% As used in the article (2022): "Cortical networks show characteristic recruitment patterns after somatosensory stimulation by pneumatically evoked repetitive hand movements in newborn infants"
% Authors: E Ahtola, S Leikos, A Tuiskula, L Haataja, E Smeds, H Piitulainen, V Jousmäki, A Tokariev, S Vanhatalo
% 
% Corticokinematic coherence (CKC) reflects coupling between cortical EEG and peripheral kinematic signals.
% This script can be used to detect CKC using inter-trial phase coherence (ITC) that yields a measure of phase-locked
% synchronisation between repeated EEG trials in relation to the stimulus events (e.g. Tallon-Baudry et al. 1996; Delorme and Makeig 2004).
%
% Here CKC is defined as an average of ITC t–f decomposition values from the response frequency band (3.56 Hz) excluding the first and last 200 ms of the epoch length.
% CKC/ITC values are tranformed into p-values using equation p = exp(sqrt(1+4*N+4*(N^2-(N*itc).^2))-(1+2*N)), where N is the trial count.
%
% Function produces two graphs:
% 1) T-f figure plate that shows ITC values of each channel as a time
% frequency plot. CKC is an average of ITC values at the response frequency (marked as a grey rectangle).  
% 2) Topoplot figure that visualizes CKC values projected on a 2D headmap.
% CKC values are masked by a threshold of p<0.01 (FDR corrected), i.e.
% highlights the statistically significant response components.
%
% INPUT:
% EEG           EEG data as a matrix (channels x data x trials)
% Fs            Sampling frequency (Hz)
% target_freq   Response frequency, we used first harmonic of the stimulation frequency (3.56 Hz)
% ChNames       EEG channel labels as a cell structure
% eloc          Data structure for EEG electrode locations as in EEGLAB (see example data). used in visualization (topoplot)
% 
% OUTPUT:
% ckc           Channel-specific CKC magnitude
% pval          p-values corresponding to the CKC response
%
% EXAMPLE DATA: ckc_csd_trial_data.mat
% EEG           EEG trials from a single CKC recording (right hand stimulation). Data is prefiltered (0.5-30Hz) and in current source density (CSD) format (Perrin et al. 1989).
% Fs            250 Hz
% target_freq   3.5714 Hz
% ChNames       Channel labels (19)
% eloc          Electrode positions
%
% Test by runnning commands:
% load ckc_csd_trial_data.mat;
% [ckc,pval] = compute_ckc_responses(EEG, Fs, target_freq, ChNames, eloc);
%
% Script requires following external functions (included):
% fdr_bh.m. (David Groppe 2022). Retrieved from MATLAB Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/27418-fdr_bh).
%
% Optional script for 2D topographic plot of CKC detection results (Not Included)
% topoplot.m. (Delorme and Makeig 2004) requires freely available EEGLAB toolbox (https://sccn.ucsd.edu/eeglab/index.php).
%
% Eero Ahtola 3.5.2022
% eero.ahtola@hus.fi
%

p_limit = 0.01;         %Alpha value, e.g. threshold for statistically significant responses
fdr = 1;                %Use false discovery rate (FDR) correction for multiple comparisons (1=yes, 0=no)? Uses the procedure by Benjamini and Hochberg (1995).
draw_tflim_edges = 1;   %Draw edges of flim and tlim values to t-f graphs?  1=yes, 0=no
create_topoplot = 1;    %Create 2D topographic plot of the results (1=yes, 0=no)? Requires installation of EEGLAB toolbox.

%Set trial start and end times (trigger at 0 s):
trial_start = -0.2;       % in s 
trial_end = 0.9;          % in s

%Set time and frequency limits for average CKC value ouput
tlim = [0 trial_end-0.2];            %exclude trial edges
flim = target_freq+[-0.2 0.2];       %peak f @ ~3.56 Hz  (Compute CKC from single frequency bin)

[Nch,~,~] = size(EEG);
[~,~,Ntrials] = size(EEG);
trial_t = trial_start:1/Fs:trial_end-(1/Fs);
T_trial = length(trial_t)/Fs;
pnts = int64(length(trial_t));
disp(['Epochs: N=' num2str(Ntrials) ', T=' num2str(T_trial), ', Ttotal=' num2str(Ntrials*T_trial)]);

%Detrend trials
EEG_mod=zeros(size(EEG));
for nch=1:Nch
 for ntrial=1:Ntrials
     epoch1 = EEG(nch,:,ntrial);             
     ep_mod = detrend(epoch1, 'constant');
     EEG_mod(nch,:,ntrial) = ep_mod;
 end
end    
EEG = EEG_mod;
clear EEG_mod;


%T-f tranformation parameters for frequency range
num_freqs = 40;  
min_freq = 0.5;
max_freq = 12;

%Adjust freq bins optimally to match the target frequency
[freqs, min_freq] = match_frex(min_freq, max_freq, num_freqs, target_freq); 

% Wavelet parameters
wave_time = -6:1/Fs:6;                       %Create time vector for wavelets (-x s -> +x s)
half_wave_size = (length(wave_time)-1)/2;
cycle_range = [4 10];
wavecycles = logspace(log10(cycle_range(1)),log10(cycle_range(end)),num_freqs);

% Format FFT variables
nWave = length(wave_time);
nWave = int64(nWave);
nData = pnts*Ntrials;
nData = int64(nData);
nConv = nWave+nData-1;               

%Loop over all EEG channels
ckc = zeros(Nch,1);
pval = zeros(Nch,1);
for nch=1:Nch 
        
    disp(['Processing ch ' num2str(nch) '/' num2str(Nch) ' : ' ChNames{nch}]);
        
    % FFT of data (doesn't change on freq iterations)
    dataX = fft(reshape(EEG(nch,:,:),1,nData),double(nConv));

    % initialize output time-frequency data
    tf = zeros(num_freqs,pnts);        
    tf_rp = zeros(num_freqs,pnts); %Rayleigh test pvalue

    % loop over frequencies
    for fi = 1:num_freqs

        % create Morlet wavelet
        s = wavecycles(fi)/(2*pi*freqs(fi));
        wavelet  = exp(2*1i*pi*freqs(fi).*wave_time) .* exp(-wave_time.^2./(2*s^2));

        %get wavelet FFT
        waveletX = fft(wavelet,double(nConv));  

        % run convolution
        as = ifft(waveletX.*dataX,double(nConv));
        as = as(half_wave_size+1:end-half_wave_size);
        as = reshape(as,pnts,Ntrials);

        % compute ITC (abs mean of phase over trials)
        tf(fi,:) = abs(mean(exp(1i*angle(as)),2));

        %Convert itc values to Rayeigh z-scores and p-values
        [~, pval_pnt] = itc2rayleigh(tf(fi,:), Ntrials);                           
        tf_rp(fi,:) = pval_pnt;
    end
        
    pick_avg_t = [trial_t>=tlim(1) & trial_t<tlim(2)];
    pick_avg_f = [freqs>=flim(1) & freqs<flim(2)];
    ckc(nch) = mean(mean(tf(pick_avg_f,pick_avg_t),2));
    pval(nch) = mean(mean(tf_rp(pick_avg_f,pick_avg_t),2));
             
    %Plot ITC t-f graphs as a subplot for figure plate
    if nch==1  %Create fig plate in first run
        figure; clf; hold on; set(gcf,'Color','w'); box off;    
        set(gcf,'position',[20 20 1400 850]);                    
        %Use 10-20 montage layout:
        scalp_mtg = {'','Fp1','','Fp2','';...
                   'F7','F3','Fz','F4','F8';...
                   'T3','C3','Cz','C4','T4';...
                   'T5','P3','Pz','P4','T6';...
                      '','O1','Oz','O2',''};
        [fpositions,fsize] = subplot_tfpos_for_scalp(ChNames,scalp_mtg); 
    end    
    subplot(fsize(1),fsize(2),fpositions(nch)); hold on;
    contourf(trial_t,freqs,tf,40,'linecolor','none')
    set(gca,'clim',[0 .6],'ydir','normal','xlim',[trial_start trial_end])
    titletext = [ChNames{nch} ', CKC=' num2str(ckc(nch),2)];        
    title(titletext,'Interpreter','none');
    xlabel('Time (s)'); 
    ylabel('Frequency (Hz)'); 
    colormap(gcf,'jet'); 
    caxis([0 0.6]);
    cb=colorbar;
    set(cb,'FontSize',9);  
    if draw_tflim_edges>0  %draw tlim and flim edges       
        rlims = rectangle('Position',[tlim(1),flim(1),diff(tlim),diff(flim)],'EdgeColor',0.6*[1 1 1],'LineWidth',1);        
    end      
    set(gca, ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'YMinorTick'  , 'off'      , ...
        'LineWidth'   , 1         );
    set(gca,'FontName','Helvetica','FontSize',9);
    set(gca,'layer','top');                    
end

% FDR correction for multiple comparisons
if fdr>0 %Use FDR correction for p-values:
    [h, crit_p, ~, ~]=fdr_bh(pval,p_limit,'pdep','no');  %Apply FDR correction. fdr_bh.m available on the Mathworks FileExchange    
    sign_mask = double(h);     
else  %FDR correction not used:
    sign_mask = double(pval <= p_limit & pval > 0); %Create binary mask               
end

% Optional: Visualize CKC results as a topoplot
% Requires topoplot.m from EEGLAB toolbox and locs information for the electrode setup (eloc).
if create_topoplot>0
    figure; hold on; 
    topoplot(ckc,eloc,'electrodes','on','pmask',sign_mask);  %no labels. ks. topoplot.m line 1177
    set(gca, 'CLim', [0, 0.6]);
    set(gcf,'Color','w');
    %Add colorbar
    c=colorbar;
    cbticks=[0: 0.2: 0.6];      
    set(c, 'YTick', cbticks); 
    set(c,'fontsize',12);
end

end


% Internal functions

function [frex_new, new_min_freq] = match_frex(min_freq, max_freq, num_frex, target)

%Help function that creates a logspaced vector from min and max freqencies (num_frex elements).
%Min freq is adjusted so that target frequency is matched exactly at certain element.

%Eero Ahtola 2019

%Try default:
frex = logspace(log10(min_freq),log10(max_freq),num_frex);

%Adjust nearest to target
[~, idx] = min(abs(frex-target));
frex1 = logspace(log10(target),log10(max_freq),num_frex-idx+1); %logspaced from target to max 
d1 = diff(log10(frex1(1:2)));                                   %Linear difference

new_min_freq = 10^(log10(target)-d1*(idx-1));                       %Min freq adjusted so that target is matched at idx
frex_new = logspace(log10(new_min_freq),log10(max_freq),num_frex);  %New frex with new min freq
end


function [fpos,fsize] = subplot_tfpos_for_scalp(ChNames,scalp_mtg)

% Help function that yields location definitions for t-f subplots
% based on given channel labels and requested scalp montage visualization.
%
% E Ahtola 2021

%   1×19 cell array
%   Columns 1 through 10
%     {'Fp1'}    {'Fp2'}    {'F3'}    {'F4'}    {'C3'}    {'C4'}    {'P3'}    {'P4'}    {'O1'}    {'O2'}
%   Columns 11 through 19
%     {'F7'}    {'F8'}    {'T3'}    {'T4'}    {'T5'}    {'T6'}    {'Fz'}    {'Cz'}    {'Pz'}

%             scalp_mtg = {'','Fp1','','Fp2','';...
%                        'F7','F3','Fz','F4','F8';...
%                        'T3','C3','Cz','C4','T4';...
%                        'T5','P3','Pz','P4','T6';...
%                           '','O1','','O2',''};

fsize=size(scalp_mtg);
nums=reshape([1:fsize(1)*fsize(2)],fsize(1),fsize(2))';

fpos=[];
for n=1:length(ChNames)       
    pos1 = nums(find(ismember(scalp_mtg,ChNames{n}))); 
    
    if ~isempty(pos1)
        fpos(n) = pos1;
    else
        disp('Electrode definition (fpos) not found!');
        ChNames{n}
    end      
end            
end


function [rz, rp] = itc2rayleigh(itc,N)

%Convert itc values to Rayeigh z scores and p-values
% N = number of trials
% Rayleigh's Z-score: rz = N*(itc.^2), valid approximation if N is high
% Rayleigh's p-value (coarse):   rp = exp(-N*(itc.^2), valid approximation if N is high
% Rayleigh's p-value (accurate): rp = exp(sqrt(1+4*N+4*(N^2-(N*itc).^2))-(1+2*N))

rz = N*(itc.^2);
rp = exp(sqrt(1+4*N+4*(N^2-(N*itc).^2))-(1+2*N));
end


























































