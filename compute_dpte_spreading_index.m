
function [SI,Zscores,Zmap] = compute_dpte_spreading_index(dPTE,dPTE_ctrl,Atlas,FidelityOperator,plimit,source_size,stim_side,response_side)

% A script for calculation of Spreading Index (SI) of corticokinematic coherence (CKC) response network from 
% directed phase transfer entropy (dPTE) interaction matrices of movement stimulation and control condition recordings.
% As used in the article (2022): "Cortical networks show characteristic recruitment patterns after somatosensory stimulation by pneumatically evoked repetitive hand movements in newborn infants"
% Authors: E Ahtola, S Leikos, A Tuiskula, L Haataja, E Smeds, H Piitulainen, V Jousmäki, A Tokariev, S Vanhatalo
% 
% Corticokinematic coherence (CKC) reflects coupling between cortical EEG and peripheral kinematic signals. 
% Spreading Index (SI) quantifies the extent of the information spread from the subset of parcels that respond most prominently to the stimulation.
% The information flow (dPTE) during the movement stimulation is compared to reference data of a control recording. 
% The connections significantly different from the surrogate data at a given alpha level are identifies using a right-tailed Z-test.
% SI ranges between 0 and 100% and is defined as the proportion between the number of significant edges and all possible edges originating from the primary sources.
%
% INPUT:
% dPTE              dPTE interaction matrix of movement stimulation
% dPTE_ctrl         dPTE interaction matrix of control condition (baseline reference)
% Atlas             Atlas information of the parcellation
% FidelityOperator  Matrix for masking unreliable connections (due to suboptimal number and layout of recording electrodes)
% plimit            Alpha-level for the statisctical testing (e.g. 0.01)
% source_size       Selected size for primary source parcels (e.g. 4)
% stim_side         Stimulation side: 'right' / 'left' (hand)
% response_side     Selected hemispere: 'contra' / 'ipsi' (in relation to the stimulation side)
% 
% OUTPUT:
% SI                Spreading Index: Percentage of significant connections from primary source parcels of the selected hemisphere
% Zscores             Z-score matrix (right/left stim. dPTE compared to control condition)
% Zmap              Binary thresholded Z-score matrix. Limit speficied in plimit.
%
% EXAMPLE DATA: dPTE_example_data.mat
% dPTE_right        dPTE interaction matrix of right hand movement stimulation
% dPTE_control      dPTE interaction matrix of control condition
% MyAtlas           Atlas information of the parcellation
% FidelityOperator  Matrix for masking unreliable connections
%
% Test by runnning commands:
% load dPTE_example_data;
% [SI,Zscores,Zmap] = compute_dpte_spreading_index(dPTE_right,dPTE_control,MyAtlas,FidelityOperator,0.01,4,'right','contra');
%
% Eero Ahtola 4.5.2022
% eero.ahtola@hus.fi
%

%Replace FidelityOperator zeros with NaNs
FidOpNaN = FidelityOperator;
FidOpNaN(FidOpNaN==0) = NaN;

%Reference interaction matrix
Ref_dPTE = dPTE_ctrl .* FidOpNaN;

%Z-score for each sample: z = (x – mean(ref)) / std(ref)
Zscores = (dPTE - nanmean(Ref_dPTE(:))) ./ nanstd(Ref_dPTE(:)); 
 
%Compare to plimit
%Statistical Test used: right-tailed z-test
above = Zscores > norminv(1-plimit);
above = above .* FidOpNaN;
Zmap = above;  %Store output matrix

%Identify contra- and ipsilateral parcels from the Atlas
pick_ipsi_prc = [];    
pick_contra_prc = [];    
if strcmpi(stim_side,'right')        
    pick_ipsi_prc = find(strcmpi(Atlas.Hemisphere,'R')); %stim right (ipsi)    
    pick_contra_prc = find(strcmpi(Atlas.Hemisphere,'L')); %stim right (contra)    
elseif strcmpi(stim_side,'left')
    pick_ipsi_prc = find(strcmpi(Atlas.Hemisphere,'L')); %stim left (ipsi)    
    pick_contra_prc = find(strcmpi(Atlas.Hemisphere,'R')); %stim left (contra)
else
    disp('Set correct input for stim_side'); keyboard;
end

%Calculate spreading index (SI)
%SI: Mean of all outbound connections from primary source parcels

% Find contralateral parcels that provide most significant edges
P_signif_edges=nanmean(above,2);  %avg

%Set ipsi-/contra-lateral side to zero
if strcmpi(response_side,'contra')    
    P_signif_edges(pick_ipsi_prc)=0;      %Set ipsilateral side to zero        
elseif strcmpi(response_side,'ipsi')
    P_signif_edges(pick_contra_prc)=0;    %Set contralateral side to zero    
else
    disp('Set correct input for response_side'); keyboard;
end      

%Identify primary source parcels
[~,prc_ind_sorted] = sort(P_signif_edges);                       %sort values
select_prc1 = sort(prc_ind_sorted(end-source_size+1:end));       %select highest (N=source_size)

%Calculate SI from the proportion between the number of edges (significant vs all possible) 
A1=above(select_prc1,:);    
a1 = nansum(A1(:));                                      %all significant outbound from select_prc1 nodes   
a2 = sum(sum(FidOpNaN(select_prc1,:)==1,2));             %all possible valid connections (from FidelityOperator)
a3 = sum(sum(FidOpNaN(select_prc1,select_prc1)==1,2))/2; %because if A(i,j)>1 -> A(j,i)=0
SI = a1/(a2-a3);   

%Scale to 100%
SI = 100*SI;

end
















