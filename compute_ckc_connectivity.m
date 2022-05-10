function [dPTE, dPTE_FidOp] = compute_ckc_connectivity(EEG, Atlas, FidelityOperator, cx_surface)

% A script for calculation of corticokinematic coherence (CKC) response network from pre-filtered parcel signals (EEG).
% As used in the article (2022): "Cortical networks show characteristic recruitment patterns after somatosensory stimulation by pneumatically evoked repetitive hand movements in newborn infants"
% Authors: E Ahtola, S Leikos, A Tuiskula, L Haataja, E Smeds, H Piitulainen, V Jousmäki, A Tokariev, S Vanhatalo
% 
% Corticokinematic coherence (CKC) reflects coupling between cortical EEG and peripheral kinematic signals.
% Connectivity between signal pairs is calculated using phase transfer entropy. We use algorithm presented by Hillebrand et al. (2016),
% where values between signals x and y (PTExy from x to y and PTEyx from y to x) are normalised together to yield a measure
% of directed phase transfer entropy (dPTE) that estimates the preferred direction of the information flow. The normalised dPTE values
% range between 0 and 1, with value of 0.5 indicating equal causality between the signals.
%
% Function produces three figures:
% 1) dPTE values presented as an interaction matrix (cf. Fig S7).
% 2) Circular graphs of CKC network showing the connections with highest dPTE (cf. Fig S1).
% 3) CKC network graph with the strongest connections projected over a 3D cortex surface model.
%
% INPUT:
% EEG               Parcel signals as a matrix (channels x data x trials). Data must be pre-filtered around the response frequency, e.g. f1 +/- 0.3Hz.
% Atlas             Atlas information of the parcellation. For visualization of the network connections.
% FidelityOperator  Matrix for masking unreliable connections (due to suboptimal number and layout of recording electrodes)  
% cx_surface        3D surface model of cortex. For visualization of the network connections.
% 
% OUTPUT:
% dPTE              dPTE matrix showing the preferred direction of the information flow measured with phase transfer entropy  
% dPTE_FidOp        dPTE matrix masked with the FidelityOperator (unreliable connections as NaNs).
%
% EXAMPLE DATA: parcel_signals_trial_data.mat
% EEG               Parcel signal trials from a single CKC recording (right hand stimulation; Fs=250 Hz). Data is filtered around the response frequency 3.571Hz +/-0.3Hz.
% MyAtlas           Atlas information of the parcel transformation
% FidelityOperator  Matrix for masking unreliable connections
% real_cx           3D surface model of cortex.
%
% Test by runnning commands:
% load parcel_signals_trial_data.mat;
% [dPTE, dPTE_FidOp] = compute_ckc_connectivity(EEG, MyAtlas, FidelityOperator, real_cx)
%
% Script requires following external functions (all included):
% PhaseTE_MF.m: Matteo Fraschini and Arjan Hillebrand 2016. Function is included in Brainstorm Toolbox (Tadel et al. 2011; https://neuroimage.usc.edu/brainstorm/).
% arrow.m: Erik A Johnson 2016. Retrieved from MATLAB Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/278-arrow).
% arrow3.m: Tom Davis 2002-2020. Retrieved from MATLAB Central File Exchange (https://www.mathworks.com/matlabcentral/fileexchange/14056-arrow3).
%
% Eero Ahtola 4.5.2022
% eero.ahtola@hus.fi
%

%Concatenate trials together
[Nch,pnts,Ntrials] = size(EEG);
pnts = int64(pnts);
parcel_signals=reshape(EEG, Nch, Ntrials*double(pnts))';

%Calculate Phase Transfer Entropy (PTE) with external function PhaseTE_MF.m.
[dPTE, PTE] = PhaseTE_MF(parcel_signals, [], 'scott');  
disp('-Phase Transfer Entropy (dPTE) calculated-');

%Create new ChNames vector (labels) based on Atlas information
Nprc = length(Atlas.Parcels);
ChNames = cellstr(num2str([1:Nprc]'));      %Generate new Channel names based on the parcel numbers
ChNames = strrep(ChNames, ' ', '');         %remove empty spaces    


% Visualizations of the network connections

% Fig1. Plot dPTE as an interaction matrix 
A = dPTE .* FidelityOperator; A(A==0)=0.5;  %Remove unreliable connections with FidelityOperator mask

figure; set(gcf,'Color','w');
imagesc(A);
colormap('jet');
caxis([0.43 0.57]);
set(gca,'XTick',[10:10:58]); set(gca,'YTick',[10:10:58]);
xlabel('Sink parcel'); ylabel('Source parcel');
set(gca,'fontsize',12.5);
c=colorbar;
ax = gca; axpos = ax.Position; cpos = c.Position;
cpos(1)=cpos(1)-0.1*cpos(3); cpos(3) = 0.3*cpos(3);
c.Position = cpos; ax.Position = axpos;
set(c,'fontsize',11);
set(gcf,'position',[122 228 468 416]);
daspect([1 1 1]);

% Fig2. Plot circular dPTE connection graphs (include connections with dPTE > 0.54)
% Thresholds for connection colour scale: >0.54 >0.57 >0.60 >0.63
plot_dpte_circular(dPTE,Atlas,FidelityOperator);
set(gcf,'position',[16, 6, 11 11]);


% Fig3. Plot network of the strongest connections over 3D brain model (dPTE>0.57)
A = dPTE .* FidelityOperator; %Remove unreliable connections with FidelityOperator mask
A = A>=0.57;                  %Include only connections with dPTE>0.57

%Specify view angle for the 3D image
az1=270; el1=90;  %View angle from top
%az1=270; el1=0;  %View angle from behind
%az1=180; el1=0;  %View angle from left
%az1=0; el1=0;    %View angle from right

plot_3D_network(A,Atlas, cx_surface);    
view(az1,el1); set(gcf,'position',[1035  227 480 420]); 

%Apply FidelityOperator for the output dPTE matrix:
FidOpNaN = FidelityOperator;
FidOpNaN(FidOpNaN==0) = NaN;   %Convert zeros to NaNs
dPTE_FidOp = dPTE .* FidOpNaN; %Mark unreliable connections with NaNs (use FidelityOperator mask)

end


% Internal functions

function plot_dpte_circular(A,Atlas,FidelityOperator)

%Create a circular graph for dPTE results

%Inputs 
% A                 dPTE matrix
% Atlas             Atlas information of the parcellation. For visualization of the network connections.
% FidelityOperator  Matrix for masking unreliable connections (due to suboptimal number and layout of recording electrodes)  

predef_pte_lims = [0.54 0.57 0.60 0.63];    %Thresholds for connection colour scale  

% Remove unreliable connections with FidelityOperator mask
A = A.*FidelityOperator;
A(A==0)=NaN;

%Select strongest p% of the connections
pte_limits = predef_pte_lims;  % use predefined values: limits: limit to colours
%pte_scale = predef_pte_lims;   % use predefined values: scale: which edges are drawn 

% N of parcels  
Np = length(Atlas.Parcels);

% Define connection color coding:
clr = [];
clr(1, :) = [0.6816 0.8714 0.9706];  %1st color
clr(2, :) = [0.2588 0.5725 0.7765];  %2nd color
clr(3, :) = [0.0314 0.1882 0.4196];  %3rd color    
clr(4, :) = [0.0126 0.1270 0.2447];  %4th color    
[Nclr,~]=size(clr);  

% Create a new figure
figure; set(gcf, 'Color', 'w'); hold on    
axis off; xlim([-0.5 0.5]); ylim([-0.5 0.5]);   

t = linspace(0, 1, 501); %vector for connection trace drawing
Nmap=200;                %for color scaling
all_clr_inds=[0:Nmap,Inf];        

for  c=1:length(all_clr_inds)-1   %Drawing order from low to high dPTE
    
  %Loop over channel pairs  
  for ChA = 1:Np
     for ChB = 1:Np                             
          clr_ind=round((A(ChA, ChB)-pte_limits(1))/(pte_limits(end)-pte_limits(1))*Nmap);  %Scale between pte_limits values (Nmap steps)
          if clr_ind>=all_clr_inds(c) && clr_ind<all_clr_inds(c+1)
              
             %Arrow parameters -> uses arrow.m
             ArrowLength = 12;
             d1=9; d2=6;  %arrow start/stop displacement values
                                         
             %1st color                 
             if A(ChA, ChB) > 0 && A(ChA, ChB) <= pte_limits(1)                  
                pts = kron((1-t).^2, Atlas.Circular_xy(ChA, :)') + kron(2*(1-t).^t, [0; 0]) + kron(t.^2, Atlas.Circular_xy(ChB, :)'); 
                plot(pts(1, :), pts(2, :), 'Color', clr(1, :), 'Linewidth', 0.75);
                Start = pts(:,end-d1)'; Stop = pts(:,end-d2)';                                                 
                arrow(Start,Stop,ArrowLength,'EdgeColor',clr(1, :),'FaceColor',clr(1, :)); %arrow function from file central
             end                                       
             %2nd color
             for nclr=2:Nclr
                 if A(ChA, ChB) >= pte_limits(nclr-1) && A(ChA, ChB) <= pte_limits(nclr)                                   
                    pts = kron((1-t).^2, Atlas.Circular_xy(ChA, :)') + kron(2*(1-t).^t, [0; 0]) + kron(t.^2, Atlas.Circular_xy(ChB, :)'); 
                    plot(pts(1, :), pts(2, :), 'Color', clr(nclr-1, :), 'Linewidth', 0.75);                                     
                    Start = pts(:,end-d1)'; Stop = pts(:,end-d2)';                                               
                    arrow(Start,Stop,ArrowLength,'EdgeColor',clr(nclr-1, :),'FaceColor',clr(nclr, :)); %arrow function from file central
                 end 
             end                 
             %3rd color
             if A(ChA, ChB) > pte_limits(end)                 
                pts = kron((1-t).^2, Atlas.Circular_xy(ChA, :)') + kron(2*(1-t).^t, [0; 0]) + kron(t.^2, Atlas.Circular_xy(ChB, :)');  
                plot(pts(1, :), pts(2, :), 'Color', clr(Nclr, :), 'Linewidth', 0.75);
                Start = pts(:,end-d1)'; Stop = pts(:,end-d2)';                                                
                arrow(Start,Stop,ArrowLength,'EdgeColor',clr(Nclr, :),'FaceColor',clr(Nclr, :)); %arrow function from file central                
             end             
          end
      end
  end
end  

%Set fig size
set(gcf, 'units', 'centimeters', 'pos', [36 15 10 10]);   
set(gca, 'Position', [0.07 0.07 0.86 0.86]);

% Draw nodes  
scatter(Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'O', 1), Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'O', 2), 30, 'ok', 'filled');
scatter(Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'C', 1), Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'C', 2), 30, 'MarkerFaceColor', [0.75 0.00 0.75], 'MarkerEdgeColor', [0.75 0.00 0.75]);
scatter(Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'F', 1), Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'F', 2), 30, 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', [0.93 0.69 0.13]);
scatter(Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'T', 1), Atlas.Circular_xy(cell2mat(Atlas.Areas) == 'T', 2), 30, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerEdgeColor', [0.47 0.67 0.19]);   

% Draw parcel numbers
for j = 1:Np
  text(1.06*Atlas.Circular_xy(j, 1), 1.06*Atlas.Circular_xy(j, 2), num2str(j), 'FontSize', 8, 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end  

% Draw Regional tags
Regions = {'O','O','T','T','C','C','F','F'};
%Positions:
Reg_xy = [-0.25 -0.52; 0.25 -0.52;...
      -0.48 -0.31; 0.48 -0.31;...
      -0.53  0.17; 0.53  0.17;...
      -0.29  0.49; 0.29  0.49];
%Colors:
Reg_clr = [0 0 0; 0 0 0; 0.47 0.67 0.19; 0.47 0.67 0.19; 0.75 0.00 0.75; 0.75 0.00 0.75; 0.93 0.69 0.13; 0.93 0.69 0.13];
for j = 1:length(Regions)        
    text(Reg_xy(j,1), Reg_xy(j,2), Regions{j}, 'FontSize', 12, 'Color', Reg_clr(j,:), 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');  
end 
    hold off
    
end


function plot_3D_network(A, Atlas, cx_surface)

%Draw 3D cortex surface model and superimpose network connections

%inputs 
% A                 Binary network array (e.g. thresholded dPTE matrix)
% Atlas             Atlas information of the parcellation.
% cx_surface        3D surface model of cortex.

figw=750;  %default fig width

% N of parcels  
Np = length(Atlas.Parcels);     
    
% Create a figure
figure; hold on
set(gcf,'position',[900 500 figw figw/1.33]); 
set(gcf, 'Color', 'w');

a_real = 0.1;  %alpha value for cx surface (transparency)

% Show brain (real)
patch('Vertices', cx_surface.Vertices, 'Faces', cx_surface.Faces, 'FaceColor', [0.6 0.6 0.6], 'FaceAlpha', a_real, 'EdgeColor', 'none');  

%Mark nodes (parcel centroids) of the network with dots
marked_nodes = cell2mat(Atlas.Areas) ~= 'X';  %All nodes (except Xs)        
    
% Mark centroids with colors according to location   
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'O' & marked_nodes>0, 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O' & marked_nodes>0, 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'O' & marked_nodes>0, 3), 40, 'ok', 'filled');  
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'T' & marked_nodes>0, 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T' & marked_nodes>0, 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'T' & marked_nodes>0, 3), 40, 'MarkerFaceColor', [0.47 0.67 0.19], 'MarkerEdgeColor', [0.47 0.67 0.19]);
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'C' & marked_nodes>0, 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C' & marked_nodes>0, 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'C' & marked_nodes>0, 3), 40, 'MarkerFaceColor', [0.75 0.00 0.75], 'MarkerEdgeColor', [0.75 0.00 0.75]);
scatter3(Atlas.Centroids(cell2mat(Atlas.Areas) == 'F' & marked_nodes>0, 1), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F' & marked_nodes>0, 2), Atlas.Centroids(cell2mat(Atlas.Areas) == 'F' & marked_nodes>0, 3), 40, 'MarkerFaceColor', [0.93 0.69 0.13], 'MarkerEdgeColor', [0.93 0.69 0.13]);
    
%Mark connections of the network with lines
% Arrow heads for asymmetric connections -> arrow3.m

%Line and arrow parameters:
ln_clr=[0 0.4470 0.7410]; %Blue color line
LW=1;                     %Line width
arrow_type = '-_b0';      %Arrow head
arwidth=0.75;             %Arrow width

%Loop over channel pairs
for ChA=1:Np
  for ChB=1:Np          
      if A(ChA,ChB)>0 %Check is the edge is part of the network  

          %Draw line
          ln1=line(Atlas.Centroids([ChA;ChB],1),Atlas.Centroids([ChA;ChB],2),Atlas.Centroids([ChA;ChB],3),'LineWidth',LW,'Color',ln_clr);                          

          %Draw arrow
          p0 = Atlas.Centroids(ChA,:);
          p2 = Atlas.Centroids(ChB,:);                            
          p1 = p0 + 0.95*(p2-p0);
          p11 = p1 + 0.05*(p2-p0);                  
          arrow3(p1,p11,arrow_type,arwidth,1.5);                                     
      end      
  end
end

% Show parcel numbers
for j = 1:Np    
  %black:
  text(1.1*Atlas.Centroids(j, 1), 1.1*Atlas.Centroids(j, 2), 1.1*Atlas.Centroids(j, 3), num2str(j), 'FontSize', 7, 'Color', 'k', 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
end  

axis off;
hold off ;  
  
end










