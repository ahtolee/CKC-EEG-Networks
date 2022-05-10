# CKC-EEG-Networks
A compilation of MATLAB scripts for performing network-based corticokinematic coherence (CKC) analysis for EEG data recorded under repetitive hand movement stimulation.
The code was developed for the study (2022): "Cortical networks show characteristic recruitment patterns after somatosensory stimulation by pneumatically evoked repetitive hand movements in newborn infants", Authors: E. Ahtola, S. Leikos, A. Tuiskula, L. Haataja, E. Smeds, H. Piitulainen, V. Jousmäki, A. Tokariev, S. Vanhatalo.

Function compute_ckc_responses.m is used for detection of CKC responses from CSD-transformed EEG-signals. CKC reflects coupling between cortical EEG and peripheral kinematic signals. The algorithm is based on inter-trial phase coherence (ITC) that yields a measure of phase-locked synchronization between repeated EEG trials in relation to the stimulus events.

Function compute_ckc_connectivity.m is used to characterize the cortical network related to a CKC response from pre-filtered parcel signals (EEG). Connectivity between signal pairs is calculated using directed phase transfer entropy that estimates the preferred direction of the information flow.

Function compute_dpte_spreading_index.m is used for calculation of Spreading Index (SI) of CKC response network from directed phase transfer entropy (dPTE) interaction matrices of movement stimulation and control condition recordings. SI quantifies the extent of the information spread from the subset of parcels that respond most prominently to the stimulation. The information flow during the movement stimulation is compared to reference data of a control recording.

Attached mat-files contain authentic test data from one CKC recording (right hand stimulation at 1.79 Hz) in a form fully compatible with the script. There is a mat-file for each analysis function separately. Instructions for the use are provided in the header part of the m-files.

External m-files required by the analysis are attached as a separate zip package.

Eero Ahtola. 2022. eero.ahtola (at) hus.fi
