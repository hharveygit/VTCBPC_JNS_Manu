DEPENDENCIES

- SPM12: https://www.fil.ion.ucl.ac.uk/spm/software/spm12/
- mnl_ieegBasics: https://github.com/MultimodalNeuroimagingLab/mnl_ieegBasics
- matmef: https://github.com/MultimodalNeuroimagingLab/matmef
- mnl_seegview: https://github.com/MultimodalNeuroimagingLab/mnl_seegview
- vistasoft: https://github.com/vistalab/vistasoft
- Freesurfer v7: https://surfer.nmr.mgh.harvard.edu/fswiki/DownloadAndInstall
- GIfTI: https://github.com/gllmflndn/gifti

*****

USAGE

i. The data should be copied to a folder named "data", in the root-level code directory 

ii. Pial, cortical, and subcortical segmentations for each subject were obtained using Freesurfer v7. The relevant Freesurfer outputs for each subject are located in data/derivatives/freesurfer


Step-by-step script sections to generate all outputs are in "main.m". Please open "main.m" and follow the instructions below.
--

1. Calculate MNI 152 coordinates for electrodes in subjects 1 through 4
	ia. Use SPM12 to generate MNI152 forward deformation field from each subject's T1-weighted nifti
		Press Batch. > Spatial -> Segment.
		Specify the subject's T1w nifti file for the "Volumes" input. Select "Forward" for "Deformation Fields". Press the play button.
		Keep the file titled "y_[t1FileName].nii". all other output files can be deleted
	This MNI 152 forward transformation was performed using the non-defaced AC-PC space T1w MRI as input for subject 1 and using the non-defaced native space T1w MRIs as inputs for subjects 2 through 4. Subject 5 was not used in consensus BPC analysis, and therefore not transformed to MNI 152 space.
	ib. Run saveMNIelectrodes.m. Outputs are saved to 'output\sub-X\MNI152\.'

2. Values for Table 1. Run determineStimSites.m to calculate the number of stimulation sites per subject (including stimulated measurement electrodes)

3. Figure 1C-E: Projection matrix, significance matrix, and NNMF illustrations. Run methodfigMatrices.m

4. Figures 2 and 6: Save subject-level BPCs and per-subject wavelet spectrogram, broadband estimates for all measurement electrodes. This will take a while to run
- Run saveBPCs.m

5. Save cortical labels to stim sites in *_BPCSNR.tsv outputs for all electrodes/subjects
- Run assignStimDestrieux.m
- LA1-LA2 stim sites in sub-1 measurement electrodes were manually corrected (in the TSV) afterward to "Left_Amygdala" by visual inspection of CT/T1w MRI overlay

6. Figure 3 and Figure 3-1. Calculate consensus BPCs and plot anatomical distributions on MNI 152 brain. Also calculates cross-validation accuracy of consensus BPC model
- Run clusterBPCs.m
- For cortical and subcortical anatomical legends in Figure 3D (right), run bpcRegions.m

7. Figure 4. Plot stim sites by consensus BPCs on slices of MNI152 T1 MRI
- Run plotBPCs2MNISlices.m, manually adjusting z-coordinate by values on line 38

8. Figure 5. Calculate spectrograms and broadband estimates for consensus BPCs, across collateral sulcus measurement electrodes in subjects 1-4
- Run fig5_exStimSite.m to get the stim site in figure 5A
- Run bpcSpectraBB.m
- Panel A corresponds to RC2-3 stim pair recorded at electrode 1 in subject 2

9. Figure 7 (BPCs) and all of Figure 8. Plot stim sites, colored by consensus BPCs, on subject inflated pial surfaces; Also plot BPCs colored and ordered by consensus identity
- Run plotInflatedBPCs.m

10. Figure 7 MRI slices. Plot stim sites by consensus BPC color on subject T1 MRI slices
- Run plotBPCs2Slice.m. Manually change sub and ch variables and click within hippocampus on coronal view. The (z, x) coordinates plotted for subjects 1 through 4 are (-15.0, -23.0), (-15.0, 23.0), (-13.0, -24.0), and (-13.0, 26.5), respectively.

11. Figure 9. Create Basis profile spectrogram (BPS) outputs and Adjusted Rand Index outputs quantifying similarity to subject BPC outputs
- Run overlapWithBPS

12. Figure 2-1. Save monopoly subject-level BPC outputs without common average rereferencing.
- Run BPCsNoCAR
- To plot example of channels that go into construction of the common average reference, run figRev1_CARex.m.

13. Comparison of line noise removal methods on smearing of stimulation artifact: spectrum interpolation method used in this paper vs Butterworth notch filter
- Run figRev2_LNremoval

14. Analysis of distance between stimulation site and measurement electrode vs. consensus BPC assignment or BPC SNR
- Run distanceVBPC

