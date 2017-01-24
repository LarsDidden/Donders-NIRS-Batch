# Donders-NIRS-Batch
Batch created to preprocess and analyse NIRS data using SPM software.

This batch was created to process NIRS-data for the FOCOM project (http://www.focom-project.net/).
The batch offers preproccesing and analyzing options given by the SPM software together with a few additions:

-Scalp Coupling Index (SCI) check to check data quality per channel and (if necessary) exclude channels with low data quality.
-Downsampling
-MARA filtering to remove movement artifacts from the data.
-Calculation of the mean (base-line corrected) NIRS signal for all conditions/channels as extra outcome measure.
-Averaging of the NIRS signal for every trial in every condition, similar to event-related potentials in EEG as extra outcome measure.
-An INFO_file is added to the batch to facilitate selecting options for the batch.
-The NIRS analysis has been automated so the (group)analysis of several subjects can be run with one click (when analysis options for different subjects stay the same).

Manual_DNB.pdf gives more information on all batch files and explains how to use the batch. 

To run the batch, all files in the Donders-NIRS-Batch repository are necessary.
In addition to the DNB files, the following files are needed to run the batch:

-SPM_nirs files (http://bispl.weebly.com/nirs-spm.html#/)
-SPM files (http://www.fil.ion.ucl.ac.uk/spm/software/spm12/)
-Oxy3convert files (only when using .oxy3 files)
-Fieldtrip files (http://www.fieldtriptoolbox.org/)

Lars Didden
Joost Wegman
