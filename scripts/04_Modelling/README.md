# Modelling

This folder contains scripts to perform modelling of Pol II occupancy, Pol II turnover, and transcription initiation.

## Order of Scripts

1. [01_data_prep_for_modelling.R](/scripts/04_Modelling/01_data_prep_for_modelling.R) 
   * Prepares SMF-derived Pol II occupancy and PRO-seq-derived Pol II turnover data.

2. [02_modelling_Fig4_Fig5F.m](/scripts/04_Modelling/02_modelling_Fig4_Fig5F.m) 
   * Performs modelling of Pol II occupancy (q), Pol II turnover rate (k~t~), and transcription initiation rate (k~i~). This corresponds to Figure 4B and 5F.

3. [03_Revision_visualising_PolII_occupancy_and_initiation_time.R](/scripts/04_Modelling/03_Revision_visualising_PolII_occupancy_and_initiation_time.R)  
   * Visualises the distribution of Pol II occupancy (q) and transcription initiation time (k~i~). This corresponds to Figure EV5.