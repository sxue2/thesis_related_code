1. 
Model 1-3 are fitted using fullscan_GM(model number).py
2.
Model4 are fitted seperately because each trait is taking different variance components estimate.
Variance component estimate are taking from Data_prep/CheckingFormatsofInputeFilesforGWAS.html
python_traitName.py is used for fitting model 4.
3.
Model 5 is fitted using python_leave_one.py
4.
Simulation for Model 4 are seperate files, each ended with a simulation.py
5. Significance test
For fixed effect model, fixedMarkerTest.py is the file for testing.
For mixed effect model, markerTestGM.py is the file for testing.