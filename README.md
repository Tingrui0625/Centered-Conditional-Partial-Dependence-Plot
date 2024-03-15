This code is for the master thesis of ccPDP only.

"1. Tuning and Model Selection.R" includes tuning NN and comparison between NN and other baseline models. "2. ccPDP functions.R" gives the 5 ccPDP calculation functions. The functions for Mplot, Loss to ALE are written here additionally. "3.1 Tuning ccPDP to ALE&PDP.R" uses bbotk package to tune the hyperparameters in the raised ccPDPs in multisession. 3.2 file tunes the hyperparameter regarding model fidelity loss in multisession. "4. Do once functions.R" records the ccPDP calculation for one example. It is prepared for multisession afterward in "5. Computing PDP". "6. Creation for mplot.R" outputs the mplot for all the examples in the thesis. "7. Plotting&LOSS.R" contains the code to run the final plots with the CSV files from the previous code. 8.1 applies our methods in the Wine quality dataset while 8.2 in the KC Housing dataset. "99.test" is for flexible testing tasks.

Here the data generation function in "5. Computing ccPDP.R" is used multiple times for in other files.

The folder "ccPDP hyperparameters" contains the tuned hyperparameters with different goals. The resulting ccPDPs are in folder "ccPDP values".



