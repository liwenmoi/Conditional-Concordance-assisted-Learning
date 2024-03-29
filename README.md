# Conditional-Concordance-assisted-Learning

Executable codes of the paper:

Li W, Li R, Yan Q, Feng Z, Ning J (2023). "Conditional Concordance-assisted Learning for Combining Biomarkers for Cancer Population Screening". Statistics In Medicine.

(1) main.R: R file contains wrapper functions used to implement the proposed method. It shows how to analyze the example dataset "sub.csv" and to validate the performance on "val_sub.csv". This is essentially one run of the simulation studies in the paper.

(2) HelperFunctions.R: helper functions that are used in main.R.

(3) LogLikeExactC.CPP: codes of the CCAF function that are used in main.R. It is written in CPP language to boost computational speed.

(2) sub.csv: example data file that contains disease status "Y", two biomarkers "X1" and "X2", matching group membership "Z", and "groupID" which shows the matched pairs.

(3) full.csv: the full dataset that is used to generate the matched case-control data set "sub.csv". Note that there is no biomarker information in this data.

(4) val_sub.csv: the large simulated data set that serves as the validation data.

(5) val_full.csv: the full dataset corresponds to "val_sub.csv".
