# DiagTestCluster
Calculating diagnostic test values in the presence of clustering.  

Advanced screening and diagnostics for binary-scale diagnostic tests  
Scripts and example datasets to calculate binary-scale diagnostic test characteristics in the context of clustering  
based on Genders et al. (2012); Kirkwood & Sterne (2003); Hujoel, Moulton, & Loesche (1990), Ying et al. (2020), etc.  

# Input: data structure in wide format (i.e., one row per patient)  
id = patient identification number (identifies the clusters i =1…I)  
TP = number of true-positive observations  
FN = number of false-negative observations  
FP = number of false-positive observations  
TN = number of true-negative observations  
Test_type (optional) = indicator variable if multiple index tests are considered (eg, 0 = test A, 1 = test B)  

# Data generated from input:  
j = number of available observations per cluster  
D = number of diseased observations  
N = number of nondiseased observations  
Test = index test result (e.g., Genders (2012): CT angiography: 1 if positive, 0 if negative)  
Inv_test = inverse result of test for specificity calculations  
Disease = presence of disease according to reference test (conventional angiography: 1 if positive, 0 if negative)  
