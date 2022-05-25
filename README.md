# DiagTestCluster
Calculating diagnostic test values in the presence of clustering.  

Advanced screening and diagnostics for binary-scale diagnostic tests  
Scripts and example datasets to calculate binary-scale diagnostic test characteristics in the context of clustering  
based on Genders et al. (2012); Hujoel, Moulton, & Loesche (1990); Kirkwood & Sterne (2003); McDonald (2019);   
Mercaldo et al. (2007); Williams (2000); and Ying et al. (2020).

## Input: data structure in wide format (i.e., one row per patient)  
id = patient identification number (identifies the clusters i =1…I)  
TP = number of true-positive observations  
FN = number of false-negative observations  
FP = number of false-positive observations  
TN = number of true-negative observations  
Test_type (optional) = indicator variable if multiple index tests are considered (eg, 0 = test A, 1 = test B)  

____________________________							

## Results generated from input:
### Patient-level analyses							
							
Patient-level contingency table							

Prevalence-independent measures:   	  						
	Type_of_method,	Sensitivity (%)	Lower_CI	Higher_CI,	Specificity (%)	Lower_CI	Higher_CI
												
Prevalence-dependent measures:   	  						
	Type_of_method,	PPV (%)	Lower_CI	Higher_CI,	NPV (%)	Lower_CI	Higher_CI
							
Patient-level likelihood ratios							
							
### Segment-level analyses							
							
Segment-level contingency table							
							
Prevalence-independent measures:	     						
	Type_of_method,	Sensitivity (%)	Lower_CI	Higher_CI,	Specificity (%)	Lower_CI	Higher_CI

Prevalence-dependent measures	     						
	Type_of_method,	PPV (%)	Lower_CI	Higher_CI,	NPV (%)	Lower_CI	Higher_CI
							
Segment-level likelihood ratios (NB: CI's not adjusted for clustering!)							

Intracluster correlation coefficients	

____________________________

## Visualizations
Forest plots: Allow for comparison of outcomes obtained by each method. 
____________________________						

Some remarks: 
The methods reported by Genders et al. (2012) are validated using their dataset and supplementary materials.  
They report slightly different CI's in the main text vs. in the supplementary materials, the latter of which   
are exactly reproduced by the R script. CI's from the mixed effects logistic regression (called logistic   
random-effects model by Genders et al., 2012) are slightly different due to the default number of quadrature  
points when approximating the integral over the random effects structure (1 in R's glmer vs. 7 in STATA's xtmelogit).

When sourcing the script, use the "print.eval=TRUE" option so that the forest plots can be generated, e.g.: source("~path/Advanced_screening_diagnostics_Tamasi_website.R", print.eval=TRUE)

  

