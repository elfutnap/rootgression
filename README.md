# Rootgression

ROOTgression is a .C macro for CERN's ROOT6 for Multiple Linear Regression calculation over an experimental dataset.
Although easy to use, the macro is able to obtain all the main information related to the datasets that are submitted to it.

## Long story short
ROOTgression is capable of performing:
 - Analysis of replicas (Mean, Variance, Experimental error, Coefficient of variation).
 - Analysis of spatial distribution for experiments of which were analyzed up to three factors.
 - Multiple Linear Regression for up to 4 factors (up to 16 coefficients).
 - ANOVA  (Model tests of linearity and significance with F-test, Sums of Squares,...)
 - Model parameters R^2, R^2 ajusted
 - Residual analysis and normal probability plot.
 
## How do I use it?
Open terminal and reach the directory of Rootgression and the input database. 
Execute ROOT6 by typing "root".
Type: ".x rootgression("myinput.ext")", where "myinput.ext" is your input database. 
Follow the instructions of the macro.

## Which input should I have?
The input dataset should be a tabulation-separated dataset. It can cointain up to 5 columns where the first 4 represents the value of the 4 factors and the last one the answer Y. 
If in the dataset there are 4 columns and you select only two factors to study, the third column will be automatically recognized as the answer column.

## May I see an example?
Of course. In the folder you can try the two given datasets "experimental_data.dat", "factor_design.dat" and "factor_design_2.dat".

For "experimental_data.dat" there are 3 replicas of the centre of the experimental domain, which values are 25, 28 and 29. Minimum for all factor is 0 and maximum is 50.

For "factor_design.dat", there are 5 replicas of the centre of the experimental domain, which values are 42, 41, 41, 39 and 40. Minimum for factor 1 is 0 and maximum 40. Minimum for factor 2 is 1 and maximum 10. Minimum for factor 3 is 5 and maximum 100. As you'll see in the spatial distribution graph, this is a 2-level factor design.

For "factor_design_2.dat", factor coordinates are already autoscaled, so for each factor the minimum is -1 and the maximum is 1. There are 4 replicas at the centre of the experimental design, which are 17.11, 16.05, 17.50, 18.00.
