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
Of course. In the folder you can try the given datasets.

## What is rootgression_quick?
Rootgression performs a rigorous matrix-based MLR which is optimized for experiments designed with the DOE techniques, such as factor designs. Parameters of the linear model are or aren't included in the model based on the standard error and the number of replicas performed at the centre of the experimental design. The calculated model considers every interaction between the factors, so the full model of a 2 factor design would be: 

y = a + b1*f1 + b2*f2 + b3*f1*f2

Rootgression quick is designed in order to be applied to every experimental situation which does not consider the experiment distribution, such as poll answers. The linear model is also calculated without considering the interaction between factors, so a 2-factor experiment will have the following shape:

y = a + b1*f1 + b2*f2

#### smoke.dat
This dataset has 4 factors and 100 experiments. It is an investigation of the aversion toward the smoke based on 4 explored parameters (exposition, risk, figure, age). This dataset is suggested to be analyzed with rootgression_quick because, being random samples and not based on an experimental design, there aren't replicas of the samples on which calculate the crytical sb that gives significativity on the model parameters. 
Execute the example just by copying smoke.dat from the sample folder to the same folder of rootgression_simplified.dat and input ".x rootgression_simplified.c ("smoke.dat"). Requested parameters will be the number of factors (4) and the number of the samples (100).

#### factor_design.dat
In this dataset there are 3 factors and 8 experiments. There are 5 replicas of the centre of the experimental domain, which y-values are 42, 41, 41, 39 and 40. As you'll see in the spatial distribution graph, this is a 2-level factor design.
