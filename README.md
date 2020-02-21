BetaMBC
-----

This is the git repostitory for my Masters Thesis project.
The code is written in R and features:
- Model based clustering
- Finite mixtures of unimodal beta distributions
- EM Algorithm
- S3 OOP

ABSTRACT
----
Diferentially methylated DNA regions are valuable biomarkers of human disease
with potential in the feld of personalised medicine. There are three common methylation
profiles described in literature: hypomethylation, hemimethylation and hypomethylation.
Currently, their detection is based on simple univariate parametric
tests that require multiple testing corrections and data transformations. Due to
the statistically undesirable features of the raw DNA methylation data, alternative
approaches to the problem that do not call for data transformation are limited. Furthermore,
DNA methylation data analysis is a high-dimensional "big data" problem
which requires a computationally efficient solution.
The aim of this work was to investigate the potential of model based clustering
of DNA methylation data to effectively and efficently differentiate hypomethylated,
hemimethylated and hypermethylated CpG sites between two sets of samples using a
fnite mixture of unimodal densities and the EM algorithm for parameter estimation.
We present a computationally eficient and novel model based clustering approach
to differential methylation analysis and test it on real DNA methylation
data. We compare our method with the conventional testing approach and shared
kernel Bayesian screening. Our approach was able to identify patterns in the data
congruent with biological theory and shows promise as a supplementary screening
tool for differential methylation region detection. Furthermore, it introduces the use
of model based clustering as a viable tool for this type of DNA methylation data
analysis.



