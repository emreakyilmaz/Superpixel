# Superpixel
Superpixel segmentation as a preprocessing step is an oversegmentation technique that groups similar neighboring pixels into
regularly organized segments with approximately the same size. As boundaries of the objects are important elements to be traced, superpixels should adhere well to
the edges. This can only be achieved by an algorithm robust to speckle noise. 

In this work, similarity ratio is first developed as a new metric  that is robust to
speckle noise. Secondly, Mahalanobis distance is used instead of Euclidian so that
the superpixel can fit better to shapes in the real world. Thirdly, the constant
determining the relative importance of radiometric and geometric terms is replaced
with an adaptive function. 

## Algorithms
1. SREP: similarity ratio with Euclidean distance 
2. SRMP: similarity ratio with Mahalanobis distance
3. SRAMP: Similarity Ratio with Mahalanobis distance with adaptive scheme

## Similarity Ratio
For superpixelling in SAR images, A similarity ratio metric as a more robust approach is proposed.
This metric has its basis on likelihood ratio. The likelihood ratio can be derived by assuming that two samples are drawn from
a population.

Let say a sample X is N(μ1 , σ1^2), and another sample Y is N(μ2 ,σ2^2). Let (x1, …, xm) be a sample space of X and (y1, …, yn) be a sample space of
Y. If these two samples are drawn from population Z whose random space is (x1,…, xm, y1, …, yn), then Z is N(μ , σ2
) as shown in Figure 1:

![Figure_1] (./Figures/Figure_1.PNG)

'Figure 1:' Illustration of Sample X, Y and Population Z


Moreover, if all random variables are assumed mutually independent then the
likelihood functions for X, Y and Z would be as:

![Figure_2] (./Figures/Figure_2.PNG)

The maximum likelihood estimates can be computed from the above likelihood
functions as:

![Figure_3] (./Figures/Figure_3.PNG)

The likelihood functions are simplified by using the maximum likelihood
estimates as:

![Figure_4] (./Figures/Figure_4.PNG)

Since the two samples, X and Y, are drawn from the population Z, then their
likelihood functions can be related as:

![Figure_5] (./Figures/Figure_5.PNG)

Substituting this equation into the likelihood formulation and after some algebraic manupilations similarity ratio can be
obtained as:

![Figure_6] (./Figures/Figure_6.PNG)

This ratio converges to unity from positive infinity as similarity between the means increases. Hence, it is a measure of the amount of deviation between the
intensities of the superpixels, i.e., it is a measure of dissimilarity.