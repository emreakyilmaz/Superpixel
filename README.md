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
SREP: similarity ratio with Euclidean distance 
SRMP: similarity ratio with Mahalanobis distance
SRAMP: Similarity Ratio with Mahalanobis distance with adaptive scheme
