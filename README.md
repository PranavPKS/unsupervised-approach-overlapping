# An Unsupervised Approach for Overlapping Cervical Cell Cytoplasm Segmentation

https://goo.gl/dWzqcX


Use:
*  The "CytoScriptExample.m" was used for train set cell segmentation.
    It saves the 'CytoGroundTruth.mat' and 'SegmentationResult.mat' in current folder
        'CytoGroundTruth.mat' : contains the segmented ground truth of cytoplasm
        'SegmentationResult.mat': contains the segmentation results.
*  The "CytoScriptExampleTestingset.m" was used for test set cell segmentation


Function Details:
fstack_mod.m -      takes the path of the image stack and returns the EDF image.
getCytoplasmGT.m -  returns the cytoplasm ground truth
preProcess.m -      prepares the EDF images for processing
ReadImgs.m -        read all images present in a folder
cytoplasmSegmentLSM.m-for cytoplasm segmentation
myotsuSTw1.m -      otsu thresholding with prior probability of one class
cropRegion.m -      crop a rectangular region from image
addBorder.m  -      to add appropriate border in the image
