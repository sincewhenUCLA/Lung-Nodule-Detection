# Lung-Noudle-Detection

This the extra-low-dose CT lung nodule detection demo code developed by "Shiwen Shen" shiwenshen@ucla.edu for CDSC project.

This file introduces the workflow and usage of the lung nodule detection pipeline.


######################################################
#
#                     Workflow
######################################################
On-line detection task:


Input: 		CT lung image stacks (Analyze file format)
Output: 	Lung nodule binary mask which could be mapped directly to the original images(ANALYZE 7.5 format)
1 segmentation (see folder segmentation): 	segment the initial nodule candidates from CT images
2 preselection (see folder preselection): 	reduce the false positive rate based on pre-defined rules
3 feature extraction (see folder feature extraction): generate 27 features for each nodule candidates
4 classification (see mainNoduleDetection.m file): output the final nodule using trained classifier

**PLEASE NOTE** Depending on the size of the input image, the pipeline may require a large amount (> 4Gb) of memory to
complete its task. As a result, the computer may freeze until the pipeline has completed its execution. We recommend that
you run the pipeline on a dedicated development environment.

#################################################
#
#                     Usage
#################################################

###############################
main lung nodule detection task

Step 1:
Add the folder of current code and its subfolders to path

Step 2:
Open "mainNoduleDetectionAnalyzeFile.m" file, this is the main function for lung nodule detection.

Step 3:
Run this main function and 3D nodule mask will be automatically shown and you can view it by scrolling the mouse 


Author:	 Shiwen Shen
Date: 	 09/28/2014
Email: 	 shiwenshen@ucla.edu
Copyright:  Medical Imaging Informatics Group, UCLA
