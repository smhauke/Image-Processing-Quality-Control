# Image-Processing-Quality-Control

This repository contains an SOP for image pre-processing and some Python scripts for quality control. 

Here is a detailed breakdown of each item in the repo:
Image Processing Pipeline - V1.0.docx: An SOP detailing image pre-processing steps prior to image analysis in Imaris. 
Running Python Scripts for Image Processing and Analysis - V1.0.docx: A basic rundown of running Python scripts for the average biologist/non-programmer. Frequently referenced by the Image Processing Pipeline document.
Microglia_2D_Preprocessing.ijm: Example Fiji macro. Takes .czi image (from Zeiss microscope) as input and outputs a .tiff file.
RawImageQC.py: First round of quality control. Reads raw .czi file and flags differences in microscope settings.
PostProcessingQC.py: Second round of quality control. Reads processed .tiff files
tiff_to_ims.py: Final pre-processing step. Manually converts .tiff images to .ims images that can be read by Imaris. This bypasses the buggy Imaris TIFF reader in Imaris v9.
ims_integrity_validator.py: Scans through files open in the Imaris Arena (found in a folder within ProgramData) and flags abnormalities.
