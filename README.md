# DRT-Peak-Analysis-
Automated extraction of R, Tau, Capacitance, and peak frequency from DRTtools distribution curves.

IMPORTANT: This tool is designed to work with outputs from DRTtools.
Link: https://github.com/ciuccislab/DRTtools

STEPS TO USE:
1. Perform deconvolution on γ(lnτ)-τ curve in the 'Peak Analysis' section of DRTtools.
2. Save the resulting figure using 'Export Results' in DRTtools.
3. Use this module to 'Load a single .fig file' OR a 'folder of .fig files'.

 Features:
 - Automated extraction of R, Capacitance, and peak tau/frequency from DRTtools distribution curves. 
 - Extract all Line/Scatter series (XData/YData)
 - Preview curves (only the last loaded file's data is previewed)
 - Save CSV: Creates individual CSV files for the full data (all curves).
 - Save Parameters CSV: Creates CSV files for peak parameters (R, Tau, Capacitance, Frequency).

 Developed by Masood Fakouri Hasanabadi
 Note: Any use of this code requires proper citation of the related publication of author Masood Fakouri Hasanabadi.

This repository includes a copy of DRTtools (developed by Ciucci Lab) for user convenience. DRTtools is licensed under the LGPL-3.0 License. Original source: https://github.com/ciuccislab/DRTtools
