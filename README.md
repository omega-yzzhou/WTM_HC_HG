WTM_HC_HG
Wuhan University Tropospheric Model - Height Correction for Horizontal Gradients

Overview
WTM_HC_HG is a professional MATLAB toolbox designed for tropospheric horizontal gradients processing and height correction. Developed at the GNSS Research Center of Wuhan University, 
this tool implements the Wuhan University Tropospheric Model (WTM) to bridge the vertical gap between grid-wise product's reference orography and target point. When using external grid
products (like GRAD), height differences between the target point and the grid often lead to significant errors. This tool provides a first-order exponential height correction to ensure 
high-precision space geodetic data processing and aircraft navigation.

Key Features
1. Precise Height Correction: Uses a specialized exponential model to adjust gradients for height differences.
2. Seamless Integration: Compatible with TU Wien 1.0°×1.0° GRAD products.
3. Global Support: Includes annual and semiannual time variation terms.

Installation
1. Clone the repository: git clone https://github.com/YourUsername/WTM_HC_HG.git.
2. Required Files: Ensure the following grid coefficient file (WTM_HC_HG.MOD) is present in your project directory.

Quick Start
1. Prepare Coordinates: Create an Example.BLH file (Format: station name, Lat_deg, Lon_deg, Hgt_m).
2. Run the Main Script: Open Cal_GRAD_HG_HC.m and configure your time range.
3. Execute: Run the script in MATLAB to generate your corrected .HG results.

Repository Structure
1. Code
--Cal_GRAD_HG_HC.m: Main entry point for processing.
--grad_grid.m: Horizontal gradient calculation by grid-wise GRAD.
--wtm_hc_hg.m: Height correction coefficient calculation.
--read_blhfile.m: Read .BLH coordinate file.
--modified_julday.m: Convert date to mjd.
2. input and output files
--Example.BLH: Example coordinate file.
--Example_Result.HG: Example result file.
3. Model and product files
--WTM_HC_HG.MOD: Height correction model file.
--2019: Directory for example grid-wise GRAD products.

Contact
Yaozong Zhou (周要宗), Postdoc, GNSS Research Center, Wuhan University
Email: zhouyaozong@whu.edu.cn
