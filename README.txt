README - Code and Data for "Towards Robust Matched Observational Studies with
         General Treatment Types: Consistency, Efficiency, and Adaptivity"

Authors: Siyu Heng, Elaine K. Chiu, and Hyunseung Kang


------------------------------------------------------------
DESCRIPTION
------------------------------------------------------------

This repository contains the R code and the post-matching data set needed to
reproduce every numerical result reported in the paper "Towards Robust Matched
Observational Studies with General Treatment Types: Consistency, Efficiency,
and Adaptivity", including all figures and tables in the main text and in the
supplementary material.

The material is organized into two parts. The /simulations/ folder contains the
simulation studies underlying Figure 2 and Tables 2-4 of the main text and
Tables S.1-S.12 of the supplement. The /data_analysis/ folder contains the real
data analysis of the 2011-2012 National Health and Nutrition Examination Survey
(NHANES), underlying Figure 3 and Table 5 of the main text and Table S.13 of
the supplement. A single file of shared routines,
matched_sensitivity_analysis_functions.R, is used throughout.

The EXHIBIT INDEX below states, for each figure and table in the paper, which
script produces it, where in that script it is produced, and how the output is
returned.


------------------------------------------------------------
REQUIREMENTS
------------------------------------------------------------

- R version 4.6.0 or higher (tested on 4.6.0)

- Required R packages:
    - sensitivitymv
    - designmatch
    - MASS
    - dplyr
    - haven
    - moments
    - mvtnorm
    - ggplot2

- Required R code file, sourced by every script:
    - matched_sensitivity_analysis_functions.R


------------------------------------------------------------
FOLDER STRUCTURE
------------------------------------------------------------

Repository root

    matched_sensitivity_analysis_functions.R
        Shared routines used by both the simulations and the data analysis:
        the conversion between the generalized sensitivity parameter on the 
        Gamma and gamma scale, the rank-based and U-statistic test statistics
        and their dose-weighted counterparts, the worst-case p-value under the
        generalized Rosenbaum sensitivity bounds, the adaptive testing
        procedure, the generalized design sensitivity, and the data-generating
        functions for the six dose-response curves.

/data_analysis/

    Male_ETS_Lung_Function.R
        Section 5 of the main text and Appendix C.4 of the supplement.
        Part 1 (data cleaning and matching) constructs the matched sample from
        the raw NHANES files; Part 2 (data analysis) reproduces Table S.13,
        Figure 3, and Table 5.

    male_40_matched_data.csv
        The post-matching data set, I = 479 matched pairs of men over 40.
        Each row is a matched pair. The columns low.dose and high.dose hold the
        treatment dose (serum cotinine) of the lower- and higher-dose subject in
        the pair; low.response and high.response hold the corresponding
        FEV1/FVC ratios; the columns prefixed low.dose. and high.dose. hold the
        measured covariates of the two subjects. For example, low.dose.age is
        the age of the subject with the lower serum cotinine in the pair.

/simulations/

    Continuous_Outcome_Generalized_Design_Bahadur_Efficiency.R
        Figure 2, Table 2, Table 3, and Table 4 of the main text, and Tables
        S.1-S.2 of the supplement (binary outcome case). Uses kappa(z) = z and
        phi(z) = z throughout.

    SampleSplittingAndAdaptiveTesting.R
        Tables S.3-S.6 of the supplement: finite-sample power of the adaptive
        test compared with sample splitting, at I = 100, 500, 1000, and 5000.

    Dose_Response_Curve_Transformation_Verification_Comparison.R
        Tables S.7-S.8 of the supplement: verification of the generalized
        design sensitivity formula under the dose transformation kappa, and the
        comparison of design sensitivities across kappa.

    Bahadur_Efficiency_Simulations_with_Transformation.R
        Tables S.9-S.12 of the supplement: generalized Bahadur-Rosenbaum
        relative efficiencies and exact slopes under different kappa. Because
        this script is the most computationally demanding, its results are
        written to CSV rather than printed.


------------------------------------------------------------
DATA AVAILABILITY
------------------------------------------------------------

The real data analysis in Section 5 of the main text uses the 2011-2012 cycle
(cycle "G") of the National Health and Nutrition Examination Survey (NHANES),
conducted by the National Center for Health Statistics, Centers for Disease
Control and Prevention. NHANES public-release files are in the public domain
and may be downloaded by anyone, without registration, application, or a data
use agreement. No component of the analysis relies on restricted-access data.

Cycle landing page:

    https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2011

The following six SAS transport files are required. Each address below serves
the data file; replacing the .xpt extension with .htm returns the codebook,
which documents every variable and its coding.

    DEMO_G.xpt    Demographic Variables and Sample Weights
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/DEMO_G.xpt

    COTNAL_G.xpt  Serum Cotinine and Urinary
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/COTNAL_G.xpt

    SPX_G.xpt     Spirometry - Pre and Post-Bronchodilator
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/SPX_G.xpt

    BMX_G.xpt     Body Measures
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/BMX_G.xpt

    MCQ_G.xpt     Medical Conditions
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/MCQ_G.xpt

    SMQ_G.xpt     Smoking - Cigarette Use
                  https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/SMQ_G.xpt

On the landing page above, these appear under Demographics Data (DEMO_G),
Examination Data (BMX_G, SPX_G), Laboratory Data (COTNAL_G), and Questionnaire
Data (MCQ_G, SMQ_G).

Variables used:

    Treatment dose       LBXCOT, serum cotinine in ng/mL, Box-Cox transformed
                         with lambda approximately -0.02
    Outcome              SPXNFEV1 and SPXNFVC, analysed as the FEV1/FVC ratio
    Matching covariates  RIDRETH1 (race), INDFMPIR (income-to-poverty ratio),
                         DMDHREDU (education), BMXBMI, BMXWAIST, MCQ010 (asthma)
    Study population     RIAGENDR = 1 and RIDAGEYR > 40, giving 1065 men

Downloading the files:

    nhanes.files = c("DEMO_G", "COTNAL_G", "SPX_G", "BMX_G", "MCQ_G", "SMQ_G")
    nhanes.base  = "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2011/DataFiles/"
    for (f in nhanes.files) {
      download.file(url = paste0(nhanes.base, f, ".xpt"),
                    destfile = paste0(f, ".xpt"), mode = "wb")
    }

The files must be placed in the working directory from which
Male_ETS_Lung_Function.R is run.

Derived analysis data:

The matching step in Male_ETS_Lung_Function.R reduces the 1065 eligible men to
I = 479 matched pairs and writes male_40_matched_data.csv. That file is included
in /data_analysis/, so Table 5, Table S.13, and Figure 3 can be reproduced
directly, without downloading any NHANES file and without repeating the
matching step.


------------------------------------------------------------
EXHIBIT INDEX
------------------------------------------------------------

Every figure and table in the paper is listed below with the script that
produces it, the location within that script, and the form of the output.
Line numbers refer to the files as distributed in this repository.

MAIN TEXT

  Figure 1     Illustrative diagram; no computation.

  Table 1      Summary of the scope of existing work; no computation.

  Figure 2     /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 36-139.
               Output: printed to the console as DoseResponseCurve.png.

  Table 2      /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 146-510.
               Output: printed to the console, one block per dose-response
               curve.

  Table 3      /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 513-1140.
               Output: printed to the console.

  Table 4      /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 1144-1730.
               Output: printed to the console.

  Figure 3     /data_analysis/Male_ETS_Lung_Function.R, lines 493-505.
               Output: printed to the console.

  Table 5      /data_analysis/Male_ETS_Lung_Function.R, from line 509-905.
               Output: printed to the console.

SUPPLEMENTARY MATERIAL

  Table S.1    /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 1743-2225 and lines 2305-2407.
               Output: printed to the console.

  Table S.2    /simulations/Continuous_Outcome_Generalized_Design_Bahadur_
               Efficiency.R, lines 1743-2300.
               Output: printed to the console.

  Table S.3    /simulations/SampleSplittingAndAdaptiveTesting.R, lines 522-601.
               Calls run_model_table() at I = 100 for each of the six
               dose-response curves, with planning sample proportions 0.10 and
               0.20. Output: printed to the console.

  Table S.4    /simulations/SampleSplittingAndAdaptiveTesting.R, lines 602-679.
               As above, at I = 500.

  Table S.5    /simulations/SampleSplittingAndAdaptiveTesting.R, lines 680-726.
               As above, at I = 1000.

  Table S.6    /simulations/SampleSplittingAndAdaptiveTesting.R, lines 727-761.
               As above, at I = 5000.

  Table S.7    /simulations/Dose_Response_Curve_Transformation_Verification_
               Comparison.R, lines 1-770.
               Design sensitivity and simulated power under the square and log
               dose-response curves for kappa(z) = z, z^2, and log(z).
               Output: printed to the console.

  Table S.8    /simulations/Dose_Response_Curve_Transformation_Verification_
               Comparison.R, lines 771-1150.
               Design sensitivities under all six dose-response curves and four
               choices of kappa. Output: printed to the console.

  Table S.9    /simulations/Bahadur_Efficiency_Simulations_with_Transformation.R,
               block RUN_PART_A, lines 271-398.
               Output: written to file as
               bahadur_relative_efficiency_verification_long.csv, and available
               in the workspace as the object releff_long.

  Tables       /simulations/Bahadur_Efficiency_Simulations_with_Transformation.R,
  S.10-S.12    block RUN_PART_B, lines 195-270.
               Output: written to file as bahadur_slopes_by_kappa_long.csv, and
               available in the workspace as the object slope_long. The three
               tables in the paper are the three parts of this single output,
               split by dose-response curve.

  Table S.13   /data_analysis/Male_ETS_Lung_Function.R, lines 443-490.
               Output: printed to the console as the object balance.table.


------------------------------------------------------------
HOW TO REPRODUCE THE RESULTS
------------------------------------------------------------

1. Set the working directory to the folder containing the script to be run, so
   that source("matched_sensitivity_analysis_functions.R") resolves. Each
   script sources this file at the top.

2. To reproduce a simulation exhibit, open the script named in the EXHIBIT
   INDEX and run it from the top through the block indicated. The scripts are
   sequential: the shared populations and configuration defined at the top of
   each file are used by all blocks below, so a block should not be run in
   isolation from a fresh session.

3. To reproduce the data analysis, either

     (a) run /data_analysis/Male_ETS_Lung_Function.R from the top, having first
         downloaded the six NHANES files described under DATA AVAILABILITY,
         which repeats the matching and then the analysis; or

     (b) begin from the analysis section, which reads the supplied
         male_40_matched_data.csv and requires no NHANES download.

4. Random seeds are fixed within the scripts, so the reported figures are
   reproduced exactly. The simulations use set.seed(100) for the population
   quantities; the Monte Carlo replications within the finite-sample power and
   minimum-sample-size searches set their seed from the replication index, so
   that the same simulated data sets are reused across the tests being
   compared. Comparisons across tests are therefore reflections of the test 
   performance rather than noise due to data generating process inconsistency.
