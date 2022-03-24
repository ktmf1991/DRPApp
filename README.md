# DRPApp
A Shiny app designed to perform simple and quick analysis and comparison of drug response profiling on multiple samples

# Pre-requisite
- If you are using the app on <ins>**Windows**</ins>, please install *rtools* following the instruction [here](https://cran.r-project.org/bin/windows/Rtools/rtools40.html).
- Have *shiny* installed first.
- Prepare input files with the same format as the example file in the repository. Each row is a description of a well in the drug assay plate and each column describes the following:
  - Well: The well on the plate
  - Treatment_Drug: Drug used on this well
  - Concentration_1	and Unit_1: The concentration of drug used and its unit. This column is not used in the analysis but for users to read the results with ease. For example, you can have 10 uM instead of 10000 nM in these two columns.
  - Concentration_2	and Unit_2: The concentration of drug used and its unit, but used in the analysis. Therefore, the **UNITS for the same drug MUST BE the SAME**.
  - Remark: Remark for anything, except for the **vehicle control**, which must be **"neg_control"**.
  - Count: Cell count, or any readout for your assay.
  - Patient: Patient, or other sample. Recommended to use integer numbers.
  - Plate: Plate. This app is designed for medium-sized chemical library and assume for drugs are fixed on the same plate for the same patient tested.

# To use the app:
1. Call *shiny* library. Run `runGitHub("DRPApp","ktmf1991","main")`.
2. Upload data by clicking *Browse*.
3. Click *Analyse*. It may take a minute depending on the speed of your computer.
4. Check the results.

# How to check my results?
(To be updated)
