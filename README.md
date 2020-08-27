# sepsis_and_cancer
Code and metadata associated with project: 
## Sepsis-associated pathways segregate cancer groups 
Himanshu Tripathi, Samanwoy Mukhopadhyay &amp; Saroj Kant Mohapatra
## Summary
Transcriptome-based systems biology approach segregates cancer into two groups ((termed Sepsis-Like Cancer, or SLC and Cancer Alone, or CA) based on similarity with SS. Host response to infection plays a key role in pathogenesis of SS and SLC. However, we hypothesize that some component of the host response is protective in both SS and SLC.
## How to set up data and code for analysis
### Step 1: Installation of the three Data Packages ssnibmg, tcnibmg and tcnibmgML
- Download the following data packages files from https://figshare.com/ using the links provided
below:
- ssnibmg_1.0.tar.gz from https://doi.org/10.6084/m9.figshare.4592800.v1
- ssgeosurv_1.0.tar.gz, tcnibmg_1.0.tar.gz and tcnibmgML_1.0.tar.gz from
https://doi.org/10.6084/m9.figshare.8118413.v3
- Change the directory to where you saved the files. Start R.
- At the R prompt, issue the following commands:
> install.packages(pkgs="ssnibmg_1.0.tar.gz", repos=NULL)
> install.packages(pkgs="ssgeosurv_1.0.tar.gz", repos=NULL)
> install.packages(pkgs="tcnibmg_1.0.tar.gz", repos=NULL)
> install.packages(pkgs="tcnibmgML_1.0.tar.gz", repos=NULL)
- Now the three data packages are installed on your computer.
- Check with the following commands:
> library("ssnibmg")
> library("ssgeosurv")
> library("tcnibmg")
> library("tcnibmgML")
### Step 2: Running the analysis
- It is assumed that you have access to a folder tcnibmgdoc
- Start R and set the working directory to tcnibmgdoc
- Run the script main.R
