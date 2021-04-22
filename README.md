CRISPRclassify <a href='https://github.com/CRISPRlab/CRISPRclassify/'><img src='crisprclassify_hex.png' align="right" height="139" /></a>
=================================================================================================================
<!--<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/logo.png" width="400"> -->

CRISPRclassify identifies loci from assembled genomic and metagenomic files (.fasta), and uses a *cas*-independent classification approach to predict the subtype of each locus. CRISPRclassify can be used via web-interface, or executed directly from the command line.

<!-- For more information, please see paper: [CRISPR Visualizer: rapid identification and visualization of CRISPR loci via an automated high-throughput processing pipeline](https://doi.org/10.1080/15476286.2018.1493332). -->

# Installation

CRISPRclassify is implemented in R and can be run in RStudio or in an R session on the command line. For **MacOS** and **Linux**, just ensure you have met the R and Java dependencies below. For **Windows** users, you can use a tool like [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10) or another command line solution to get access to the Linux command line.

### Software requirements:
- #### **R** >= version 3.6.3

`R --version`
```
R version 4.0.2 (2020-06-22) -- "Taking Off Again"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin17.0 (64-bit)
```

  There are many tutorials for installing or updating [R](https://www.r-project.org/) to the minimum required version. One that covers Mac, Linux, and Windows is [listed here](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu).
<br/>

- #### **Java** >= version 1.8

`Java -version`
```
java version "1.8.0_202-ea"
Java(TM) SE Runtime Environment (build 1.8.0_202-ea-b03)
Java HotSpot(TM) 64-Bit Server VM (build 25.202-b03, mixed mode)
```


  Mac users -> [Install Java on MacOS](https://opensource.com/article/20/7/install-java-mac)

  Linux users -> [Install Java on Linux](https://opensource.com/article/19/11/install-java-linux)


---
### Install Directly from GitHub (Recommended):
1. Open the RStudio console, or start an R session on the command line by typing `R`.

2. Install and load the devtools package, which enables installing packages from GitHub:

  `install.packages("devtools")`

  `library(devtools)`

3. Install the CRISPRclassify package:

  `install_github("CRISPRlab/CRISPRclassify")`

  This command will detect any packages with more recent versions available, and display the prompt below:

  <img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/packageUpdates.png" width="400">

  It is recommended to at least update the following packages if listed: **dplyr, readr, magrittr, stringr, stringi, DT, xgboost, memoise, shiny, stringdist, ggplot2, and tidyr.** The installation will fail if one of the required packages does not meet the minimum version requirements:

  ```
  Error: package ‘tidyr’ 1.0.0 was found, but >= 1.1.2 is required by ‘CRISPRclassify’
  Execution halted
  ```

  To update a specific package (tidyr, for example), run:  `install.packages('tidyr')`. Run `install_github("CRISPRlab/CRISPRclassify")` again once all required package versions have been updated.

4. Load the package and launch the app!

  `library(CRISPRclassify)`

  `CRISPRclassify::launchApp()`


Running launchApp() will automatically launch CRISPRclassify in a local web-browser window, ready for use!

**Tech Note:** The port is randomly assigned at runtime. If you close the page or are just curious about the URL, you can see the port listed in the command line output:

```
Listening on http://127.0.0.1:6923
```
<br/>

### Or Download Package from GitHub:
Alternatively, you can download the package to your local machine, and install from there.
1. Download the **.zip** from GitHub.

2. Extract the .zip file to a folder, for example, on your Desktop. The folder will be named **CRISPRclassify-master**.

3. Install the package via command line:

  `cd ~/Desktop` (assuming the extracted folder is on your desktop)

  `R CMD INSTALL CRISPRclassify-master`

  The installation will fail if one of the required packages does not meet the minimum version requirements:

  ```
  Error: package ‘tidyr’ 1.0.0 was found, but >= 1.1.2 is required by ‘CRISPRclassify’
  Execution halted
  ```

  To update a specific package (tidyr, for example), run:  `install.packages('tidyr')` in RStudio or in an R session (start a new R session via command line by simply running the command `R`). Run `R CMD INSTALL CRISPRclassify-master` on the command line again once all required package versions have been updated.

4. Open the RStudio console, or an R session on the command line by typing `R`.

5. Load the CRISPRclassify package and launch the app:

  `library(CRISPRclassify)`

  `CRISPRclassify::launchApp()`

<br/>
Once the app is running, you will see this page:

<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/start.png" width="800">


---
# Command Line Methods
To run analyses without going through the UI, you can use the following methods. These methods can be executed directly on the command line, just make sure to update the path to your .fasta or repeat.txt file. **Note:** if you are executing this command in the same directory as the file you want to classify, you only need to include the name of the file. However, we recommend providing the full file path to avoid '**File Not Found**' errors. Both examples are illustrated below.

#### Output:

Both **classifyFasta()** and **classifyRepeats()** output a **.crclass** file with classification results, located in the same directory as the file that was classified.

## Classify FASTA files with classifyFasta()
Extract and classify repeats from a genome in .fasta format:
#### Command line:
`Rscript -e 'library(CRISPRclassify); CRISPRclassify::classifyFasta("/Users/yourusername/Desktop/example_DGCC7710.fasta")'`

or using only the file name without the path:

`Rscript -e 'library(CRISPRclassify); CRISPRclassify::classifyFasta("example_DGCC7710.fasta")'`

#### RStudio (requires full file path):
`library(CRISPRclassify)`

`CRISPRclassify::classifyFasta("/Users/yourusername/Desktop/example_DGCC7710.fasta")`



## Classify repeats directly with classifyRepeats()
To use repeats derived from another CRISPR detection tool, first format them into a **.txt** or **.csv** file, with one repeat per line, without spaces or punctuation marks (commas, quotes, etc.). The file should look like the this in a text editor:

<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/repeat_format.png" width="300">

#### Command line:
`Rscript -e 'library(CRISPRclassify); CRISPRclassify::classifyRepeats("/Users/yourusername/Desktop/test_repeats.txt")'`

or using only the file name without the path:

`Rscript -e 'library(CRISPRclassify); CRISPRclassify::classifyRepeats("test_repeats.txt")'`

#### RStudio (requires full file path):
`library(CRISPRclassify)`

`CRISPRclassify::classifyRepeats("/Users/yourusername/Desktop/test_repeats.txt")`



---
# An Example Analysis

#### Importing a file:
Included in the **/example_files** directory is the genome of *Streptococcus thermophilus* DGCC 7710 to use as an example. These contigs are also available via NCBI: NZ_AWVZ01000001.1, NZ_AWVZ01000002.1, and NZ_AWVZ01000003.1. Click the **Browse** button in CRISPRclassify and select example_DGCC7710.fasta, then click Classify. The results should match the image below:

<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/example.png" width="800">


The calculated results display distinct repeats from each of the four loci. Loci 1 and 2 are predicted to be subtype II-A, while locus 3 is classified as III-A and locus 4 is predicted to be subtype I-E. The distinct k-mers that significantly contribute to the classification of each subtype are highlighted in blue (forward) and yellow (reverse complement). If we did not know the genus or species that these loci belonged to (for example, processing a metagenomic data set with potentially several genera), the **Closest Strain** column can give us some idea of what organism we're working with.

---

# Output
<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/columns.png" width="800">
The grid output contains the following columns:

 - **Locus:**  Numeric designation for each detected CRISPR array
 - **Contig:** Name of the contig where the CRISPR array is found in the source .fasta genome file
 - **Range:** Location in the contig that contains the entire detected array
 - **Subtype:** Predicted subtype for this locus
 - **Repeat:** Lists all distinct repeats for each locus (duplicates are excluded). Blue and yellow highlights in the repeat sequence indicate significant k-mers that were used for prediction of the subtype
 - **Repeat Count:** How many of this particular repeat are found at this locus (since duplicates are excluded from the grid)
 - **Probability:** The probability that this locus is the predicted subtype. Higher is better. Predictions are considered optimal when the probability is .85 or higher.
 - **Closest Strain:** Lists the strain with the closest matching repeat sequence from the original training set. Useful when working with unknown sequences in metagenomes.
 - **Edit Dist:** The Levenshtein distance calculated between the detected repeat and the corresponding repeat from the closest strain. An edit distance of 0 means a perfect match. An edit distance of 3, for example, indicates 3 nucleotides are mismatched between the repeats. Any matching repeats with more than 20 mismatched nucleotides returns NA.


 More information on the closest matching strain (in terms of repeat similarity) from the training data can be found in the **Download** file after hitting the **Download** button:
  - Name of the strain with the closest matching repeat sequence
  - Accession of the matching strain
  - Location of matching repeat locus in genome
  - CRISPR subtype of the matching repeat locus, previously determined by [Makarova et al.](https://pubmed.ncbi.nlm.nih.gov/31857715/)
  - Matching repeat sequence

---

# Pro tips
**Tech Note:** the current size limit is set limit uploaded file size to 200GB (gigabytes, not bases). If you wish to process a file larger than this, you'll need to download the package from GitHub, extract it to a folder, then set the GB_size flag in **app.R** to a higher number that will suit your needs. Then remove, reinstall, and reload the package.
```
GB_size <- 200
```
```
{in RStudio console or R prompt}
remove.packages('CRISPRclassify')

{on the command line}
R CMD INSTALL CRISPRclassify-master

{in RStudio console or R prompt}
library(CRISPRclassify)
CRISPRclassify::launchApp()
```
---
# Troubleshooting
#### 1) Error: dependencies ‘httr’, ‘rvest’, ‘xml2’ are not available for package ‘tidyverse’. installation of package ‘tidyverse’ had non-zero exit status
 - Tidyverse has some non-R package dependencies required in order to function correctly. In Linux, you can run the following to ensure they are installed:
 ```
 sudo apt install libcurl4-openssl-dev libssl-dev libxml2-dev
 ```
 Additionally, if you see
 `could not find function ggplot; str_c; arrange;`, this is an indicator that the tidvyverse package was not installed correctly.

#### 2) Error in readRDS(file) :  cannot read workspace version 3 written by R 3.6.3; need R 3.5.0 or newer
  - Ensure your version of R is >= 3.6.3. Instructions on how to do so are [listed here](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu).

#### 3) Generic 'could not find function ...'
 - The dependencies for this project include **shiny, tidyverse, DT, xgboost, memoise and stringdist**. Ensure each of these packages is installed. To check the package version, in R run:

 `packageVersion("packageName")`

 Ex.
```
packageVersion("DT")
[1] ‘0.15’
```

  To see all attached and loaded packages, run:

  `sessionInfo()`
