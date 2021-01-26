<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/logo.png" width="400">

CRISPRclassify identifies loci from assembled genomic and metagenomic files (.fasta), and uses a *cas*-independent classification approach to predict the subtype of each locus.

<!-- For more information, please see paper: [CRISPR Visualizer: rapid identification and visualization of CRISPR loci via an automated high-throughput processing pipeline](https://doi.org/10.1080/15476286.2018.1493332). -->

# Installation

CRISPRclassify is implemented in R and is run locally from the command line. For **MacOS** and **Linux**, just ensure you have met the R and Java dependencies below. For **Windows** users, you can use a tool like [WSL](https://docs.microsoft.com/en-us/windows/wsl/install-win10) or another command line solution to get access to the Linux command line.

#### Software requirements:
- **R** >= version 4.0.2

  Many users install R by downloading the [RStudio IDE](https://rstudio.com/). If you prefer a manual R installation, see [here](https://www.r-project.org/).

- **Java** >= version 1.8

  Mac users -> [Install Java on MacOS](https://opensource.com/article/20/7/install-java-mac)

  Linux users -> [Install Java on Linux](https://opensource.com/article/19/11/install-java-linux)


---
#### Download and Run:
1. Download or clone this repository into a folder, for example, on your desktop:

  `cd ~/Desktop`

  `git clone https://github.com/CRISPRlab/CRISPRclassify.git`

2. CD into the newly created folder, and run the start script:

  `cd CRISPRclassify/`

  `./start.r`


The first time the **start.r** script is run, it will automatically install the required R library dependencies. The Shiny library in particular takes some time to download, so it may take a few minutes to complete. Running **start.r** will automatically launch CRISPRclassify in a local web-browser window, ready for use!

**Tech Note:** the server port is automatically assigned at runtime. The browser window will automatically launch with the assigned port in the URL, but if you accidentally close the page or are just curious, you can see the port in the terminal output:

```
Listening on http://127.0.0.1:3014
```

Once the app is running, you will see this page:

<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/start.png" width="800">



---
## An Example Analysis

#### Importing a file:
Included in the **/example_files** directory is the genome of *Streptococcus thermophilus* DGCC 7710 to use as an example. These contigs are also available via NCBI: NZ_AWVZ01000001.1, NZ_AWVZ01000002.1, and NZ_AWVZ01000003.1. Click the **Browse** button in CRISPRclassify and select example_DGCC7710.fasta, then click Classify. The results should match the image below:

<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/example.png" width="800">


The calculated results display distinct repeats from each of the four loci. Loci 1 and 2 are predicted to be subtype II-A, while locus 3 is classified as III-A and locus 4 is predicted to be subtype I-E. The distinct k-mers that significantly contribute to the classification of each subtype are highlighted in blue (forward) and yellow (reverse complement). If we did not know the genus or species that these loci belonged to (for example, processing a metagenomic data set with potentially several genera), the **Closest Strain** column can give us some idea of what organism we're working with.

---

## Output
<img src="https://github.com/CRISPRlab/CRISPRclassify/blob/master/img/columns.png" width="800">
The grid output contains the following columns:

 - **Locus:**  Numeric designation for each detected CRISPR array
 - **Contig:** Name of the contig where the CRISPR array is found in the source .fasta genome file
 - **Range:** Location in the contig that contains the entire detected array
 - **Subtype:** Predicted subtype for this locus
 - **Repeat:** Lists all distinct repeats for each locus (duplicates are excluded). Blue and yellow highlights in the repeat sequence indicate significant k-mers that were used for prediction of the subtype
 - **Repeat Count:** How many of this particular repeat are found at this locus (since duplicates are excluded from the grid)
 - **Probability:** The probability that this locus is the predicted subtype. Higher is better.
 - **Closest Strain:** Lists the strain with the closest matching repeat sequence from the original training set. Useful when working with metagenomes. More information on this strain from the training data can be found in the Download file after hitting the Download button.

---

### Pro tips and Troubleshooting
**Tech Note:** the current size limit is set limit uploaded file size to 135GB (gigabytes, not bases). If you wish to process a file larger than this, you'll need to set the GB_size flag in app.R to a higher number that will suit your needs:
```
GB_size <- 135
```
