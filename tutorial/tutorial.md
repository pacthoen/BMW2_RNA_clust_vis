---
layout: tutorial_hands_on

title: "RNA-Seq data analysis,clustering and visualisation tutorial"
subtopic: introduction
priority: 2

tags:
    - bulk
    - rna-seq
    - collections
    - mouse
    - clustering
    - PCA
level: Introductory
zenodo_link: 
questions:
    - How to identify differentially expressed genes across multiple experimental conditions?
    - What are the biological functions impacted by the differential expression of genes?
    - How to visualise high-dimensional data
    - How to cluster similar samples and genes?
objective:
    - To learn the principles of the analysis and visualisation of a multidimensional data analysis. We will use RNA-seq data as an example of a multidimensional -omics dataset. 
time_estimation: 4h
key_points:
    - A spliced mapping tool should be used on eukaryotic RNA-Seq data
    - Numerous factors should be taken into account when running a differential gene expression analysis
contributions:
  authorship:
    - pacthoen
    - casper937
    - charlotteradboud
---


In this computer assignment, we will analyze the results of an RNA-seq dataset. This is one type of high-dimensional -omics data that biomedical scientists frequently use. The raw RNA-seq data are sequence reads, but we will use the processed data from an RNA-seq experiment, often referred to as a count matrix. A count matrix contains for every sample (in the columns of the matrix) and every gene or transcript (in the rows of the matrix) the number of sequencing reads representing that gene or transcript (an integer, i.e., 0, 1, 2, 3, ...). The expression level of that gene is derived from the read count by correcting (normalization) for the total number of sequencing reads in a particular sample.

RNA-seq is often used to explore the function of a gene or the effect of a genetic variant. Since RNA-seq can survey the expression levels of all genes, it is well suited to obtain a global idea of the effect on the physiology of a cell and the biological processes affected. In this computer assignment, we will explore a dataset from mouse B-cells from TP53-/- knockout and wild-type mice. It is a two-factor experiment, because cells of both genotypes have been analyzed with or without being exposed to ionizing radiation (IR); thus, genotype and exposure to IR represent different factors in the experiment. The research paper describing this dataset is [Torelli et al., 2015, Oncotarget 6: 24611-24626](https://doi.org/10.18632/oncotarget.5232). The TP53 gene coding for the p53 protein is a well-known tumor suppressor gene. p53 is involved in the induction of apoptosis, senescence, cell cycle arrest and metabolic reprogramming. p53 is activated by double strand DNA breaks and routes cells to apoptosis when the DNA damage cannot be repaired anymore. TP53 mutation is one the most common driver of mutations in cancer, because TP53-deficient cells continue to proliferate, even after significant DNA damage. IR induces double strand DNA breaks and cell death (apoptosis) and is a well-known activator of p53. The authors used TP53-/- cells and wild-type cells to explore which part of the response of B-cells to IR is dependent on p53.

> <comment-title>Full data</comment-title>
>
> The original data are available at NCBI Gene Expression Omnibus (GEO) under accession number [GSE71176](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE71176). The raw RNA-Seq reads have been extracted from the Sequence Read Archive (SRA) files and converted into FASTQ files and counts.
>


# Data upload

In the first part of this tutorial we will upload the count matrices to your Galaxy environment.

> <hands-on-title>Data upload</hands-on-title>
>
> 1. Create a new history for this RNA-Seq exercise
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the count files from our [github repository](https://github.com/pacthoen/BMW2_RNA_clust_vis) using the Data Upload menu (top left) and the button _Paste/Fetch Data_. [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20115748.png). _Note_ The files with null in the filename are for the p53 knockout (p53-/-) cells and the files with p53 in the filename are for the wild-type cells
>    ```text
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_mock_1.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_mock_2.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_mock_3.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_mock_4.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_IR_1.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_IR_2.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_IR_3.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/p53_IR_4.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/null_mock_1.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/null_mock_2.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/null_IR_1.csv
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/null_IR_2.csv 
>    ``` 
>
> 3. Inspect some of the files by hitting the icon with the _eye_ button. 
>
>
> 4. Upload also the gene annotation file and choose _gtf.gz_ as file type
>    ```text
>    https://raw.githubusercontent.com/pacthoen/BMW2_RNA_clust_vis/refs/heads/main/data/Mus_musculus.NCBIM37.65.gtf.gz
>    ```

_Note:_ The numbers in the count files are the number of sequence reads to a specific genes. Only mRNAs that were detected in all samples were selected during the preprocessing of the data.

>    > <question-title></question-title>
>    >
>    > 1. How many mRNAs were detected in all samples?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1.By clicking on the dataset you can see the attributes of the dataset. There are 6,882 lines and a comment line. Thus, there are 6,882 detectable mRNAs in the dataset
>    > >
>    > {: .solution}
>    {: .question}


# Analysis of the differential gene expression

## Identification of the differentially expressed features

To be able to identify differential gene expression between wild-type and p53-null cells and between IR- and mock-treated cells, you might think we can just compare the count values in the files directly and calculate the extent of differential gene expression. 
 However, it is not that simple.

Let's imagine we have RNA-Seq counts from 3 samples for a genome with 4 genes:

Gene | Sample 1 Counts | Sample 2 Counts | Sample 3 Counts
--- | --- | --- | ---
A (2kb) | 10 | 12 | 30
B (4kb) | 20 | 25 | 60
C (1kb) | 5 | 8 | 15
D (10kb) | 0 | 0 | 1

Sample 3 has more reads than the other replicates, regardless of the gene. It has a higher sequencing depth than the other replicates. Gene B is twice as long as gene A: it might explain why it has twice as many reads, regardless of replicates.

The number of sequenced reads mapped to a gene therefore depends on:

- The **sequencing depth** of the samples

    Samples sequenced with more depth will have more reads mapping to each genes

- The **length of the gene**

    Longer genes will have more reads mapping to them

To compare samples or gene expressions, the gene counts need to be normalized. We could use TPM (Transcripts Per Kilobase Million).

> <details-title>RPKM, FPKM and TPM?</details-title>
>
> These three metrics are used to normalize count tables for:
>
> - sequencing depth (the "Million" part)
> - gene length (the "Kilobase" part)
>
> Let's use the previous example to explain RPKM, FPKM and TPM.
>
> For **RPKM** (Reads Per Kilobase Million),
>
> 1. Compute the "per million" scaling factor: sum up the total reads in a sample and divide that number by 1,000,000.
>
>    Gene | Sample 1 Counts | Sample 2 Counts | Sample 3 Counts
>    --- | --- | --- | ---
>    A (2kb) | 10 | 12 | 30
>    B (4kb) | 20 | 25 | 60
>    C (1kb) | 5 | 8 | 15
>    D (10kb) | 0 | 0 | 1
>    **Total reads** | 35 | 45 | 106
>    **Scaling factor** | 3.5 | 4.5 | 10.6
>
>    *Because of the small values in the example, we are using "per tens" instead of "per million" and therefore divide the sum by 10 instead of 1,000,000.*
>
> 2. Divide the read counts by the "per million" scaling factor
>
>    This normalizes for sequencing depth, giving reads per million (RPM)
>
>    Gene | Sample 1 RPM | Sample 2 RPM | Sample 3 RPM
>    --- | --- | --- | ---
>    A (2kb) | 2.86 | 2.67 | 2.83
>    B (4kb) | 5.71 | 5.56 | 5.66
>    C (1kb) | 1.43 | 1.78 | 1.43
>    D (10kb) | 0 | 0 | 0.09
>
>    *In the example we use the "per tens" scaling factor and we get reads per tens*
>
> 3. Divide the RPM values by the length of the gene, in kilobases.
>
>    Gene | Sample 1 RPKM | Sample 2 RPKM | Sample 3 RPKM
>    --- | --- | --- | ---
>    A (2kb) | 1.43 | 1.33 | 1.42
>    B (4kb) | 1.43 | 1.39 | 1.42
>    C (1kb) | 1.43 | 1.78 | 1.42
>    D (10kb) | 0 | 0 | 0.009
>
> **FPKM** (Fragments Per Kilobase Million) is very similar to RPKM. RPKM is used for single-end RNA-seq, while FPKM is used for paired-end RNA-seq. With single-end, every read corresponds to a single fragment that was sequenced. With paired-end RNA-seq, two reads of a pair are mapped from a single fragment, or if one read in the pair did not map, one read can correspond to a single fragment (in case we decided to keep these). FPKM keeps tracks of fragments so that one fragment with 2 reads is counted only once.
>
>
> **TPM** (Transcripts Per Kilobase Million) is very similar to RPKM and FPKM, except the order of the operation
>
> 1. Divide the read counts by the length of each gene in kilobases
>
>    This gives the reads per kilobase (RPK).
>
>    Gene | Sample 1 RPK | Sample 2 RPK | Sample 3 RPK
>    --- | --- | --- | ---
>    A (2kb) | 5 | 6 | 15
>    B (4kb) | 5 | 6.25 | 15
>    C (1kb) | 5 | 8 | 15
>    D (10kb) | 0 | 0 | 0.1
>
> 2. Compute the "per million" scaling factor: sum up all the RPK values in a sample and divide this number by 1,000,000
>
>    Gene | Sample 1 RPK | Sample 2 RPK | Sample 3 RPK
>    --- | --- | --- | ---
>    A (2kb) | 5 | 6 | 15
>    B (4kb) | 5 | 6.25 | 15
>    C (1kb) | 5 | 8 | 15
>    D (10kb) | 0 | 0 | 0.1
>    **Total RPK** | 15 | 20.25 | 45.1
>    **Scaling factor** | 1.5 | 2.03 | 4.51
>
>    *As above, because of the small values in the example, we use "per tens" instead of "per million" and therefore divide the sum by 10 instead of 1,000,000.*
>
> 3. Divide the RPK values by the "per million" scaling factor
>
>    Gene | Sample 1 TPM | Sample 2 TPM | Sample 3 TPM
>    --- | --- | --- | ---
>    A (2kb) | 3.33 | 2.96 | 3.33
>    B (4kb) | 3.33 | 3.09 | 3.33
>    C (1kb) | 3.33 | 3.95 | 3.33
>    D (10kb) | 0 | 0 | 0.1
>
>
> Unlike RPKM and FPKM, when calculating TPM, we normalize for gene length first, and then normalize for sequencing depth second. However, the effects of this difference are quite profound, as we already saw with the example.
>
> The sums of each column are very different:
>
> 1. RPKM
>
>    Gene | Sample 1 RPKM | Sample 2 RPKM | Sample 3 RPKM
>    --- | --- | --- | ---
>    A (2kb) | 1.43 | 1.33 | 1.42
>    B (4kb) | 1.43 | 1.39 | 1.42
>    C (1kb) | 1.43 | 1.78 | 1.42
>    D (10kb) | 0 | 0 | 0.009
>    **Total** | 4.29 | 4.5 | 4.25
>
> 2. TPM
>
>    Gene | Sample 1 TPM | Sample 2 TPM | Sample 3 TPM
>    --- | --- | --- | ---
>    A (2kb) | 3.33 | 2.96 | 3.33
>    B (4kb) | 3.33 | 3.09 | 3.33
>    C (1kb) | 3.33 | 3.95 | 3.33
>    D (10kb) | 0 | 0 | 0.1
>    **Total** | 10 | 10 | 10
>
> The sum of all TPMs in each sample are the same. This makes it easier to compare the proportion of reads that mapped to a gene in each sample. In contrast, with RPKM and FPKM, the sum of the normalized reads in each sample may be different, and this makes it harder to compare samples directly.
>
> In the example, TPM for gene A in Sample 1 is 3.33 and in sample 2 is 3.33. The same proportion of total reads maps then to gene A in both samples (0.33 here). Indeed, the sum of the TPMs in both samples adds up to the same number (10 here), the denominator required to calculate the proportions is then the same regardless of the sample, and so the proportion of reads for gene A (3.33/10 = 0.33) for both sample.
>
> With RPKM or FPKM, it is harder to compare the proportion of total reads because the sum of normalized reads in each sample can be different (4.29 for Sample 1 and 4.25 for Sample 2). Thus, if RPKM for gene A in Sample 1 is 1.43 and in Sample B is 1.43, we do not know if the same proportion of reads in Sample 1 mapped to gene A as in Sample 2.
>
> Since RNA-Seq is all about comparing relative proportion of reads, TPM seems more appropriate than RPKM/FPKM.
{: .details}

RNA-Seq is often used to compare one tissue type to another, for example, muscle vs. epithelial tissue. And it could be that there are a lot of muscle specific genes transcribed in muscle but not in the epithelial tissue. We call this a **difference in library composition**.

It is also possible to see a difference in library composition in the same tissue type after the knock out of a transcription factor.

Let's imagine we have RNA-Seq counts from 2 samples (same library size: 635 reads), for a genome with 6 genes. The genes have the same expression in both samples, except one: only Sample 1 transcribes gene D, at a high level (563 reads). As the library size is the same for both samples, sample 2 has 563 extra reads to be distributed over genes A, B, C, E and F.

Gene | Sample 1 | Sample 2
--- | --- | --- | ---
A | 30 | 235
B | 24 | 188
C | 0 | 0
D | 563 | 0
E | 5 | 39
F | 13 | 102
**Total** | 635 | 635

As a result, the read count for all genes except for genes C and D is really high in Sample 2. Nonetheless, the only differentially expressed gene is gene D.

TPM, RPKM or FPKM do not deal with these differences in library composition during normalization, but more complex tools, like DESeq2, do.

[**DESeq2**](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) ({% cite love2014moderated %}) is a great tool for dealing with RNA-seq data and running Differential Gene Expression (DGE) analysis. It takes read count files from different samples, combines them into a big table (with genes in the rows and samples in the columns) and applies normalization for **sequencing depth** and **library composition**. We do not need to account for gene length normalization does because we are comparing the counts between sample groups for the same gene.

> <details-title>Normalization in DESeq2</details-title>
>
> Let's take an example to illustrate how DESeq2 scales the different samples:
>
> Gene | Sample 1 | Sample 2 | Sample 3
> A | 0 | 10 | 4
> B | 2 | 6 | 12
> C | 33 | 55 | 200
>
> The goal is to calculate a scaling factor for each sample, which takes read depth and library composition into account.
>
> 1. Take the log$$_e$$ of all the values:
>
>     Gene | log(Sample 1) | log(Sample 2) | log(Sample 3)
>     A | -Inf | 2.3 | 1.4
>     B | 0.7 | 1.8 | 2.5
>     C | 3.5 | 4.0 | 5.3
>
> 2. Average each row:
>
>     Gene | Average of log values
>     A | -Inf
>     B | 1.7
>     C | 4.3
>
>     The average of the log values (also known as the geometric average) is used here because it is not easily impacted by outliers (e.g. gene C with its outlier for Sample 3).
>
> 3. Filter out genes with which have infinity as a value.
>
>     Gene | Average of log values
>      |
>     B | 1.7
>     C | 4.3
>
>     Here we filter out genes with no read counts in at least 1 sample, e.g. genes only transcribed in one tissue like gene D in the previous example. This helps to focus the scaling factors on genes transcribed at similar levels, regardless of the condition.
>
> 4. Subtract the average log value from the log counts:
>
>     Gene | log(Sample 1) | log(Sample 2) | log(Sample 3)
>      |  |  |
>     B | -1.0 | 0.1 | 0.8
>     C | -0.8 | -0.3 | 1.0
>
>     $$log(\textrm{counts for gene X}) - average(\textrm{log values for counts for gene X}) = log(\frac{\textrm{counts for gene X}}{\textrm{average for gene X}})$$
>
>     This step compares the ratio of the counts in each sample to the average across all samples.
>
> 5. Calculate the median of the ratios for each sample:
>
>     Gene | log(Sample 1) | log(Sample 2) | log(Sample 3)
>      |  |  |
>     B | -1.0 | 0.1 | 0.8
>     C | -0.8 | -0.3 | 1.0
>     **Median** | -0.9 | -0.1 | 0.9
>
>     The median is used here to avoid extreme genes (most likely rare ones) from swaying the value too much in one direction. It helps to put more emphasis on moderately expressed genes.
>
> 6. Compute the scaling factor by taking the exponential of the medians:
>
>     Gene | Sample 1 | Sample 2 | Sample 3
>     **Median** | -0.9 | -0.1 | 0.9
>     **Scaling factors** | 0.4 | 0.9 | 2.5
>
> 7. Compute the normalized counts: divide the original counts by the scaling factors:
>
>     Gene | Sample 1 | Sample 2 | Sample 3
>     A | 0 | 11.11 | 1.6
>     B | 5 | 6.67 | 4.8
>     C | 83 | 61.11 | 80
>
> *This explanation is a transcription and adaptation of the [StatQuest video explaining Library Normalization in DESEq2](https://www.youtube.com/watch?v=UFB993xufUU&t=35s)*.
>
{: .details}

DESeq2 also runs the Differential Gene Expression (DGE) analysis, which has two basic tasks:

- Estimate the biological variance using the replicates for each condition
- Estimate the significance of expression differences between any two conditions

This expression analysis is estimated from read counts and attempts are made to correct for variability in measurements using replicates, that are absolutely essential for accurate results. 

> <details-title>Technical vs biological replicates</details-title>
>
> A technical replicate is an experiment which is performed once but measured several times (e.g. multiple sequencing of the same library). A biological replicate is an experiment performed (and also measured) several times.
>
> In our data, we have 2-4 biological replicates (here called samples) per condition.
>    > <question-title></question-title>
>    >
>    > 1. How many conditions are there in this experiment and what is the number of replicates in each of the conditions?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1.There are 4 conditions
>    > > - Wild-type cells not exposed to IR - 4 replicates
>    > > - Wild-type cells exposed to IR - 4 replicates
>    > > - TP53-/- cells not exposed to IR - 2 replicates
>    > > - TP53-/- cells exposed to IR - 2 replicates
>    > >
>    > {: .solution}
>    {: .question}

> We recommend to combine the count tables for different technical replicates (but not for biological replicates) before a differential expression analysis (see [DESeq2 documentation](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#collapsing-technical-replicates))
{: .details}

Multiple factors with several levels can then be incorporated in the analysis describing known sources of variation (e.g. treatment, tissue type, gender, batches), with two or more levels representing the conditions for each factor. After normalization we can compare the response of the expression of any gene to the presence of different levels of a factor in a statistically reliable way.

In our example, we have samples with two varying factors that can contribute to differences in gene expression:

- Genotype (either wild-type or knock-out)
- Treatment (either IR or mock)

The effect of Ionizing Radiation (IR) is our primary interest. Later on, we will explore the differences in IR response between wild-type and knock-out cells and which genes are dependent on p53 for their response to IR.
We will now assign each sample to the right group

{% include _includes/cyoa-choices.html option1="Basic" option2="Tag-based" option3="Collection split" default="Basic" text="Which approach would you prefer to use?" disambiguation="deseq"%}

<div class="Basic" markdown="1">

We can now run **DESeq2**:

> <hands-on-title>Determine differentially expressed features</hands-on-title>
>
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} with the following parameters:
>    - *"how"*: `Select datasets per level`
>        - In *"Factor"*:
>           - *"Specify the factor name"*: `Treatment`
>           - In *"1: Factor level"*:
>               - *"Specify the factor level"*: `IR`
>               - In *"Count file(s)"*: `Select all the IR treated count files`. _Note1:_ Use the _Switch to column select_ option to select the files. _Note2:_ The first factor level is compared to the second factor level. In this case IR vs. mock. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20123708.png)
>           - In *"2: Factor level"*:
>               - *"Specify the factor level"*: `Mock`
>               - In *"Count file(s)"*: `Select all the mock treated count files`
>       - {% icon param-repeat %} *"Insert Factor"*
>           - *"Specify the factor name"*: `Genotype`
>               - In *"Factor level"*:
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify the factor level"*: `KO`
>                        - In *"Count file(s)"*: `Select all the KO count files.`  i.e. those with null in the file name
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify the factor level"*: `WT`
>                        - In *"Count file(s)"*: `Select all the WT count files` i.e. those with p53 in the file name
>    - *"Files have header?"*: `Yes`
>    - *"Choice of Input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - In *"Advanced options"*:
>        - *"Use beta priors"*: `Yes`
>    - In *"Output options"*:
>        - *"Output selector"*: `Generate plots for visualizing the analysis results`, `Output normalised counts`, `Output VST normalized table`
>
{: .hands_on}


**DESeq2** generated 4 outputs:

- A table with the normalized counts for each gene (rows) in the samples (columns)
- A graphical summary of the results, useful to evaluate the quality of the experiment:

    1. A plot of the first 2 dimensions from a principal component analysis 

    2. Heatmap of the sample-to-sample distance matrix (with clustering) based on the normalized counts.

        The heatmap gives an overview of similarities and dissimilarities between samples: the color represents the distance between the samples. Dark blue means shorter distance, i.e. closer samples given the normalized counts.


    3. Dispersion estimates: gene-wise estimates (black), the fitted values (red), and the final maximum a posteriori estimates used in testing (blue)

        This dispersion plot is typical, with the final estimates shrunk from the gene-wise estimates towards the fitted estimates. Some gene-wise estimates are flagged as outliers and not shrunk towards the fitted value. The amount of shrinkage can be more or less than seen here, depending on the sample size, the number of coefficients, the row mean and the variability of the gene-wise estimates.

    4. Histogram of *p*-values for the genes in the comparison between the 2 levels of the 1st factor

    5. An [MA plot](https://en.wikipedia.org/wiki/MA_plot):

        This displays the global view of the relationship between the expression change of conditions (log ratios, M), the average expression strength of the genes (average mean, A), and the ability of the algorithm to detect differential gene expression. The genes that passed the significance threshold (adjusted p-value < 0.1) are colored in red.

- A summary file with the following values for each gene:

    1. Gene identifiers
    2. Mean normalized counts, averaged over all samples from both conditions
    3. Fold change in log2 (logarithm base 2)

       The log2 fold changes are based on the primary factor level 1 vs factor level 2, hence the input order of factor levels is important. Here, DESeq2 computes fold changes of 'treated' samples against 'untreated' from the first factor 'Treatment', *i.e.* the values correspond to up- or downregulation of genes in treated samples.

    4. Standard error estimate for the log2 fold change estimate
    5. [Wald](https://en.wikipedia.org/wiki/Wald_test) statistic
    6. *p*-value for the statistical significance of this change
    7. *p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure, which controls false discovery rate ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

    > <tip-title>What are p-values and what are they used for?</tip-title>
    >
    > The p-value is a measure often used to determine whether or not a particular observation possesses statistical significance. Strictly speaking, the p-value is the probability that the data could have arisen randomly, assuming that the null hypothesis is correct. In the concrete case of RNA-Seq, the null hypothesis is that there is no differential gene expression. So a p-value of 0.13 for a particular gene indicates that, for that gene, assuming it is not differentially expressed, there is a 13% chance that any apparent differential expression could simply be produced by random variation in the experimental data.
    >
    > 13% is still quite high, so we cannot really be confident differential gene expression is taking place. The most common way that scientists use p-values is to set a threshold (commonly 0.05, sometimes other values such as 0.01) and reject the null hypothesis only for p-values below this value. Thus, for genes with p-values less than 0.05, we can feel safe stating that differential gene expression plays a role. It should be noted that any such threshold is arbitrary and there is no meaningful difference between a p-value of 0.049 and 0.051, even if we only reject the null hypothesis in the first case.
    >
    > Unfortunately, p-values are often heavily misused in scientific research, enough so that Wikipedia provides a [dedicated article](https://en.wikipedia.org/wiki/Misuse_of_p-values) on the subject. See also [this article](https://fivethirtyeight.com/features/not-even-scientists-can-easily-explain-p-values/) (aimed at a general, non-scientific audience).
    {: .tip}

For more information about **DESeq2** and its outputs, you can have a look at the [**DESeq2** documentation](https://www.bioconductor.org/packages/release/bioc/manuals/DESeq2/man/DESeq2.pdf).

- An output file with Normalized counts
- An output file with VST-Normalized counts (VST stands for Variance Stabilization and Transformation; this is the preferred normalization method in DESeq2 and we will use these counts later on for PCA and clustering

Use the **filter** tool on the summary file to identify the genes with absolute fold-change greater than 2. Use the condition: abs(c3)>2. Number of header lines to skip: 1. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20130228.png). And use the **sort** tool to identify the genes with the largest fold-change (most upregulated). See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20132050.png) 
> <question-title></question-title>
>
> 1. What does it mean when an mRNA has an absolute 2logFC greater than 2?
> 2. How many mRNAs have an absolute 2logFC greater than 2?
> 3. What is the mRNA that is most upregulated in and what is its false discovery rate (the chance that this is a false positive)? 
>
> > <solution-title></solution-title>
> >
> > 1. That the genes is 2^2 = 4 times higher or lower expressed in the cells treated with IR than in mock-treated cells
> > 2. 666 (!) 
> > 3. ENSMUSG00000000308; FDR is in p-adj column. It is 1.3*10^-26
> >    
> {: .solution}
{: .question}

## Annotation of the DESeq2 results

The ID for each gene is something like ENSMUSG00000026581, which is an ID from the Ensembl gene database. These IDs are unique but sometimes we prefer to have the gene names, even if they may not reference an unique gene (e.g. duplicated after re-annotation). But gene names may hint already to a function or they help you to search for desired candidates. We would also like to display the location of these genes within the genome. We can extract such information from the annotation file which we uploaded before.

> <hands-on-title>Annotation of the DESeq2 results</hands-on-title>
>
>
> 1. {% tool [Annotate DESeq2/DEXSeq output tables](toolshed.g2.bx.psu.edu/repos/iuc/deg_annotate/deg_annotate/1.1.0) %} with:
>    - {% icon param-file %} *"Tabular output of DESeq2/edgeR/limma/DEXSeq"*: the `DESeq2 result file` (output of **DESeq2** {% icon tool %})
>    - *"Input file type"*: `DESeq2/edgeR/limma`
>    - {% icon param-file %} *"Reference annotation in GFF/GTF format"*: imported gtf `Annotation file`
>    - See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20132922.png)
{: .hands_on}

The generated output is an extension of the previous file:

1. Gene identifiers
2. Mean normalized counts over all samples
3. Log2 fold change
4. Standard error estimate for the log2 fold change estimate
5. Wald statistic
6. *p*-value for the Wald statistic
7. *p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure for the Wald statistic
8. Chromosome
9. Start
10. End
11. Strand
12. Feature
13. Gene name

> <question-title></question-title>
>
> 1. What is the Gene Symbol of the most overexpressed gene
> 2. On which chromosome is the most over-expressed gene located?
>  
> > <solution-title></solution-title>
> >
> > 1. ENSMUSG00000000308 (the top-ranked gene with the highest positive log2FC value)
> > 2. It is located on the reverse strand of chromosome X, between 10,778,953 bp and 10,786,907 bp.
> {: .solution}
{: .question}

The annotated table contains no column names, which makes it difficult to read. We would like to add them before going further.

> <hands-on-title>Add column names</hands-on-title>
>
> 1. Create a new file (`header`) from the following (header line of the DESeq2 output)
>
>    ```text
>    GeneID	Base mean	log2(FC)	StdErr	Wald-Stats	P-value	P-adj	Chromosome	Start	End	Strand	Feature	Gene name
>    ```
>
>    {% snippet faqs/galaxy/datasets_create_new_file.md name="header" format="tabular" %}
>    See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20134800.png)
>
> 3. {% tool [Concatenate datasets](cat1) %} to add this header line to the **Annotate** output:
>    - {% icon param-file %} *"Concatenate Dataset"*: the `header` dataset
>    - *"Dataset"*
>       - Click on {% icon param-repeat %} *"Insert Dataset"*
>         - {% icon param-file %} *"select"*: output of **Annotate** {% icon tool %}
>
> 4. Rename the output to `Annotated DESeq2 results`
{: .hands_on}

## Create Volcano plot to visualise differentially expressed genes

Select the tool _Volcano Plot_
> In a Volcano plot, the 2logFC (x-axis) is plotted against the -log10(P-value).
> > Choose the right column numbers. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20162602.png). Change the significant threshold to 0 and LogFC 2 as as thresholds to colour.
> >  Put in correct titles under _Plot options_.
> >
Repeat the DESeq2 analysis but now with `Genotype` as first factor (contrasting KO and WT cells). Repeat the Volcano plot. 

> <question-title></question-title>
>
> 1. Is IR causing more up- or more downregulation of gene expression?
> 2. Does the Treatment or the Genotype have a bigger effect on gene expression
>  
> > <solution-title></solution-title>
> >
> > 1. There are more mRNA upregulated than downregulated. The upregulated genes can be found on the right side of the plot. 
> > 2. There are more mRNA with high FC and strong significance in the comparison between IR and mock-treated cells than between KO and WT cells.
> {: .solution}
{: .question}



# Functional enrichment analysis of the DE genes

We have extracted genes that are differentially expressed in IR- vs. mock-treated samples. Now, we would like to know if the differentially expressed genes are enriched transcripts of genes which belong to more common or specific categories in order to identify biological functions that might be impacted.

## Gene Ontology analysis

[Gene Ontology (GO)](http://www.geneontology.org/) analysis is widely used to reduce complexity and highlight biological processes in genome-wide expression studies. However, standard methods give biased results on RNA-Seq data due to over-detection of differential expression for long and highly-expressed transcripts.

[**goseq**](https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf) ({% cite young2010gene %}) provides methods for performing GO analysis of RNA-Seq data while taking length bias into account. **goseq** could also be applied to other category-based tests of RNA-Seq data, such as KEGG pathway analysis, as discussed in a further section.

**goseq** needs 2 files as inputs:

- A tabular file with the differentially expressed genes from all genes assayed in the RNA-Seq experiment with 2 columns:
  - the Gene IDs (unique within the file), in uppercase letters
  - a boolean indicating whether the gene is differentially expressed or not (`True` if differentially expressed or `False` if not). We will use abs(2LogFC) > 2 as a criterion
- A file with information about the length of a gene to correct for potential length bias in differentially expressed genes


> 1. Use the _Gene length and GC content_ tool on the _Annotation file_ (gtf format). See this [screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20210701.png)
> 2. Merge the gene length and the annotated DESeq2 ouput file **DESeq2 result file** with primary factor Treatment using the _Join two datasets side by side_ tool. In *"Keep the header lines"*: `No`. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20214135.png)
> 3. Use the _Compute on rows_ tool to create a column with TRUE and FALSE using the following expression: `bool(float(c3)>2)`. In the *"Error handling"* choose in *"Autodetect column types"* `No` and *"Fail on references to non-existent columns"* `No` and *"If an expression cannot be computed for a row"* choose `Fill in a replacement value` and Replacement value `False`. See this [screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20222441.png)
> 4. {% tool [Cut](Cut1) %} columns from a table with the following parameters:
>    - *"Cut columns"*: `c1,c9`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the output generated in step 3
> 5. Rename the output to `Gene IDs and differential expression - IR`
> 6. {% tool [Cut](Cut1) %} columns from a table with the following parameters:
>    - *"Cut columns"*: `c1,c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the output generated in step 3
> 7. Rename the output to `Gene length differential expression - IR`
{: .hands_on}

We have now the two required input files for goseq.

> <hands-on-title>Perform GO analysis</hands-on-title>
>
> 1. {% tool [goseq](toolshed.g2.bx.psu.edu/repos/iuc/goseq/goseq/1.50.0+galaxy0) %} with
>    - *"Differentially expressed genes file"*: `Gene IDs and differential expression`
>    - *"Gene lengths file"*: `Gene IDs and length`
>    - *"Gene categories"*: `Get categories`
>       - *"Select a genome to use"*: `mouse (mm10)`
>       - *"Select Gene ID format"*: `Ensembl Gene ID`
>       - *"Select one or more categories"*: `GO: Cellular Component`, `GO: Biological Process`, `GO: Molecular Function`
>    - In *"Output Options"*
>      - *"Output Top GO terms plot?"*: `Yes`
>      - *"Extract the DE genes for the categories (GO/KEGG terms)?"*: `Yes`
>
{: .hands_on}

**goseq** generates with these parameters 3 outputs:

1. A table (`Ranked category list - Wallenius method`) with the following columns for each GO term:

    1. `category`: GO category
    2. `over_rep_pval`: *p*-value for over-representation of the term in the differentially expressed genes
    3. `under_rep_pval`: *p*-value for under-representation of the term in the differentially expressed genes
    4. `numDEInCat`: number of differentially expressed genes in this category
    5. `numInCat`: number of genes in this category
    6. `term`: detail of the term
    7. `ontology`: MF (Molecular Function - molecular activities of gene products), CC (Cellular Component - where gene products are active), BP (Biological Process - pathways and larger processes made up of the activities of multiple gene products)
    8. `p.adjust.over_represented`: *p*-value for over-representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure
    9. `p.adjust.under_represented`: *p*-value for under-representation of the term in the differentially expressed genes, adjusted for multiple testing with the Benjamini-Hochberg procedure

    To identify categories significantly enriched/unenriched below some p-value cutoff, it is necessary to use the adjusted *p*-value.

    > <question-title></question-title>
    >
    > 1. How many GO terms are over-represented with an adjusted P-value < 0.05? How many are under-represented?
    > 2. Do you find evidence that p53 plays a role in the response to IR?
    >
    > > <solution-title></solution-title>
    > >
    > > 1. 5 GO terms (4 from the Biological Processes category, 1 from the Molecular Function category) are over-represented.
    > > 2. Signal transduction by p53 class mediator is second in rank and contains 19 differentially expressed genes.
    > >    {% tool [Filter data on any column using simple expressions](Filter1) %} on c8 (adjusted p-value for over-represented GO terms) and c9 (adjusted p-value for under-represented GO terms)
    > >
    > > 2. For over-represented, 50 BP, 5 CC and 5 MF and for under-represented, 5 BP, 2 CC and 0 MF
    > >
    > >    {% tool [Group data](Grouping1) %} on column 7 (category) and count on column 1 (IDs)
    > >
    > {: .solution}
    {: .question}

2. A graph with the top 10 over-represented GO terms

    > <question-title></question-title>
    >
    > ![Top over-represented GO terms](../../images/ref-based/top_over-represented_go_terms.png)
    >
    > What is the x-axis? How is it computed?
    >
    > > <solution-title></solution-title>
    > >
    > > The x-axis is the percentage of genes in the category that have been identified as differentially expressed: $$100 \times \frac{numDEInCat}{numInCat}$$
    > >
    > {: .solution}
    {: .question}

3. A table with the differentially expressed genes (from the list we provided) associated to the GO terms (`DE genes for categories (GO/KEGG terms)`)

> <comment-title>Advanced tutorial on enrichment analysis</comment-title>
>
> In this tutorial, we covered GO enrichment analysis with **goseq**. You can also use the same method to identify metabolic and signal transduction pathways from the KEGG database (do not select plotting the Top GO pathways in that case). To learn other methods and tools for gene set enrichment analysis, please have a look at the ["RNA-Seq genes to pathways"]({% link topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.md %}) tutorial.
{: .comment}


# Visualisation of high-dimensional datasets
In the second part of this computer assignment we will use commonly used techniques to identify the main sources of variation and to evaluate which samples and sample groups are most similar. We will continue using the same dataset but first filter for the most variable genes. This is for computational efficiency and to reduce the effect of noise.


## Filter for most variable genes
For filtering the most variable genes, we will compute the coefficient of variation (CV). This is defined as the standard deviation divided by the mean. 
> 1. Use the the _Table compute_ tool to calculate the standard deviation (std)
>    - *"Table"*: `VST normalized data`
>    - *"Input data has"*: `Column names on the first row` and `Row names on the first column`
>    - *"Type of table operation "*: `Compute expression across rows or columns`
>       - *"Calculate"*: `Custom`
>       - *"Custom function on vec"*: `vec.std`
>       - *"For each"*: `row`
>    - See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-23%20085500.png)
> 2. Use the the _Table compute_ tool to calculate the mean (mean)
>    - *"Table"*: `VST normalized data`
>    - *"Input data has"*: `Column names on the first row` and `Row names on the first column`
>    - *"Type of table operation "*: `Compute expression across rows or columns`
>       - *"Calculate"*: `Custom`
>       - *"Custom function on vec"*: `vec.mean`
>       - *"For each"*: `row`
> 3. Use the _Join two files_ tool to merge the normalised count data with the standard deviation and mean datasets (two steps)
>    -  *"First line is a header line"*: `Yes`
> 4. Use the _Compute on row_ tool to calculate the CV. Expression `c14/c15`. Label the new row as `CV`
> 5. Use the _Sort_ tool to sort the genes based on CV in descending order. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-23%20104444.png)
>    - *"Number of header lines"*: `1`
>    - *"Column selections"*
>       - *"on column"*: `16`
>       - *"in"*: `decsending order`
>       - *"Flavour"*: `Fast numeric sort (-n)`     
> 6. Select the first 2000 rows using the _Select first_ tool. See this [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-23%20105700.png)
> 7. Select only the columns with data (and the first column with the gene identifiers; columns 1-13) using the _Cut_ tool
> 8. Relabel the file as `Input data for PCA and clustering` 
>    
>    - 


## Principal Component Analysis (PCA)
PCA is a dimension reduction technique to identify the main sources of variation. 
For this part of the assignment we will use the _multivariate_ tool. There is a separate [tutorial](https://training.galaxyproject.org/topics/metabolomics/tutorials/lcms-dataprocessing/tutorial.html#bibliography) for the _multivariate_ tool applied to metabolomics data.   

The tool needs a sampleMetadata file where the different groups / factors are defined

In the sampleMetadata file provided in this tutorial, you can find a column named *Osmo* that is meant to correspond to measures of osmolality,
one way to represent the overall concentration of urine samples.
Let's try to colour the PCA score plots according to that variable. For this, we will use the **Multivariate** {% icon tool %} Galaxy module.

> <hands-on-title>Using <b>Multivariate</b> to visualise the two first components of a PCA</hands-on-title>
>
> 1. {% tool [Multivariate](toolshed.g2.bx.psu.edu/repos/ethevenot/multivariate/Multivariate/2.3.10) %} with the following parameters:
>    - {% icon param-file %} *"Data matrix file"*: `Filtered_dataMatrix` (output of the last **Generic_Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Sample metadata file"*: `Filtered_sampleMetadata` (output of the last **Generic_Filter** {% icon tool %} job)
>    - {% icon param-file %} *"Variable metadata file"*: `Filtered_variableMetadata` (output of the last **Generic_Filter** {% icon tool %} job)
>    - *"Number of predictive components"*: `2`
>    - *"Advanced graphical parameters"*: `Full parameter list`
>        - *"Sample colors"*: `Osmo`
>        - *"Amount by which plotting text should be magnified relative to the default"*: `0.4`
>    - *"Advanced computational parameters"*: `Use default`
>
>
>    > <tip-title>Comment to W4M users</tip-title>
>    >
>    > In the [GTN_LCMSprocessing](https://workflow4metabolomics.usegalaxy.fr/u/peteram/h/gtnlcmsdataprocessing) history, this step corresponds to the datasets number 33 to 36.
>    {: .tip}
>
{: .hands_on}

Among the four output files generated by the tool, the one we are interested in is the *Multivariate_figure.pdf* file.
The PCA scores plot representing the projection of samples on the two first components is located at the bottom left corner:

![The picture represents the PCA scores plot of components t1 (32% of total variability) and t2 (9% of total variability) with sample labels coloured according to the values from the Osmo column of the sampleMetadata file. We observe a gradation in the colours of samples following the first bisector of the graphic.](../../images/lcmsproc_osmo.png "PCA scores plot")

Taking into account that the PCA was computed on Unit-Variance-scaled data (see the **Multivariate** {% icon tool %} tool's help section for more details),
this result confirms the impact of total concentration on our dataset.
Since this effect in data is independant of the question of interest that we suppose we investigate in the study (in the case of this tutorial),
we may want to get rid of this effect to prevent a power reduction in further statistical analysis.



## Clustering


# Conclusion

In this computer assignment, we have analyzed real RNA sequencing data to extract useful information, such as which genes are up or downregulated by IR. Your complete workflow can be extracted using the _Worklows_ button. Do this. You can always rerun and adapt your entire workflow.

