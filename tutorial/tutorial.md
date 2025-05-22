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
    - cdevisser
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
> 2. Import the count files from our [github repository](https://github.com/pacthoen/BMW2_RNA_clust_vis) using the Data Upload menu (top left) and the button _Paste/Fetch Data_ and choosing _Tabular_ as file type. [Screenshot](https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/screenshots/Screenshot%202025-05-22%20105100.png)
>    ```text
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_mock_1.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_mock_2.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_mock_3.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_mock_4.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_IR_1.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_IR_2.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_IR_3.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/p53_IR_4.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/null_mock_1.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/null_mock_2.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/null_IR_1.csv
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/null_IR_2.csv
>    ``` 
>
> 3. Inspect some of the files by hitting the icon with the _eye_ button. 
>
>
> 4. Upload also the gene annotation file and choose _gtf.gz_ as file type
>    ```text
>    https://github.com/pacthoen/BMW2_RNA_clust_vis/blob/main/data/Mus_musculus.NCBIM37.65.gtf.gz
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

To be able to identify differential gene expression induced by PS depletion, all datasets (3 treated and 4 untreated) must be analyzed following the same procedure. To save time, we have run the previous steps for you. We then obtain 7 files with the counts for each gene of *Drosophila* for each sample.

> <hands-on-title>Import all count files</hands-on-title>
>
> 1. Create a **new empty history**
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the seven count files from [Zenodo]({{ page.zenodo_link }}) or the Shared Data library:
>
>    - `GSM461176_untreat_single_featureCounts.counts`
>    - `GSM461177_untreat_paired_featureCounts.counts`
>    - `GSM461178_untreat_paired_featureCounts.counts`
>    - `GSM461179_treat_single_featureCounts.counts`
>    - `GSM461180_treat_paired_featureCounts.counts`
>    - `GSM461181_treat_paired_featureCounts.counts`
>    - `GSM461182_untreat_single_featureCounts.counts`
>
>    ```text
>    {{ page.zenodo_link }}/files/GSM461176_untreat_single_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461177_untreat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461178_untreat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461179_treat_single_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461180_treat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461181_treat_paired_featureCounts.counts
>    {{ page.zenodo_link }}/files/GSM461182_untreat_single_featureCounts.counts
>    ```
>
{: .hands_on}

You might think we can just compare the count values in the files directly and calculate the extent of differential gene expression. However, it is not that simple.

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

This expression analysis is estimated from read counts and attempts are made to correct for variability in measurements using replicates, that are absolutely essential for accurate results. For your own analysis, we advise you to use at least 3, but preferably 5 biological replicates per condition. It is possible to have different numbers of replicates per condition.

> <details-title>Technical vs biological replicates</details-title>
>
> A technical replicate is an experiment which is performed once but measured several times (e.g. multiple sequencing of the same library). A biological replicate is an experiment performed (and also measured) several times.
>
> In our data, we have 4 biological replicates (here called samples) without treatment and 3 biological replicates with treatment (*Pasilla* gene depleted by RNAi).
>
> We recommend to combine the count tables for different technical replicates (but not for biological replicates) before a differential expression analysis (see [DESeq2 documentation](http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#collapsing-technical-replicates))
{: .details}

Multiple factors with several levels can then be incorporated in the analysis describing known sources of variation (e.g. treatment, tissue type, gender, batches), with two or more levels representing the conditions for each factor. After normalization we can compare the response of the expression of any gene to the presence of different levels of a factor in a statistically reliable way.

In our example, we have samples with two varying factors that can contribute to differences in gene expression:

- Treatment (either treated or untreated)
- Sequencing type (paired-end or single-end)

Here, treatment is the primary factor that we are interested in. The sequencing type is further information we know about the data that might affect the analysis. Multi-factor analysis allows us to assess the effect of the treatment, while taking the sequencing type into account too.

> <comment-title></comment-title>
>
> We recommend that you add all factors you think may affect gene expression in your experiment. It can be the sequencing type like here, but it can also be the manipulation (if different persons are involved in the library preparation), other batch effects, etc...
{: .comment}

If you have only one or two factors with few number of biological replicates, the basic setup of **DESeq2** is enough. In the case of a complex experimental setup with a large number of biological replicates, tag-based collections are appropriate. Both approaches give the same results. The Tag-based approach requires a few additional steps before running the **DESeq2** tool but it will payoff when working with a complex experimental setup.

{% include _includes/cyoa-choices.html option1="Basic" option2="Tag-based" option3="Collection split" default="Basic" text="Which approach would you prefer to use?" disambiguation="deseq"%}

<div class="Basic" markdown="1">

We can now run **DESeq2**:

> <hands-on-title>Determine differentially expressed features</hands-on-title>
>
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} with the following parameters:
>    - *"how"*: `Select datasets per level`
>        - In *"Factor"*:
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Treatment`
>           - In *"1: Factor level"*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `treated`
>               - In *"Count file(s)"*: `Select all the treated count files (GSM461179, GSM461180, GSM461181)`
>           - In *"2: Factor level"*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `untreated`
>               - In *"Count file(s)"*: `Select all the untreated count files (GSM461176, GSM461177, GSM461178, GSM461182)`
>       - {% icon param-repeat %} *"Insert Factor"*
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Sequencing`
>               - In *"Factor level"*:
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `PE`
>                        - In *"Count file(s)"*: `Select all the paired-end count files (GSM461177, GSM461178, GSM461180, GSM461181)`
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `SE`
>                        - In *"Count file(s)"*: `Select all the single-end count files (GSM461176, GSM461179, GSM461182)`
>    - *"Files have header?"*: `Yes`
>    - *"Choice of Input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - In *"Advanced options"*:
>        - *"Use beta priors"*: `Yes`
>    - In *"Output options"*:
>        - *"Output selector"*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
>
{: .hands_on}

</div>

<div class="Tag-based" markdown="1">

DESeq2 requires to provide for each factor, counts of samples in each category. We will thus use tags on our collection of counts to easily select all samples belonging to the same category. For more information about alternative ways to set group tags, please see [this tutorial]({% link topics/galaxy-interface/tutorials/group-tags/tutorial.md %}).

> <hands-on-title>Add tags to your collection for each of these factors</hands-on-title>
>
> 1. Create a collection list with all these counts that you label `all counts`. Name each item so it only has the GSM id, the treatment and the library, for example, `GSM461176_untreat_single`.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> 2. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} with the following parameters:
>    - {% icon param-collection %} *"Dataset collection"*: `all counts`
>
>    We will now extract from the names the factors:
>
> 3. {% tool [Replace Text in entire line](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_replace_in_line/9.3+galaxy1) %}
>      - {% icon param-file %} *"File to process"*: output of **Extract element identifiers** {% icon tool %}
>      - In *"Replacement"*:
>         - In *"1: Replacement"*
>            - *"Find pattern"*: `(.*)_(.*)_(.*)`
>            - *"Replace with"*: `\1_\2_\3\tgroup:\2\tgroup:\3`
>
>     This step creates 2 additional columns with the type of treatment and sequencing that can be used with the {% tool [Tag elements](__TAG_FROM_FILE__) %} tool
>
> 4. Change the datatype to `tabular`
>
>    {% snippet faqs/galaxy/datasets_change_datatype.md datatype="tabular" %}
>
> 5. {% tool [Tag elements](__TAG_FROM_FILE__) %}
>      - {% icon param-collection %} *"Input Collection"*: `all counts`
>      - {% icon param-file %} *"Tag collection elements according to this file"*: output of **Replace Text** {% icon tool %}
>
> 6. Inspect the new collection
>
>    > <tip-title>You cannot see the changes?</tip-title>
>    >
>    > You may not see it at first glance as the names are the same. However if you click on one and click on {% icon galaxy-tags %} **Edit dataset tags**, you should see 2 tags which start with 'group:'. This keyword will allow to use these tags in **DESeq2**.
>     >
>     {: .tip}
>
{: .hands_on}

We can now run **DESeq2**:

> <hands-on-title>Determine differentially expressed features</hands-on-title>
>
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} with the following parameters:
>    - *"how"*: `Select group tags corresponding to levels`
>        - {% icon param-collection %} *"Count file(s) collection"*: output of **Tag elements** {% icon tool %}
>        - In *"Factor"*:
>            - {% icon param-repeat %} *"Insert Factor"*
>                - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Treatment`
>                - In *"Factor level"*:
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `treated`
>                        - *"Select groups that correspond to this factor level"*: `Tags: treat`
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `untreated`
>                        - *"Select groups that correspond to this factor level"*: `Tags: untreat`
>            - {% icon param-repeat %} *"Insert Factor"*
>                - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Sequencing`
>                - In *"Factor level"*:
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `PE`
>                        - *"Select groups that correspond to this factor level"*: `Tags: paired`
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `SE`
>                        - *"Select groups that correspond to this factor level"*: `Tags: single`
>    - *"Files have header?"*: `Yes`
>    - *"Choice of Input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - In *"Advanced options"*:
>        - *"Use beta priors"*: `Yes`
>    - In *"Output options"*:
>        - *"Output selector"*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
>
{: .hands_on}

</div>
<div class="Collection-split" markdown="1">

DESeq2 requires to provide for each factor, counts of samples in each category. We will thus use patterns on the name of our samples to easily select all samples belonging to the same category.

> <hands-on-title>Generate a collection of each category</hands-on-title>
>
> 1. Create a collection list with all these counts that you label `all counts`. Name each item so it only has the GSM id, the treatment and the library, for example, `GSM461176_untreat_single`.
>
>    {% snippet faqs/galaxy/collections_build_list.md %}
>
> 2. {% tool [Extract element identifiers](toolshed.g2.bx.psu.edu/repos/iuc/collection_element_identifiers/collection_element_identifiers/0.0.2) %} with the following parameters:
>    - {% icon param-collection %} *"Dataset collection"*: `all counts`
>
>    We will now split the collection by treatment. We need to find a pattern which will be present into only one of the 2 categories. We will use the word `untreat`:
>
> 3. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.3+galaxy1) %} (grep) with the following parameters:
>    - *"Select lines from"*: `Extract element identifiers on data XXX` (output of  **Extract element identifiers** {% icon tool %})
>    - *"that"*: `Match`
>    - *"Regular Expression"*: `untreat`
>
> 4. {% tool [Filter collecion](__FILTER_FROM_FILE__) %} with the following parameters:
>    - *"Input collection"*: `all counts`
>    - *"How should the elements to remove be determined"*: `Remove if identifiers are ABSENT from file`
>        - *"Filter out identifiers absent from"*: `Search in textfiles on data XXX` (output of  **Search in textfiles** {% icon tool %})
>
> 5. Rename both collections `untreated` (the filtered collection) and `treated` (the discarded collection).
>
> We will repeat the same process using `single`
>
> 6. {% tool [Search in textfiles](toolshed.g2.bx.psu.edu/repos/bgruening/text_processing/tp_grep_tool/9.3+galaxy1) %} (grep) with the following parameters:
>    - *"Select lines from"*: `Extract element identifiers on data XXX` (output of  **Extract element identifiers** {% icon tool %})
>    - *"that"*: `Match`
>    - *"Regular Expression"*: `single`
>
> 7. {% tool [Filter collecion](__FILTER_FROM_FILE__) %} with the following parameters:
>    - *"Input collection"*: `all counts`
>    - *"How should the elements to remove be determined"*: `Remove if identifiers are ABSENT from file`
>        - *"Filter out identifiers absent from"*: `Search in textfiles on data XXX` (output of  **Search in textfiles** {% icon tool %})
>
> 8. Rename both collections `single` (the filtered collection) and `paired` (the discarded collection).
{: .hands_on}

We can now run **DESeq2**:

> <hands-on-title>Determine differentially expressed features</hands-on-title>
>
> 1. {% tool [DESeq2](toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/2.11.40.8+galaxy0) %} with the following parameters:
>    - *"how"*: `Select datasets per level`
>        - In *"Factor"*:
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Treatment`
>           - In *"1: Factor level"*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `treated`
>               - {% icon param-collection %} *"Count file(s)"*: Select the collection `treated`
>           - In *"2: Factor level"*:
>               - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `untreated`
>               - {% icon param-collection %} *"Count file(s)"*: Select the collection `untreated`
>       - {% icon param-repeat %} *"Insert Factor"*
>           - *"Specify a factor name, e.g. effects_drug_x or cancer_markers"*: `Sequencing`
>               - In *"Factor level"*:
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `PE`
>                        - {% icon param-collection %} *"Count file(s)"*: Select the collection `paired`
>                    - {% icon param-repeat %} *"Insert Factor level"*
>                        - *"Specify a factor level, typical values could be 'tumor', 'normal', 'treated' or 'control'"*: `SE`
>                        - {% icon param-collection %} *"Count file(s)"*: Select the collection `single`
>    - *"Files have header?"*: `Yes`
>    - *"Choice of Input data"*: `Count data (e.g. from HTSeq-count, featureCounts or StringTie)`
>    - In *"Advanced options"*:
>        - *"Use beta priors"*: `Yes`
>    - In *"Output options"*:
>        - *"Output selector"*: `Generate plots for visualizing the analysis results`, `Output normalised counts`
>
{: .hands_on}
</div>

**DESeq2** generated 3 outputs:

- A table with the normalized counts for each gene (rows) in the samples (columns)
- A graphical summary of the results, useful to evaluate the quality of the experiment:

    1. A plot of the first 2 dimensions from a principal component analysis ([PCA](https://en.wikipedia.org/wiki/Principal_component_analysis)), run on the normalized counts of the samples

        > <details-title>What is a PCA?</details-title>
        >
        > Let's imagine we have some beer bottles standing here on the table. We can describe each beer by its colour, its foam, by how strong it is, and so on. We can compose a whole list of different characteristics of each beer in a beer shop. But many of them will measure related properties and so will be redundant. If so, we should be able to summarize each beer with fewer characteristics. This is what PCA or principal component analysis does.
        >
        > With PCA, we do not just select some interesting characteristics and discard the others. Instead, we construct some new characteristics that summarize our list of beers well. These new characteristics are constructed using the old ones. For example, a new characteristic might be computed, e.g. foam size minus beer pH. They are linear combinations.
        >
        > In fact, PCA finds the best possible characteristics, the ones that summarize the list of beers. These characteristics can be used to find similarities between beers and group them.
        >
        > Going back to read counts, the PCA is run on the normalized counts for all the samples. Here, we would like to describe the samples based on the expression of the genes. So the characteristics are the number of reads mapped on each genes. We use them and linear combinations of them to represent the samples and their similarities.
        >
        > *The beer analogy has been adapted from [an answer on StackExchange](https://stats.stackexchange.com/questions/2691/making-sense-of-principal-component-analysis-eigenvectors-eigenvalues)*.
        >
        {: .details}

        It shows the samples in the 2D plane spanned by their first two principal components. Each replicate is plotted as an individual data point. This type of plot is useful for visualizing the overall effect of experimental covariates and batch effects.

        > <question-title></question-title>
        >
        > ![DESeq PCA](../../images/ref-based/deseq2_pca.png "Principal component plot of the samples")
        >
        > 1. What is the first dimension (PC1) separating?
        > 2. And the second dimension (PC2)?
        > 3. What can we conclude about the DESeq design (factors, levels) we choose?
        >
        > > <solution-title></solution-title>
        > >
        > > 1. The first dimension is separating the treated samples from the untreated sample.
        > > 2. The second dimension is separating the single-end datasets from the paired-end datasets.
        > > 3. The datasets are grouped following the levels of the two factors. No hidden effect seems to be present on the data. If there is unwanted variation present in the data (e.g. batch effects), it is always recommended to correct for this, which can be achieved in DESeq2 by including in the design any known batch variables.
        > {: .solution}
        {: .question}

    2. Heatmap of the sample-to-sample distance matrix (with clustering) based on the normalized counts.

        The heatmap gives an overview of similarities and dissimilarities between samples: the color represents the distance between the samples. Dark blue means shorter distance, i.e. closer samples given the normalized counts.

        > <question-title></question-title>
        >
        > ![Heatmap of the sample-to-sample distances](../../images/ref-based/deseq2_sample_sample_distance_heatmap.png "Heatmap of the sample-to-sample distances")
        >
        > How are the samples grouped?
        >
        > > <solution-title></solution-title>
        > >
        > > They are first grouped by the treatment (the first factor) and secondly by the sequencing type (the second factor), as in the PCA plot.
        > >
        > {: .solution}
        {: .question}

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

> <question-title></question-title>
>
> 1. Is the FBgn0003360 gene differentially expressed because of the treatment? If yes, how much?
> 2. Is the *Pasilla* gene (ps, FBgn0261552) downregulated by the RNAi treatment?
> 3. We could also hypothetically be interested in the effect of the sequencing (or other secondary factors in other cases). How would we know the differentially expressed genes because of sequencing type?
> 4. We would like to analyze the interaction between the treatment and the sequencing. How could we do that?
>
> > <solution-title></solution-title>
> >
> > 1. FBgn0003360 is differentially expressed because of the treatment: it has a significant adjusted p-value ($$2.8 \cdot 10^{-171} << 0.05$$). It is less expressed (`-` in the log2FC column) in treated samples compared to untreated samples, by a factor ~8 ($$2^{log2FC} = 2^{2.99542318410271}$$).
> >
> > 2. You can manually check for the `FBgn0261552` in the first column or run {% tool [Filter data on any column using simple expressions](Filter1) %}
> >   - {% icon param-file %} *"Filter"*: the `DESeq2 result file` (output of **DESeq2** {% icon tool %})
> >   - *"With following condition"*: `c1 == "FBgn0261552"`
> >
> >    The log2 fold-change is negative so it is indeed downregulated and the adjusted p-value is below 0.05 so it is part of the significantly changed genes.
> >
> > 3. DESeq2 in Galaxy returns the comparison between the different levels for the 1st factor, after
> >    correction for the variability due to the 2nd factor. In our current case, treated against untreated for any sequencing type. To compare sequencing types, we should run DESeq2 again switching factors: factor 1 (treatment) becomes factor 2 and factor 2 (sequencing) becomes factor 1.
> > 4. To add the interaction between two factors (e.g. treated for paired-end data vs untreated for single-end), we should run DESeq2 another time but with only one factor with the following 4 levels:
> >    - treated-PE
> >    - untreated-PE
> >    - treated-SE
> >    - untreated-SE
> >
> >    By selection *"Output all levels vs all levels of primary factor (use when you have >2 levels for primary factor)"* to `Yes`, we can then compare treated-PE vs untreated-SE.
> >
> {: .solution}
{: .question}

## Annotation of the DESeq2 results

The ID for each gene is something like FBgn0003360, which is an ID from the corresponding database, here Flybase ({% cite thurmond2018flybase %}). These IDs are unique but sometimes we prefer to have the gene names, even if they may not reference an unique gene (e.g. duplicated after re-annotation). But gene names may hint already to a function or they help you to search for desired candidates. We would also like to display the location of these genes within the genome. We can extract such information from the annotation file which we used for mapping and counting.

> <hands-on-title>Annotation of the DESeq2 results</hands-on-title>
>
> 1. Import the Ensembl gene annotation for *Drosophila melanogaster* (`Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`) from the previous history, or from the Shared Data library or from Zenodo:
>
>    ```text
>    {{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz
>    ```
>
> 2. {% tool [Annotate DESeq2/DEXSeq output tables](toolshed.g2.bx.psu.edu/repos/iuc/deg_annotate/deg_annotate/1.1.0) %} with:
>    - {% icon param-file %} *"Tabular output of DESeq2/edgeR/limma/DEXSeq"*: the `DESeq2 result file` (output of **DESeq2** {% icon tool %})
>    - *"Input file type"*: `DESeq2/edgeR/limma`
>    - {% icon param-file %} *"Reference annotation in GFF/GTF format"*: imported gtf `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
>
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
> 1. Where is the most over-expressed gene located?
> 2. What is the name of the gene?
> 3. Where is the *Pasilla* gene located (FBgn0261552)?
>
> > <solution-title></solution-title>
> >
> > 1. FBgn0025111 (the top-ranked gene with the highest positive log2FC value) is located on the reverse strand of chromosome X, between 10,778,953 bp and 10,786,907 bp.
> > 2. From the table, we got the gene symbol: Ant2. After some search on the [online biological databases](https://www.ncbi.nlm.nih.gov/gene/32008), we find that Ant2 corresponds to adenine nucleotide translocase 2.
> > 3. The *Pasilla* gene is located on the forward strand of chromosome 3R, between 9,417,939 bp and 9,455,500 bp.
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
>
> 2. {% tool [Concatenate datasets](cat1) %} to add this header line to the **Annotate** output:
>    - {% icon param-file %} *"Concatenate Dataset"*: the `header` dataset
>    - *"Dataset"*
>       - Click on {% icon param-repeat %} *"Insert Dataset"*
>         - {% icon param-file %} *"select"*: output of **Annotate** {% icon tool %}
>
> 3. Rename the output to `Annotated DESeq2 results`
{: .hands_on}

## Extraction and annotation of differentially expressed genes

Now we would like to extract the most differentially expressed genes due to the treatment with a fold change > 2 (or < 1/2).

> <hands-on-title>Extract the most differentially expressed genes</hands-on-title>
>
> 1. {% tool [Filter data on any column using simple expressions](Filter1) %} to extract genes with a significant change in gene expression (adjusted *p*-value below 0.05) between treated and untreated samples:
>    - {% icon param-file %} *"Filter"*: `Annotated DESeq2 results`
>    - *"With following condition"*: `c7<0.05`
>    - *"Number of header lines to skip"*: `1`
>
> 2. Rename the output `Genes with significant adj p-value`.
>
>    > <question-title></question-title>
>    >
>    > How many genes have a significant change in gene expression between these conditions?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > We get 966 (967 lines including a header) genes (4.04%) with a significant change in gene expression between treated and untreated samples.
>    > >
>    > {: .solution}
>    {: .question}
>    >
>    > <comment-title></comment-title>
>    >
>    > The file with the independently filtered results can be used for further downstream analysis as it excludes genes with only a few read counts, as these genes will not be considered as significantly differentially expressed.
>    {: .comment}
>
>    We will now select only the genes with a fold change (FC) > 2 or FC < 0.5. Note that the DESeq2 output file contains $$log_{2} FC$$, rather than FC itself, so we filter for $$abs(log_{2} FC) > 1$$ (which implies FC > 2 or FC < 0.5).
>
> 3. {% tool [Filter data on any column using simple expressions](Filter1) %} to extract genes with an $$abs(log_{2} FC) > 1$$:
>    - {% icon param-file %} *"Filter"*: `Genes with significant adj p-value`
>    - *"With following condition"*: `abs(c3)>1`
>    - *"Number of header lines to skip"*: `1`
>
> 4. Rename the output `Genes with significant adj p-value & abs(log2(FC)) > 1`.
>
>    > <question-title></question-title>
>    >
>    > 1. How many genes have been conserved?
>    > 2. Can the *Pasilla* gene (ps, FBgn0261552) be found in this table?
>    >
>    > > <solution-title></solution-title>
>    > >
>    > > 1. We get 113 genes (114 lines including a header), or 11.79% of the significantly differentially expressed genes.
>    > > 2. The *Pasilla* gene can be found with a quick Search (or even using  {% tool [Filter data on any column using simple expressions](Filter1) %} )
>    > {: .solution}
>    {: .question}
>
{: .hands_on}

We now have a table with 113 lines and a header corresponding to the most differentially expressed genes. For each gene, we have its ID, its mean normalized counts (averaged over all samples from both conditions), its $$log_{2} FC$$ and other information including gene name and position.

## Visualization of the expression of the differentially expressed genes

We could plot the $$log_{2} FC$$ for the extracted genes, but here we would like to look at a heatmap of expression for these genes in the different samples. So we need to extract the normalized counts for these genes.

We proceed in several steps:

- Extract and plot the normalized counts for these genes for each sample with a heatmap, using the normalized count file generated by DESeq2
- Compute, extract and plot the Z-score of the normalized counts

> <comment-title>Advanced tutorials on visualization</comment-title>
>
> In this tutorial, we will quickly explain some possible visualization. For more details, please have a look in the extra tutorials on visualization of RNA-Seq results:
>
> - [Visualization of RNA-Seq results with heatmap2]({% link topics/transcriptomics/tutorials/rna-seq-viz-with-heatmap2/tutorial.md %})
> - [Visualization of RNA-Seq results with Volcano Plot]({% link topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot/tutorial.md %})
>
{: .comment}

### Visualization of the normalized counts

To extract the normalized counts for the interesting genes, we join the normalized count table generated by DESeq2 with the table we just generated. We will then keep only the columns corresponding to the normalized counts.

> <hands-on-title>Extract the normalized counts of the most differentially expressed genes</hands-on-title>
>
> 1. {% tool [Join two Datasets side by side on a specified field](join1) %} with the following parameters:
>    - {% icon param-file %} *"Join"*: the `Normalized counts` file (output of **DESeq2** {% icon tool %})
>    - *"using column"*: `Column: 1`
>    - {% icon param-file %} *"with"*: `Genes with significant adj p-value & abs(log2(FC)) > 1`
>    - *"and column"*: `Column: 1`
>    - *"Keep lines of first input that do not join with second input"*: `No`
>    - *"Keep the header lines"*: `Yes`
>
>    The generated file has more columns than we need for the heatmap: mean normalized counts, $$log_{2} FC$$ and other annotation information. We need to remove the extra columns.
>
> 2. {% tool [Cut](Cut1) %} columns from a table with the following parameters to extract the columns with the gene IDs and normalized counts:
>    - *"Cut columns"*: `c1-c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the joined dataset (output of **Join two Datasets** {% icon tool %})
>
> 3. Rename the output to `Normalized counts for the most differentially expressed genes`
>
{: .hands_on}

We now have a table with 114 lines (the 113 most differentially expressed genes and a header) and the normalized counts for these genes across the 7 samples.

> <hands-on-title>Plot the heatmap of the normalized counts of these genes for the samples</hands-on-title>
>
> 1. {% tool [heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.1.3.1+galaxy0) %} to plot the heatmap:
>    - {% icon param-file %} *"Input should have column headers"*: `Normalized counts for the most differentially expressed genes`
>    - *"Data transformation"*: `Log2(value+1) transform my data`
>    - *"Enable data clustering"*: `Yes`
>    - *"Labeling columns and rows"*: `Label columns and not rows`
>    - *"Type of colormap to use"*: `Gradient with 2 colors`
>
{: .hands_on}

You should obtain something similar to:

![Heatmap with the normalized counts for the most differentially expressed genes](../../images/ref-based/heatmap2_normalized_counts.png "Normalized counts for the most differentially expressed genes")

> <question-title></question-title>
>
> 1. What does the X-axis of the heatmap represent? What about the Y axis?
> 2. Do you observe anything in the clustering of the samples and the genes?
> 3. What changes if you regenerate the heatmap, this time selecting `Plot the data as it is` in *"Data transformation"*?
> 4. Why cannot we use `Log2(value) transform my data` in *"Data transformation"*?
> 5. How could you generate a heatmap of normalized counts for all up-regulated genes with fold change > 2?
>
> > <solution-title></solution-title>
> >
> > 1. The X-axis shows the 7 samples, together with a dendrogram representing the similarity between their levels of gene expression. The Y-axis shows the 113 differentially expressed genes, likewise with a dendrogram representing the similarity between the levels of gene expression.
> > 2. The samples are clustering by treatment.
> > 3. The scale changes and we only see few genes.
> > 4. Because the normalized expression of the gene `FBgn0013688` in `GSM461180_treat_paired` is at `0`.
> > 5. Extract the genes with $$log_{2} FC$$ > 1 (filter for genes with `c3>1` on the summary of the differentially expressed genes) and run **heatmap2** {% icon tool %} on the generated table.
> {: .solution}
{: .question}

### Visualization of the Z-score

To compare the gene expression over samples, we could also use the Z-score, which is often represented in publications.

The Z-score gives the number of standard-deviations that a value is away from the mean of all the values in the same group, here the same gene. A Z-score of -2 for the gene X in sample A means that this value is 2 standard-deviations lower than the mean of the values for gene X in all the samples (A, B, C, etc).

The Z-score $$z_{i,j}$$ for a gene $$i$$ in a sample $$j$$ given the normalized count $$x_{i,j}$$ is computed as $$z_{i,j} = \frac{x_{i,j}- \overline{x_i}}{s_i}$$ with $$\overline{x_i}$$ the mean and $$s_i$$ the standard deviation of the normalized counts for the gene $$i$$ over all samples.

> <details-title>Compute the Z-score for all genes</details-title>
>
> We often need the Z-score for some visualisations. To compute the Z-score, we break the process into 2 steps:
>
> 1. Substract each value by the mean of values in the row (i.e. $$x_{i,j}- \overline{x_i}$$) using the normalized count table
> 2. Divide the previous values by the standard deviation of values of row, using 2 tables (the normalized counts and the table computed in the  previous step)
>
> > <hands-on-title>Compute the Z-score of all genes</hands-on-title>
> >
> > 1. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters to >  first substract the mean values per row
> >    - *"Input Single or Multiple Tables"*: `Single Table`
> >      - {% icon param-file %} *"Table"*: `Normalized counts file on data ... and others` (output of **DESeq2** {% icon tool %})
> >      - *"Type of table operation"*: `Perform a full table operation`
> >        - *"Operation"*: `Custom`
> >          - *"Custom expression on 'table', along 'axis' (0 or 1)"*: `table.sub(table.mean(1), 0)`
> >
> >            The `table.mean(1)` expression computes the mean for each row (here the genes) and `table.sub(table.mean(1), 0)` substracts each value by the mean of the row (computed with `table.mean(1)`)
> >
> > 2. {% tool [Table Compute](toolshed.g2.bx.psu.edu/repos/iuc/table_compute/table_compute/1.2.4+galaxy0) %} with the following parameters:
> >    - *"Input Single or Multiple Tables"*: `Multiple Table`
> >      - Click on {% icon param-repeat %} *"Insert Tables"*
> >      - In *"1: Tables"*:
> >        - {% icon param-file %} *"Table"*: `Normalized counts file on data ... and others` (output of **DESeq2** {% icon tool %})
> >      - Click on {% icon param-repeat %} *"Insert Tables"*
> >      - In *"2: Tables"*:
> >        - {% icon param-file %} *"Table"*: output of the first **Table Compute** {% icon tool %}
> >      - *"Custom expression on 'tableN'"*: `table2.div(table1.std(1),0)`
> >
> >        The `table1.std(1)` expression computes the standard deviations of each row on the 1st table (normalized counts) and `table2.div` divides the values of 2nd table (previously computed) by these standard deviations.
> >
> > 3. Rename the output to `Z-scores`
> > 4. Inspect the output file
> {: .hands_on}
>
> We now have a table with the Z-score for all genes in the 7 samples.
>
> > <question-title></question-title>
> >
> > 1. What is the range for the Z-score?
> > 2. Why some rows are empty?
> > 3. What can we say about the Z-scores for the differentially expressed genes (for example, `FBgn0037223`)?
> > 4. Can we use the Z-score to estimate the strength of the differential expression of a gene?
> >
> > > <solution-title></solution-title>
> > >
> > > 1. The Z-score ranges from -3 standard deviations up to +3 standard deviations. It can be placed on a normal distribution curve: -3 being the far left of the normal distribution curve and +3 the far right of the normal distribution curve
> > > 2. If all counts are identicals (usually to 0), the standard deviation is 0, Z-score cannot be computed for these genes.
> > > 3. When a gene is differentially expressed between two groups (here treated and untreated), the Z-scores for this gene will be (mostly) positive for the samples in one group and (mostly) negative for the samples in the other group.
> > > 4. The Z-score is a signal-to-noise ratio. Large absolute Z-scores, i.e. large positive or negative values, is not a direct estimate of the effect, i.e. the strength of the differential expression. A same large Z-score can have different meanings, depending on the noise:
> > >    - with large noise: a very large effect
> > >    - with some noise: a rather large effect
> > >    - with only little noise: a rather small effect
> > >    - with almost no noise: a tiny effect
> > >
> > >    The problem is that "noise" here is not only noise from the measure. It can also be linked to the "tightness" of the gene regulation control. Not tightly controlled genes, i.e. whose expression may vary in a wide range over samples, can be considerably induced or repressed. Their absolute Z-score will be small as the variations over samples is big. In contrast, genes that are tightly controlled may have only very small changes in their expression, without any biological impact. The absolute Z-score will be large for these genes.
> > >
> > {: .solution}
> {: .question}
{: .details}

We would like now to plot a heatmap for the Z-scores:

![Heatmap with the Z-score counts for the most differentially expressed genes](../../images/ref-based/z-score-heatmap.png "Z-scores for the most differentially expressed genes")

> <hands-on-title>Plot the Z-score of the most differentially expressed genes</hands-on-title>
>
> 1. {% tool [heatmap2](toolshed.g2.bx.psu.edu/repos/iuc/ggplot2_heatmap2/ggplot2_heatmap2/3.1.3.1+galaxy0) %} to plot the heatmap:
>    - {% icon param-file %} *"Input should have column headers"*: `Normalized counts for the most differentially expressed genes`
>    - *"Data transformation"*: `Plot the data as it is`
>    - *"Compute z-scores prior to clustering"*: `Compute on rows`
>    - *"Enable data clustering"*: `Yes`
>    - *"Labeling columns and rows"*: `Label columns and not rows`
>    - *"Type of colormap to use"*: `Gradient with 3 colors`
{: .hands_on}

# Functional enrichment analysis of the DE genes

We have extracted genes that are differentially expressed in treated (PS gene-depleted) samples compared to untreated samples. Now, we would like to know if the differentially expressed genes are enriched transcripts of genes which belong to more common or specific categories in order to identify biological functions that might be impacted.

## Gene Ontology analysis

[Gene Ontology (GO)](http://www.geneontology.org/) analysis is widely used to reduce complexity and highlight biological processes in genome-wide expression studies. However, standard methods give biased results on RNA-Seq data due to over-detection of differential expression for long and highly-expressed transcripts.

[**goseq**](https://bioconductor.org/packages/release/bioc/vignettes/goseq/inst/doc/goseq.pdf) ({% cite young2010gene %}) provides methods for performing GO analysis of RNA-Seq data while taking length bias into account. **goseq** could also be applied to other category-based tests of RNA-Seq data, such as KEGG pathway analysis, as discussed in a further section.

**goseq** needs 2 files as inputs:

- A tabular file with the differentially expressed genes from all genes assayed in the RNA-Seq experiment with 2 columns:
  - the Gene IDs (unique within the file), in uppercase letters
  - a boolean indicating whether the gene is differentially expressed or not (`True` if differentially expressed or `False` if not)
- A file with information about the length of a gene to correct for potential length bias in differentially expressed genes

> <hands-on-title>Prepare the first dataset for goseq</hands-on-title>
>
> 1. {% tool [Compute](toolshed.g2.bx.psu.edu/repos/devteam/column_maker/Add_a_column1/2.0) %} on rows with the following parameters:
>    - {% icon param-file %} *"Input file"*: the `DESeq2 result file` (output of **DESeq2** {% icon tool %})
>    - In *"Expressions"*:
>      - {% icon param-text %} *"Add expression"*: `bool(float(c7)<0.05)`
>      - {% icon param-select %} *"Mode of the operation?"*: `Append`
>    - Under *"Error handling"*:
>      - {% icon param-toggle %} *"Autodetect column types"*: `No`
>      - {% icon param-select %} *"If an expression cannot be computed for a row"*: `Fill in a replacement value`
>      - {% icon param-select %} *"Replacement value"*: `False`
>
> 2. {% tool [Cut](Cut1) %} columns from a table with the following parameters:
>    - *"Cut columns"*: `c1,c8`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: the output of the **Compute** {% icon tool %}
>
> 3. {% tool [Change Case](ChangeCase) %} with
>    - {% icon param-file %} *"From"*: the output of the previous **Cut** {% icon tool %}
>    - *"Change case of columns"*: `c1`
>    - *"Delimited by"*: `Tab`
>    - *"To"*: `Upper case`
>
> 4. Rename the output to `Gene IDs and differential expression`
{: .hands_on}

We just generated the first input for **goseq**. As second input for **goseq** we need the gene lengths. We can use here the gene lengths generated by **featureCounts**  or **Gene length and GC content** and format the gene IDs.

> <hands-on-title>Prepare the gene length file</hands-on-title>
>
> <div class="featureCounts" markdown="1">
> 1. Copy the feature length collection previously generated by **featureCounts** {% icon tool %} into this history
>
>    {% snippet faqs/galaxy/histories_copy_dataset.md %}
>
> 2. {% tool [Extract Dataset](__EXTRACT_DATASET__) %} with:
>    - {% icon param-collection %} *"Input List"*: `featureCounts on collection N: Feature lengths`
>    - *"How should a dataset be selected?"*: `The first dataset`
>
> </div>
>
> <div class="STAR" markdown="1">
> 1. Copy the output of **Gene length and GC content** {% icon tool %} (`Gene length`) into this history
>
>    {% snippet faqs/galaxy/histories_copy_dataset.md %}
> </div>
>
> 2. {% tool [Change Case](ChangeCase) %} with the following parameters:
>
>    - {% icon param-file %} *"From"*: <span class="featureCounts" markdown="1">`GSM461177_untreat_paired` (output of **Extract Dataset** {% icon tool %})</span><span class="STAR" markdown="1">`Gene length`</span>
>    - *"Change case of columns"*: `c1`
>    - *"Delimited by"*: `Tab`
>    - *"To"*: `Upper case`
>
> 3. Rename the output to `Gene IDs and length`
{: .hands_on}

We have now the two required input files for goseq.

> <hands-on-title>Perform GO analysis</hands-on-title>
>
> 1. {% tool [goseq](toolshed.g2.bx.psu.edu/repos/iuc/goseq/goseq/1.50.0+galaxy0) %} with
>    - *"Differentially expressed genes file"*: `Gene IDs and differential expression`
>    - *"Gene lengths file"*: `Gene IDs and length`
>    - *"Gene categories"*: `Get categories`
>       - *"Select a genome to use"*: `Fruit fly (dm6)`
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
    > 2. How are the over-represented GO terms divided into MF, CC and BP? And for under-represented GO terms?
    >
    > > <solution-title></solution-title>
    > >
    > > 1. 60 GO terms (0.50%) are over-represented and 7 (0.07%) under-represented.
    > >
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
> In this tutorial, we covered GO enrichment analysis with **goseq**. To learn other methods and tools for gene set enrichment analysis, please have a look at the ["RNA-Seq genes to pathways"]({% link topics/transcriptomics/tutorials/rna-seq-genes-to-pathways/tutorial.md %}) tutorial.
{: .comment}

## KEGG pathways analysis

**goseq** can also be used to identify interesting KEGG pathways. The KEGG pathway database is a collection of pathway maps representing current knowledge of molecular interaction, reaction and relation networks. A map can integrate many entities including genes, proteins, RNAs, chemical compounds, glycans, and chemical reactions, as well as disease genes and drug targets.

For example, the pathway `dme00010` represents the glycolysis process (conversion of glucose into pyruvate with generation of small amounts of ATP and NADH) for Drosophila melanogaster:

![dme00010 KEGG pathway](../../images/ref-based/dme00010_empty.png)

> <hands-on-title>Perform KEGG pathway analysis</hands-on-title>
>
> 1. {% tool [goseq](toolshed.g2.bx.psu.edu/repos/iuc/goseq/goseq/1.50.0+galaxy0) %} with
>    - *"Differentially expressed genes file"*: `Gene IDs and differential expression`
>    - *"Gene lengths file"*: `Gene IDs and length`
>    - *"Gene categories"*: `Get categories`
>       - *"Select a genome to use"*: `Fruit fly (dm6)`
>       - *"Select Gene ID format"*: `Ensembl Gene ID`
>       - *"Select one or more categories"*: `KEGG`
>    - In *"Output Options"*
>      - *"Output Top GO terms plot?"*: `No`
>      - *"Extract the DE genes for the categories (GO/KEGG terms)?"*: `Yes`
>
{: .hands_on}

**goseq** generates with these parameters 2 outputs:

1. A large table with the KEGG terms and some statistics

    > <question-title></question-title>
    >
    > 1. How many KEGG pathways terms have been identified?
    > 2. How many KEGG pathways terms are over-represented with an adjusted P value < 0.05?
    > 3. What are the over-represented KEGG pathways terms?
    > 4. How many KEGG pathways terms are under-represented with an adjusted P value < 0.05?
    >
    > > <solution-title></solution-title>
    > >
    > > 1. The file has 128 lines including an header, so 127 KEGG pathways have been identified.
    > > 2. 2 KEGG pathways (2.34%) are over-represented, using {% tool [Filter data on any column using simple expressions](Filter1) %} on c6 (adjusted p-value for over-represented KEGG pathways)
    > > 3. The 2 KEGG pathways over-represented are `01100` and `00010`. By searching on the [KEGG database](https://www.genome.jp/kegg/kegg2.html) for them, we can find more information about these pathways: `01100` corresponds to all metabolic pathways and `00010` to pathway for Glycolysis / Gluconeogenesis.
    > > 4. No KEGG pathway is under-represented, using {% tool [Filter data on any column using simple expressions](Filter1) %} on c7 (adjusted p-value for under-represented KEGG pathways)
    > {: .solution}
    {: .question}

2. A table with the differentially expressed genes (from the list we provided) associated with the KEGG pathways (`DE genes for categories (GO/KEGG terms)`)

We could investigate which genes are involved in which pathways by looking at the second file generated by **goseq**. However, this can be cumbersome and we would like to see the pathways as represented in the previous image. **Pathview** ({% cite luo2013pathview %}) can help to automatically generate similar images to the previous one, while also adding extra information about the genes (e.g. expression) in our study.

This tool needs 2 main inputs:

- Pathway ID(s) to plot, either as just one ID or as a file with one column with the pathway IDs
- A tabular file with the genes in the RNA-Seq experiment with 2 (or more) columns:
  - the gene IDs (unique within the file)
  - some information about the genes

    This can be for example a p-value or a fold change. This information will be added to the pathway plot: the node of the corresponding gene will be colored given the value. If there are different columns, the different information will be plotted side by side on the node.

Here we would like to visualize the 2 KEGG pathways: the over-represented `00010` (Glycolysis / Gluconeogenesis) and the most under-represented (yet not significantly) `03040` (Spliceosome). We would like the gene nodes to be colored by Log2 Fold Change for the differentially expressed genes because of the treatment.

> <hands-on-title>Overlay log2FC on KEGG pathway</hands-on-title>
>
> 1. {% tool [Cut](Cut1) %} columns from a table with the following parameters:
>    - *"Cut columns"*: `c1,c3`
>    - *"Delimited by"*: `Tab`
>    - {% icon param-file %} *"From"*: `Genes with significant adj p-value`
>
> 2. Rename to `Genes with significant adj p-value and their Log2 FC`
>
>    We extracted the ID and Log2 Fold Change for the genes that have a significant adjusted p-value.
>
> 3. Create a new tabular file from the following (IDs of the pathways to plot) named `KEGG pathways to plot`
>
>    ```text
>    00010
>    03040
>    ```
>
> 4. {% tool [Pathview](toolshed.g2.bx.psu.edu/repos/iuc/pathview/pathview/1.34.0+galaxy0) %} with
>    - *"Number of pathways to plot"*: `Multiple`
>      - {% icon param-file %} *"KEGG pathways"*: `KEGG pathways to plot`
>      - *"Does the file have header (a first line with column names)?"*: `No`
>    - *"Species to use"*: `Fly`
>    - *"Provide a gene data file?"*: `Yes`
>      - {% icon param-file %} *"Gene data"*: `Genes with significant adj p-value and their Log2 FC`
>      - *"Does the file have header (a first line with column names)?"*: `Yes`
>      - *"Format for gene data"*: `Ensembl Gene ID`
>    - *"Provide a compound data file?"*: `No`
>    - In *"Output Options"*
>      - *"Output for pathway"*: `KEGG native`
>        - *"Plot on same layer?"*: `Yes`
{: .hands_on}

**Pathview** generates a collection with the KEGG visualization: one file per pathway.

> <question-title></question-title>
>
> `dme00010` KEGG pathway from **Pathview**
>
> ![KEGG pathway](../../images/ref-based/dme00010.png)
>
> 1. What are the colored boxes?
> 2. What is the color code?
>
> > <solution-title></solution-title>
> >
> > 1. The colored boxes are genes in the pathway that are differentially expressed
> > 2. Pay attention that the color code is counterintuitive: green is for value below 0, it means for genes with a log2FC < 0 and red for genes with a log2FC > 0.
> >
> {: .solution}
{: .question}

{% comment %}

# Inference of the differential exon usage

Next, we would like to know the differential exon usage between treated (PS depleted) and untreated samples using RNA-Seq exon counts. We will rework the mapping results we generated previously.

We will use [DEXSeq](https://www.bioconductor.org/packages/release/bioc/html/DEXSeq.html). DEXSeq detects high sensitivity genes, and in many cases exons, that are subject to differential exon usage. But first, as for the differential gene expression, we need to count the number of reads mapping to the exons.

## Count the number of reads per exon

This step is similar to the step of [counting the number of reads per annotated gene](#count-the-number-of-reads-per-annotated-gene) except that, instead of HTSeq-count, we are using DEXSeq-Count.

> <hands-on-title>Counting the number of reads per exon</hands-on-title>
>
> 1. {% tool [DEXSeq-Count](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq_count/1.28.1.0) %}: Use the **DEXSeq-Count** to prepare the *Drosophila* annotations to extract only exons with corresponding gene ids
>     - *"Mode of operation"*: `Prepare annotation`
>       - {% icon param-file %} *"GTF file"*: `Drosophila_melanogaster.BDGP6.32.109_UCSC.gtf.gz`
>
>    The output is again a GTF file that is ready to be used for counting
>
> 2. {% tool [DEXSeq-Count](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq_count/1.28.1.0) %}: Count reads using **DEXSeq-Count** with
>     - *"Mode of operation"*: `Count reads`
>       - {% icon param-files %} *"Input bam file"*: the `BAM` files generated by **RNA STAR**
>       - {% icon param-file %} *"DEXSeq compatible GTF file"*: the GTF file generated by **DEXSeq-Count**
>       - *"Is library paired end?"*: `Yes`
>       - *"Is library strand specific?*: `No`
>       - *"Skip all reads with alignment quality lower than the given minimum value"*:  `10`
>
{: .hands_on}

DEXSeq generates a count table similar to the one generated by featureCounts, but with counts for exons.

> <question-title></question-title>
>
> 1. Which exon has the most reads mapped to it for both samples?
> 2. Which gene does this exon belong to?
> 3. Is there a connection to the previous result obtained with featureCounts?
>
> > <solution-title></solution-title>
> >
> > FBgn0284245:005 is the exon with the most reads mapped to it for both samples. It is part of FBgn0284245, the feature with the most reads mapped on it (from featureCounts).
> >
> {: .solution}
>
{: .question}

## Differential exon usage

DEXSeq usage is similar to DESeq2. It uses similar statistics to find differentially used exons.

As for DESeq2, in the previous step, we counted only reads that mapped to exons on chromosome 4 and for only one sample. To be able to identify differential exon usage induced by PS depletion, all datasets (3 treated and 4 untreated) must be analyzed following the same procedure. To save time, we did that for you. The results are available on [Zenodo]({{ page.zenodo_link }}):

- The results of running DEXSeq-count in 'Prepare annotation' mode
- Seven count files generated in 'Count reads' mode

> <hands-on-title></hands-on-title>
>
> 1. Create a **new empty history**
>
>    {% snippet faqs/galaxy/histories_create_new.md %}
>
> 2. Import the seven count files from [Zenodo]({{ page.zenodo_link }}) or the Shared Data library (if available):
>
>    - `Drosophila_melanogaster.BDGP6.87.dexseq.gtf`
>    - `GSM461176_untreat_single.exon.counts`
>    - `GSM461177_untreat_paired.exon.counts`
>    - `GSM461178_untreat_paired.exon.counts`
>    - `GSM461179_treat_single.exon.counts`
>    - `GSM461180_treat_paired.exon.counts`
>    - `GSM461181_treat_paired.exon.counts`
>    - `GSM461182_untreat_single.exon.counts`
>
>    ```text
>    {{ page.zenodo_link }}/files/Drosophila_melanogaster.BDGP6.87.dexseq.gtf
>    {{ page.zenodo_link }}/files/GSM461176_untreat_single.exon.counts
>    {{ page.zenodo_link }}/files/GSM461177_untreat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461178_untreat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461179_treat_single.exon.counts
>    {{ page.zenodo_link }}/files/GSM461180_treat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461181_treat_paired.exon.counts
>    {{ page.zenodo_link }}/files/GSM461182_untreat_single.exon.counts
>    ```
>
> 3. {% tool [DEXSeq](toolshed.g2.bx.psu.edu/repos/iuc/dexseq/dexseq/1.28.1+galaxy1) %}: Run **DEXSeq** with
>    - {% icon param-file %} *"GTF file created from DEXSeq-Count tool"*: `Drosophila_melanogaster.BDGP6.87.dexseq.gtf`
>    - In *"Factor"*:
>       - In "1: Factor"
>           - *"Specify a factor name"*: `condition`
>           - In *"Factor level"*:
>               - In *"1: Factor level"*:
>                   - *"Specify a factor level"*: `treated`
>                   - {% icon param-files %} *"Counts file(s)"*: the 3 exon count files with `treated` in their name
>               - In *"2: Factor level"*:
>                   - *"Specify a factor level"*: `untreated`
>                   - {% icon param-files %} *"Counts file(s)"*: the 4 exon count files with `untreated` in their name
>       - Click on *"Insert Factor"* (not on "Insert Factor level")
>       - In "2: Factor"
>           - "Specify a factor name" to `sequencing`
>           - In *"Factor level"*:
>               - In *"1: Factor level"*:
>                   - *"Specify a factor level"*: `PE`
>                   - {% icon param-files %} *"Counts file(s)"*: the 4 exon count files with `paired` in their name
>               - In *"2: Factor level"*:
>                   - *"Specify a factor level"*: `SE`
>                   - {% icon param-files %} *"Counts file(s)"*: the 3 exon count files with `single` in their name
>
>    > <comment-title></comment-title>
>    >
>    > Unlike DESeq2, DEXSeq does not allow flexible primary factor names. Always use your primary factor name as "condition"
>    {: .comment}
{: .hands_on}

Similarly to DESeq2, DEXSeq generates a table with:

1. Exon identifiers
2. Gene identifiers
3. Exon identifiers in the Gene
4. Mean normalized counts, averaged over all samples from both conditions
5. Logarithm (to basis 2) of the fold change

   The log2 fold changes are based on primary factor level 1 vs. factor level 2. The order of factor levels is then important. For example, for the factor 'Condition', DESeq2 computes fold changes of 'treated' samples against 'untreated', *i.e.* the values correspond to up- or downregulations of genes in treated samples.

6. Standard error estimate for the log2 fold change estimate
7. *p*-value for the statistical significance of this change
8. *p*-value adjusted for multiple testing with the Benjamini-Hochberg procedure which controls false discovery rate ([FDR](https://en.wikipedia.org/wiki/False_discovery_rate))

> <hands-on-title></hands-on-title>
>
> 1. {% tool [Filter data on any column using simple expressions](Filter1) %} to extract exons with a significant differential usage (adjusted *p*-value equal or below 0.05) between treated and untreated samples
>
> > <question-title></question-title>
> >
> > How many exons show a significant change in usage between these conditions?
> >
> > > <solution-title></solution-title>
> > >
> > > We get 38 exons (12.38%) with a significant usage change between treated and untreated samples.
> > >
> > {: .solution}
> {: .question}
{: .hands_on}

{% endcomment %}

# Conclusion



In this tutorial, we have analyzed real RNA sequencing data to extract useful information, such as which genes are up or downregulated by depletion of the *Pasilla* gene, but also which GO terms or KEGG pathways they are involved in. To answer these questions, we analyzed RNA sequence datasets using a reference-based RNA-Seq data analysis approach. This approach can be summarized with the following scheme:

![Summary of the analysis pipeline used](../../images/ref-based/tutorial-scheme.png "Summary of the analysis pipeline used")
