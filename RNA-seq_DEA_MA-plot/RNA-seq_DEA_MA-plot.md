# Visualise differential expression analysis results with MA plot

I-Hsuan Lin

University of Manchester

February 08, 2021

## 1. Introduction

This notebook shows readers how to use `ggplot2` package to create various plots to show the magnitude of the change in expression between two groups (M) in comparison to the magnitude of average expression over all samples (A).

### About this dataset

The RNA-seq dataset [E-MTAB-8411](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411) contains 5 wild-type (GK1, GK3, GK5, GK7 and GK9) and 4 macrophage-specific *Bmal1* knockout samples (GK2, GK4, GK6 and GK10). The fastq files were retreived from [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-8411/samples/) and processed with the follow steps:

1. Adapter and QC trimming with `BBDuk`
2. Mapping to Mouse genome (assemble MM10 and [GENCODE M24](https://www.gencodegenes.org/mouse/release_M24.html) annotation) with `STAR`
3. Differential expression analysis with `DESeq2`. *The log2 fold change shrinkage was not applied.*

## 2. Loading required libraries


```R
library(data.table)
library(ggplot2)
library(ggrepel)
library(scales)
```

## 3. Set output parameters


```R
# Set width
options(width = 110)

# Set output image size
options(repr.plot.width = 12, repr.plot.height = 8, repr.plot.res = 150)
```

## 4. Retrieve example results

The `data.frame` should contain at least the average expression and log2 fold change (or fold change) information in order for one to produce MA plot.


```R
data <- fread("https://raw.githubusercontent.com/ycl6/Simple-R-gists/main/RNA-seq_DEA_MA-plot/DESeq2_DEG.txt")
head(data, 10)
```

## 5. Create plots

### Prepare plot settings


```R
# Set alpha
alpha = 0.05

# Evaluate if adjusted P-values is smaller than alpha 
data$isDE = ifelse(is.na(data$padj), FALSE, data$padj < alpha)
table(isDE = data$isDE)

# Subset result to genes with detectable expression
# Genes with 0 expression will not be able to show on log-scaled MA plot, removed to prevent warning message
data = subset(data, basemean != 0)

# Show log2 fold change distribution
print("log2 fold change distribution:")
round(quantile(data$log2fc),2)

# Use log2 fold change distribution to set an appropriate Y-axis limits
# Set limit at 99th quantile to ensure there are some genes exceeding the limits,
# and they will be shown as different shapes
ylim = quantile(abs(data$log2fc[is.finite(data$log2fc)]), probs = 0.99)
# Use floor() to set lower bound and create upper and lower limits
ylim = c(-1,1) * (floor(ylim * 2) / 2)

print(paste("Y-axis lower bound:", ylim[1]))
print(paste("Y-axis upper bound:", ylim[2]))

# Re-scale fold change to within ylim upper and lower limits
data$log2fc2 = pmax(ylim[1], pmin(ylim[2], data$log2fc))

# Determine shape type
# A: If smaller than lower limit - down-pointing triangle
# B: If larger than upper limit - up-pointing triangle
# C: Within upper and lower limits - round circle
data$shape = ifelse(data$log2fc < ylim[1], "A", ifelse(data$log2fc > ylim[2], "B", "C"))
```


```R
head(data)
```

### Create MA plot

Here, we use the *re-scaled* log2 fold change values for `y` aesthetics, and all the data points will be within the specified Y-axis limits. Data points that originally exceed the Y-axis limits will be shown as either up-pointing or down-pointing triangles.


```R
# Set log10 breaks and labels
log10break = trans_breaks("log10", function(x) 10^x)
log10label = trans_format("log10", math_format(10^.x))
                          
ggplot(data, aes(x = basemean, y = log2fc2, color = isDE, size = isDE, alpha = isDE, shape = shape)) + 
    geom_point() + theme_classic(base_size = 18) + 
    geom_hline(yintercept = 0, size = 1, color = "darkorange", linetype = "dashed") +
    scale_x_log10(breaks = log10break, labels = log10label) + 
    scale_y_continuous(limits = ylim) +
    scale_color_manual(name = "Is DEG", values = setNames(c("red", "gray50"), c(TRUE, FALSE))) +
    scale_alpha_manual(guide = FALSE, values = setNames(c(0.6, 0.3), c(TRUE, FALSE))) + 
    scale_size_manual(guide = FALSE, values = setNames(c(2, 1), c(TRUE, FALSE))) +
    scale_shape_manual(guide = FALSE, values = setNames(c(6, 2, 16), c("A", "B", "C"))) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    ggtitle("MA plot") + xlab("mean of normalized counts") + ylab("log2 fold change")
```


```R
# Show up-regulated genes beyong Y-axis limits
data[data$isDE == TRUE & data$log2fc > ylim[2],]

# Show down-regulated genes beyong Y-axis limits
data[data$isDE == TRUE & data$log2fc < ylim[1],]
```

Use `geom_text_repel` to label genes that exceed the Y-axis limits.


```R
ggplot(data, aes(x = basemean, y = log2fc2, color = isDE, size = isDE, alpha = isDE, shape = shape)) + 
    geom_point() + theme_classic(base_size = 18) + 
    geom_hline(yintercept = 0, size = 1, color = "darkorange", linetype = "dashed") +
    scale_x_log10(breaks = log10break, labels = log10label) + 
    scale_y_continuous(limits = ylim*1.01) +
    scale_color_manual(name = "Is DEG", values = setNames(c("red", "gray50"), c(TRUE, FALSE))) +
    scale_alpha_manual(guide = FALSE, values = setNames(c(0.6, 0.3), c(TRUE, FALSE))) + 
    scale_size_manual(guide = FALSE, values = setNames(c(2, 1), c(TRUE, FALSE))) +
    scale_shape_manual(guide = FALSE, values = setNames(c(6, 2, 16), c("A", "B", "C"))) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    geom_text_repel(aes(label = ifelse(isDE == TRUE & abs(log2fc) > ylim[2], GeneSymbol, "")), # add labels
                    color = "blue", seed = 123, min.segment.length = unit(0, "lines"), 
                    box.padding = 0.5, size = 4, max.overlaps = Inf) +
    ggtitle("MA plot") + xlab("mean of normalized counts") + ylab("log2 fold change")
```

## 6. Notes

Differential expression analysis frameworks such as [limma](https://bioconductor.org/packages/limma), [edgeR](https://bioconductor.org/packages/edgeR) and [DESeq2](https://bioconductor.org/packages/DESeq2) each have its own function to produce a MA plot, and requires the input to be an object created as part of the analysis workflow. For example, DESeq2's `plotMA` function requires an a `DESeqResults` object produced by `results`; or a `DESeqDataSet` processed by `DESeq`, or the individual functions `nbinomWaldTest` or `nbinomLRT`.

You can use these build-in functions to create MA plot as part of the analysis. By using `ggplot2`, it can give you more flexibility and allows you to customised the layout and appearance.

## Session info


```R
sessionInfo()
```
