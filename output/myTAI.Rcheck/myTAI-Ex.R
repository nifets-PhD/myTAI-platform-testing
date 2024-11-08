pkgname <- "myTAI"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "myTAI-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('myTAI')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CollapseReplicates")
### * CollapseReplicates

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CollapseReplicates
### Title: Combine Replicates in an ExpressionSet
### Aliases: CollapseReplicates

### ** Examples

data(PhyloExpressionSetExample)

# combine the expression levels of the 2 replicates (const) per stage
# using mean as window function and rename new stages: "S1","S2","S3"
CollapseReplicates(ExpressionSet = PhyloExpressionSetExample[1:5,1:8], 
                   nrep          = 2, 
                   FUN           = mean, 
                   stage.names   = c("S1","S2","S3"))

# combine the expression levels of the 2 replicates (stage one), 2 replicates (stage two),
# and 3 replicates (stage three) using mean as window function 
# and rename new stages: "S1","S2","S3"
CollapseReplicates(ExpressionSet = PhyloExpressionSetExample[1:5,1:9], 
                   nrep          = c(2,2,3), 
                   FUN           = mean, 
                   stage.names   = c("S1","S2","S3"))





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CollapseReplicates", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("CombinatorialSignificance")
### * CombinatorialSignificance

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: CombinatorialSignificance
### Title: Compute the Statistical Significance of Each Replicate
###   Combination
### Aliases: CombinatorialSignificance

### ** Examples

## Not run: 
##D # load a standard PhyloExpressionSet
##D data(PhyloExpressionSetExample)
##D 
##D # we assume that the PhyloExpressionSetExample 
##D # consists of 3 developmental stages 
##D # and 2 replicates for stage 1, 3 replicates for stage 2, 
##D # and 2 replicates for stage 3
##D # FOR REAL ANALYSES PLEASE USE: permutations = 1000 or 10000
##D # BUT NOTE THAT THIS TAKES MUCH MORE COMPUTATION TIME
##D p.vector <- CombinatorialSignificance(ExpressionSet = PhyloExpressionSetExample, 
##D                                       replicates    = c(2,3,2), 
##D                                       TestStatistic = "FlatLineTest", 
##D                                       permutations  = 1000, 
##D                                       parallel      = FALSE)
## End(Not run)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("CombinatorialSignificance", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("DiffGenes")
### * DiffGenes

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: DiffGenes
### Title: Differential Gene Expression Analysis
### Aliases: DiffGenes

### ** Examples


data(PhyloExpressionSetExample)

# Detection of DEGs using the fold-change measure
DEGs <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[ ,1:8],
                  nrep          = 2,
                  comparison    = "below",
                  method        = "foldchange",
                  stage.names   = c("S1","S2","S3"))


head(DEGs)


# Detection of DEGs using the log-fold-change measure
# when choosing method = "log-foldchange" it is assumed that
# your input expression matrix stores log2 expression levels 
log.DEGs <- DiffGenes(ExpressionSet = tf(PhyloExpressionSetExample[1:5,1:8],log2),
                      nrep          = 2,
                      comparison    = "below",
                      method        = "log-foldchange",
                      stage.names   = c("S1","S2","S3"))


head(log.DEGs)


# Remove fold-change values < 2 from the dataset:

## first have a look at the range of fold-change values of all genes 
apply(DEGs[ , 3:8],2,range)

# now remove genes undercutting the alpha = 2 threshold
# hence, remove genes having p-values <= 0.05 in at
# least one sample comparison
DEGs.alpha <- DiffGenes(ExpressionSet = PhyloExpressionSetExample[1:250 ,1:8],
                        nrep          = 2,
                        method        = "t.test",
                        alpha         = 0.05,
                        comparison    = "above",
                        filter.method = "n-set",
                        n             = 1,
                        stage.names   = c("S1","S2","S3"))

# now again have a look at the range and find
# that fold-change values of 2 are the min value
apply(DEGs.alpha[ , 3:5],2,range)

# now check whether each example has at least one stage with a p-value <= 0.05
head(DEGs.alpha)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("DiffGenes", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("EarlyConservationTest")
### * EarlyConservationTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EarlyConservationTest
### Title: Perform Reductive Early Conservation Test
### Aliases: EarlyConservationTest

### ** Examples


data(PhyloExpressionSetExample)

# perform the early conservation test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
# stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
# module 3 = late
EarlyConservationTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7), 
                       permutations = 1000)


# use your own permutation matrix based on which p-values (EarlyConservationTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

EarlyConservationTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7), 
                       custom.perm.matrix = custom_perm_matrix)
                       



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EarlyConservationTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("EnrichmentTest")
### * EnrichmentTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: EnrichmentTest
### Title: Phylostratum or Divergence Stratum Enrichment of a given Gene
###   Set based on Fisher's Test
### Aliases: EnrichmentTest

### ** Examples


data(PhyloExpressionSetExample)

set.seed(123)
test_set <- sample(PhyloExpressionSetExample[ , 2],1000)

E.Result <- EnrichmentTest(ExpressionSet = PhyloExpressionSetExample,
                           test.set      = test_set ,
                           measure       = "log-foldchange")
                           
# get the log-fold change table
E.Result$enrichment.matrix

# get P-values for the enrichment significance for each Phylostratum
E.Result$p.values




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("EnrichmentTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("Expressed")
### * Expressed

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Expressed
### Title: Filter for Expressed Genes
### Aliases: Expressed

### ** Examples

data(PhyloExpressionSetExample)

# remove genes that have an expression level below 8000 
# in at least one developmental stage
FilterConst <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
                         cut.off       = 8000, 
                         method        = "const",
                         comparison    = "below")
                              
dim(FilterConst) # check number of retained genes

# remove genes that have an expression level below 8000 
# in at least 3 developmental stages 
# (in this case: ceiling(7/2) = 4 stages fulfilling the cut-off criteria)
FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
                          cut.off       = 8000, 
                          method        = "min-set",
                          comparison    = "below")
                               
dim(FilterMinSet) # check number of retained genes

# remove genes that have an expression level below 8000 
# in at least 5 developmental stages (in this case: n = 2 stages fulfilling the criteria)
FilterNSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
                        cut.off       = 8000, 
                        method        = "n-set",
                        comparison    = "below",
                        n             = 2)
                               
dim(FilterMinSet) # check number of retained genes



# remove expression levels that exceed the cut.off criteria
FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
                          cut.off       = 12000, 
                          method        = "min-set",
                          comparison    = "above")
                               
dim(FilterMinSet) # check number of retained genes


# remove expression levels that undercut AND exceed the cut.off criteria
FilterMinSet <- Expressed(ExpressionSet = PhyloExpressionSetExample, 
                          cut.off       = c(8000,12000), 
                          method        = "min-set",
                          comparison    = "both")
                               
dim(FilterMinSet) # check number of retained genes





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Expressed", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("FlatLineTest")
### * FlatLineTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: FlatLineTest
### Title: Perform Flat Line Test
### Aliases: FlatLineTest

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)

# example PhyloExpressionSet using 100 permutations
FlatLineTest(PhyloExpressionSetExample,
             permutations  = 100,
             plotHistogram = FALSE)

# use your own permutation matrix based on which p-values (FlatLineTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

FlatLineTest(PhyloExpressionSetExample,
             custom.perm.matrix = custom_perm_matrix)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("FlatLineTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("GroupDiffs")
### * GroupDiffs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: GroupDiffs
### Title: Quantify the significant differences between gene expression
###   distributions of age groups
### Aliases: GroupDiffs

### ** Examples


data(PhyloExpressionSetExample)
# perform a Wilcoxon Rank Sum test to statistically quantify the
# difference between PS-Group 1 expression levels versus PS-Group 2
# expression levels
GroupDiffs(ExpressionSet = PhyloExpressionSetExample,
           Groups       = list(group_1 = 1:3,group_2 = 4:12),
           legendName   = "PS")

# quantify the significant difference of a selected set of genes
set.seed(123)
ExampleGeneSet <- sample(PhyloExpressionSetExample[ , 2],5000)  
             
GroupDiffs(ExpressionSet = PhyloExpressionSetExample,
           Groups       = list(group_1 = 1:3,group_2 = 4:12),
           legendName   = "PS",
           gene.set     = ExampleGeneSet)               





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("GroupDiffs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("LateConservationTest")
### * LateConservationTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: LateConservationTest
### Title: Perform Reductive Late Conservation Test
### Aliases: LateConservationTest

### ** Examples


data(PhyloExpressionSetExample)

# perform the late conservation test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
# stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
# module 3 = late
LateConservationTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7), 
                       permutations = 1000)


# use your own permutation matrix based on which p-values (LateConservationTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

LateConservationTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7), 
                       custom.perm.matrix = custom_perm_matrix)
                       



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("LateConservationTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("MatchMap")
### * MatchMap

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: MatchMap
### Title: Match a Phylostratigraphic Map or Divergence Map with a
###   ExpressionMatrix
### Aliases: MatchMap

### ** Examples

        
# load a standard PhyloExpressionSet
data(PhyloExpressionSetExample)
        
# in a standard PhyloExpressionSet, 
# column one and column two denote a standard 
# phylostratigraphic map
PhyloMap <- PhyloExpressionSetExample[ , 1:2]
        
# look at the phylostratigraphic map standard
head(PhyloMap)
        
# in a standard PhyloExpressionSet, column two combined 
# with column 3 - N denote a standard ExpressionMatrix
ExpressionMatrixExample <- PhyloExpressionSetExample[ , c(2,3:9)]
        
# these two data sets shall illustrate an example 
# phylostratigraphic map that is returned
# by a standard phylostratigraphy run, and a expression set 
# that is the result of expression data analysis 
# (background correction, normalization, ...)
        
# now we can use the MatchMap function to merge both data sets
# to obtain a standard PhyloExpressionSet
        
PES <- MatchMap(PhyloMap, ExpressionMatrixExample)
        
# note that the function returns a head() 
# of the matched gene ids to enable
# the user to find potential mis-matches
        
# the entire procedure is analogous to merge() 
# with two data sets sharing the same gene ids 
# as column (primary key)
PES_merge <- merge(PhyloMap, ExpressionMatrixExample)
        
        
        



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("MatchMap", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PairwiseTest")
### * PairwiseTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PairwiseTest
### Title: Perform Pairwise Difference Test
### Aliases: PairwiseTest

### ** Examples


data(PhyloExpressionSetExample)

# perform the pairwise difference test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to contrast1,
# stages 3-7 correspond to contrast 2.
# We test whether TAI in contrast1 is greater than contrast 2.
PairwiseTest(PhyloExpressionSetExample,
         modules = list(contrast1 = 1:2, contrast2 = 3:7), 
         altHypothesis = "greater",
         permutations = 1000)

# We can also test whether TAI in contrast1 is less than contrast 2.
PairwiseTest(PhyloExpressionSetExample,
         modules = list(contrast1 = 1:2, contrast2 = 3:7), 
         altHypothesis = "less",
         permutations = 1000)

# if we only want to test whether TAI in stage 1 (contrast 1) is greater than stage 3 (contrast 2).
PairwiseTest(PhyloExpressionSetExample,
         modules = list(contrast1 = 1, contrast2 = 2), 
         altHypothesis = "greater",
         permutations = 1000)

# use your own permutation matrix based on which p-values (PairwiseTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)
# We test whether TAI in contrast1 is greater than contrast 2.
PairwiseTest(PhyloExpressionSetExample,
         modules = list(contrast1 = 1:2, contrast2 = 3:7), 
         altHypothesis = "greater",
         custom.perm.matrix = custom_perm_matrix)
                       



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PairwiseTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotBarRE")
### * PlotBarRE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotBarRE
### Title: Plot Mean Relative Expression Levels as Barplot
### Aliases: PlotBarRE

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
PlotBarRE(ExpressionSet = PhyloExpressionSetExample,
          Groups        = list(c(1:3), c(4:12)))


# example DivergenceExpressionSet
PlotBarRE(ExpressionSet = DivergenceExpressionSetExample,
          Groups        = list(c(1:5), c(6:10)))


# Perform PlotBarRE() with p-value adjustment method Benjamini & Hochberg (1995)
PlotBarRE(ExpressionSet   = PhyloExpressionSetExample,
          Groups          = list(c(1:3), c(4:12)),
          p.adjust.method = "BH")
       
             
# Example: plot ratio
# the ratio curve visualizes the ratio between bar 1 / bar 2
# the z - axis shows the corresponding ratio value of bar 1 / bar 2
PlotBarRE(ExpressionSet = PhyloExpressionSetExample,
          Groups        = list(c(1:3), c(4:12)),
          ratio         = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotBarRE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCIRatio")
### * PlotCIRatio

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotCIRatio
### Title: Plot Transcriptome Index using bootstrapping and confidence
###   intervals
### Aliases: PlotCIRatio

### ** Examples

data("PhyloExpressionSetExample")
PlotCIRatio(PhyloExpressionSetExample,"TAI",5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCIRatio", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCategoryExpr")
### * PlotCategoryExpr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotCategoryExpr
### Title: Plot the Expression Levels of each Age or Divergence Category as
###   Boxplot, Violinplot, or Dotplot
### Aliases: PlotCategoryExpr

### ** Examples


data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

## Not run: 
##D 
##D # category-centered visualization of PS specific expression level distributions (log-scale)
##D PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
##D                      legendName    = "PS",
##D                      test.stat     = TRUE,
##D                      type          = "category-centered",
##D                      distr.type    = "boxplot",
##D                      log.expr      = TRUE)
##D                      
##D 
##D # stage-centered visualization of PS specific expression level distributions (log-scale)
##D PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
##D                      legendName    = "PS",
##D                      test.stat     = TRUE,
##D                      distr.type    = "boxplot",
##D                      type          = "stage-centered",
##D                      log.expr      = TRUE)
##D 
##D                      
##D                                                                
##D # category-centered visualization of PS specific expression level distributions (log-scale)
##D # as violoin plot
##D PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
##D                      legendName    = "PS",
##D                      test.stat     = TRUE,
##D                      distr.type    = "violin",
##D                      type          = "stage-centered",
##D                      log.expr      = TRUE)
##D 
##D 
##D 
##D 
##D # analogous for DivergenceExpressionSets
##D PlotCategoryExpr(ExpressionSet = DivergenceExpressionSetExample,
##D                      legendName    = "DS",
##D                      test.stat     = TRUE,
##D                      type          = "category-centered",
##D                      distr.type    = "boxplot",
##D                      log.expr      = TRUE)
##D 
##D 
##D # visualize the expression levels of 500 example genes
##D set.seed(234)
##D example.gene.set <- PhyloExpressionSetExample[sample(1:25260,500) , 2]
##D 
##D PlotCategoryExpr(ExpressionSet = PhyloExpressionSetExample,
##D                  legendName    = "PS",
##D                  test.stat     = TRUE,
##D                  type          = "category-centered",
##D                  distr.type    = "boxplot",
##D                  log.expr      = TRUE,
##D                  gene.set      = example.gene.set)
##D                  
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCategoryExpr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotContribution")
### * PlotContribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotContribution
### Title: Plot Cumulative Transcriptome Index
### Aliases: PlotContribution

### ** Examples


 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)
 
 # visualize phylostratum contribution to global TAI
 PlotContribution(PhyloExpressionSetExample, legendName = "PS")
 
 # visualize divergence stratum contribution to global TDI
 PlotContribution(DivergenceExpressionSetExample, legendName = "DS")
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotContribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotCorrelation")
### * PlotCorrelation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotCorrelation
### Title: Plot the Correlation Between Phylostrata and Divergence Strata
### Aliases: PlotCorrelation

### ** Examples

# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)
 
# plot the PS and DS correlation
PlotCorrelation(PhyloExpressionSetExample, 
                DivergenceExpressionSetExample, 
                method      = "pearson", 
                linearModel = TRUE)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotCorrelation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotDistribution")
### * PlotDistribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotDistribution
### Title: Plot the Frequency Distribution of Phylostrata or Divergence
###   Strata
### Aliases: PlotDistribution

### ** Examples


# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# plot the phylostratum distribution of a PhyloExpressionSet
PlotDistribution(PhyloExpressionSetExample)

# plot the relative frequency distribution of a PhyloExpressionSet
PlotDistribution(PhyloExpressionSetExample, as.ratio = TRUE)


# a example for visualizing the PS distribution for a subset of genes
PlotDistribution(PhyloExpressionSetExample[sample(20000,5000) , ], as.ratio = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotDistribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotEnrichment")
### * PlotEnrichment

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotEnrichment
### Title: Plot the Phylostratum or Divergence Stratum Enrichment of a
###   given Gene Set
### Aliases: PlotEnrichment

### ** Examples


data(PhyloExpressionSetExample)

set.seed(123)
test_set <- sample(PhyloExpressionSetExample[ , 2],10000)

## Examples with complete.bg = TRUE
## Hence: the entire background set of the input ExpressionSet is considered 
## when performing Fisher's exact test 

# measure: log-foldchange
PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
               test.set      = test_set , 
               legendName    = "PS", 
               measure       = "log-foldchange")
               
               
# measure: foldchange
PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
               test.set      = test_set , 
               legendName    = "PS", 
               measure       = "foldchange")
   
               
## Examples with complete.bg = FALSE
## Hence: the test.set genes are excluded from the background set before
## Fisher's exact test is performed
     
                                       
# measure: log-foldchange
PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
               test.set      = test_set ,
                complete.bg  = FALSE,
               legendName    = "PS", 
               measure       = "log-foldchange")
               
               
# measure: foldchange
PlotEnrichment(ExpressionSet = PhyloExpressionSetExample,
               test.set      = test_set , 
               complete.bg   = FALSE,
               legendName    = "PS", 
               measure       = "foldchange")     
               



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotEnrichment", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotGeneSet")
### * PlotGeneSet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotGeneSet
### Title: Plot the Expression Profiles of a Gene Set
### Aliases: PlotGeneSet

### ** Examples

data(PhyloExpressionSetExample)

# the best parameter setting to visualize this plot:
# png("test_png.png",700,400)
PlotGeneSet(ExpressionSet = PhyloExpressionSetExample, 
            gene.set      = PhyloExpressionSetExample[1:5, 2], 
            type          = "l", 
            lty           = 1, 
            lwd           = 4,
            xlab          = "Ontogeny",
            ylab          = "Expression Level")

# dev.off()

# In case you would like to work with the expression levels
# of selected genes you can specify the 'get.subset' argument:

PlotGeneSet(ExpressionSet = PhyloExpressionSetExample, 
            gene.set      = PhyloExpressionSetExample[1:5, 2], 
            get.subset    = TRUE)


# get a gene subset using only a phylostratihraphic map
ExamplePSMap <- PhyloExpressionSetExample[ , 1:2]

PlotGeneSet(ExpressionSet = ExamplePSMap, 
            gene.set      = PhyloExpressionSetExample[1:5, 2], 
            get.subset    = TRUE,
            use.only.map  = TRUE)
            



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotGeneSet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotGroupDiffs")
### * PlotGroupDiffs

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotGroupDiffs
### Title: Plot the significant differences between gene expression
###   distributions of PS or DS groups
### Aliases: PlotGroupDiffs

### ** Examples


data(PhyloExpressionSetExample)

PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
               Groups        = list(group_1 = 1:3,group_2 = 4:12),
               legendName    = "PS",
               type          = "b",
               lwd           = 6,
               xlab          = "Ontogeny")
               
               
# only receive the p-values without the corresponding plot               
PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
               Groups        = list(group_1 = 1:3,group_2 = 4:12),
               legendName    = "PS",
               plot.p.vals   = FALSE,
               type          = "b",
               lwd           = 6,
               xlab          = "Ontogeny")
               
               
# quantify the significant difference of a selected set of genes
# only receive the p-values without the corresponding plot
set.seed(123)
ExampleGeneSet <- sample(PhyloExpressionSetExample[ , 2],5000)
                              
PlotGroupDiffs(ExpressionSet = PhyloExpressionSetExample,
               Groups        = list(group_1 = 1:3,group_2 = 4:12),
               legendName    = "PS",
               plot.p.vals   = FALSE,
               gene.set      = ExampleGeneSet)                 


# plot differences as boxplot for each developmental stage
PlotGroupDiffs(ExpressionSet = tf(PhyloExpressionSetExample,log2),
               Groups        = list(group_1 = 1:3,group_2 = 4:12),
               legendName    = "PS",
               plot.type     = "boxplot")





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotGroupDiffs", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotMeans")
### * PlotMeans

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotMeans
### Title: Plot Mean Expression Profiles
### Aliases: PlotMeans

### ** Examples

### Example using a PhyloExpressionSet
### and DivergenceExpressionSet
# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# load PhyloExpressionSet
data(DivergenceExpressionSetExample)

# plot evolutionary old PS (PS1-3) vs evolutionary young PS (PS4-12)
PlotMeans(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = TRUE)

# if users wish to not adjust the y-axis scale when 
# 2 groups are selected they can specify: adjust.range = FALSE
PlotMeans(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = FALSE)
          
          
# plot conserved DS (DS1-5) vs divergent DS (PS6-10)
# NOTE: DS are always defined in the range 1, 2, ... , 10.
# Hence, make sure that your groups are within this range!
PlotMeans(DivergenceExpressionSetExample,
          Groups = list(c(1:5), c(6:10)), 
          legendName = "DS",
          adjust.range = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotMeans", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotMedians")
### * PlotMedians

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotMedians
### Title: Plot Median Expression Profiles
### Aliases: PlotMedians

### ** Examples

### Example using a PhyloExpressionSet
### and DivergenceExpressionSet
# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# load PhyloExpressionSet
data(DivergenceExpressionSetExample)

# plot evolutionary old PS (PS1-3) vs evolutionary young PS (PS4-12)
PlotMedians(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = TRUE)

# if users wish to not adjust the y-axis scale when 
# 2 groups are selected they can specify: adjust.range = FALSE
PlotMedians(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = FALSE)
          
          
# plot conserved DS (DS1-5) vs divergent DS (PS6-10)
# NOTE: DS are always defined in the range 1, 2, ... , 10.
# Hence, make sure that your groups are within this range!
PlotMedians(DivergenceExpressionSetExample,
          Groups = list(c(1:5), c(6:10)), 
          legendName = "DS",
          adjust.range = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotMedians", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotPattern")
### * PlotPattern

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotPattern
### Title: Plot the Transcriptome Age Index or Transcriptome Divergence
###   Index
### Aliases: PlotPattern

### ** Examples


# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# only visualize the TAI profile without any test statistics...
# this is equavalent to performing: plot(TAI(PhyloExpressionSetExample), type = "l", lwd = 6)
PlotPattern(ExpressionSet = PhyloExpressionSetExample,
            TestStatistic = NULL,
            type          = "l",
            xlab          = "Ontogeny",
            ylab          = "TAI",
            lwd           = 9)

# the simplest example of plotting the TAI profile of a given PhyloExpressionSet:
# In this case (default) the FlatLineTest will be performed to quantify
# the statistical significance of the present TAI pattern and will be drawn as 'p = ... '
# in the plot

PlotPattern(ExpressionSet = PhyloExpressionSetExample, 
            TestStatistic = "FlatLineTest",
            permutations  = 100, 
            type          = "l", 
            xlab          = "Ontogeny", 
            ylab          = "TAI", 
            lwd           = 9)

# an example performing the ReductiveHourglassTest and printing the p-value
# and shaded area of the presumptive phylotypic period into the plot
# Here the 'p = ...' denotes the p-value that is returned by the ReductiveHourglassTest

PlotPattern(
            ExpressionSet = PhyloExpressionSetExample,
            TestStatistic = "ReductiveHourglassTest",
            modules       = list(early = 1:2,mid = 3:5,late = 6:7), 
            permutations  = 100, 
            p.value       = TRUE, 
            shaded.area   = TRUE, 
            xlab          = "Ontogeny", 
            ylab          = "TAI", 
            type          = "l", 
            lwd           = 9)

# testing for early conservation model 
PlotPattern( ExpressionSet = PhyloExpressionSetExample,
             TestStatistic = "EarlyConservationTest",
            modules        = list(early = 1:2,mid = 3:5,late = 6:7), 
            permutations   = 100,
            p.value        = TRUE, 
            shaded.area    = TRUE, 
            xlab           = "Ontogeny", 
            ylab           = "TAI", 
            type           = "l", 
            lwd            = 9)
            

# use your own permutation matrix
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

PlotPattern(ExpressionSet      = PhyloExpressionSetExample, 
            TestStatistic      = "FlatLineTest",
            custom.perm.matrix = custom_perm_matrix, 
            type               = "l", 
            xlab               = "Ontogeny", 
            ylab               = "TAI", 
            lwd                = 9)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotPattern", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotRE")
### * PlotRE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotRE
### Title: Plot Relative Expression Levels
### Aliases: PlotRE

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
PlotRE(PhyloExpressionSetExample,
       Groups = list(c(1:3), c(4:12)), 
       legendName = "PS")


# or you can choose any combination of groups
PlotRE(PhyloExpressionSetExample,
       Groups = list(c(1,7,9), c(2:6,8,10:12)),
       legendName = "PS")

    
# example DivergenceExpressionSet
PlotRE(DivergenceExpressionSetExample,
       Groups = list(c(1:5), c(6:10)), 
       legendName = "DS")


  



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotRE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotReplicateQuality")
### * PlotReplicateQuality

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotReplicateQuality
### Title: Plot the Quality of Biological Replicates
### Aliases: PlotReplicateQuality

### ** Examples


data(PhyloExpressionSetExample)

# visualize log(var(x)) between replicates for each gene and developmental stage 
PlotReplicateQuality(ExpressionSet = PhyloExpressionSetExample[1:5000 , 1:8],
                     nrep          = 2,
                     legend.pos   = "topright",
                     ylim          = c(0,0.2),
                     lwd           = 6)
                     




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotReplicateQuality", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSelectedAgeDistr")
### * PlotSelectedAgeDistr

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSelectedAgeDistr
### Title: Plot the PS or DS distribution of a selected set of genes
### Aliases: PlotSelectedAgeDistr

### ** Examples

data(PhyloExpressionSetExample)

# generate an example gene set
set.seed(123)
ExGeneSet <- sample(PhyloExpressionSetExample[ , 2], 5000)

# gene count example
PlotSelectedAgeDistr(ExpressionSet = PhyloExpressionSetExample,
                     gene.set      = ExGeneSet,
                     legendName    = "PS",
                     as.ratio      = TRUE)

# relative gene count example
PlotSelectedAgeDistr(ExpressionSet = PhyloExpressionSetExample,
                     gene.set      = ExGeneSet,
                     legendName    = "PS",
                     as.ratio      = FALSE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSelectedAgeDistr", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSignature")
### * PlotSignature

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSignature
### Title: Plot evolutionary signatures across transcriptomes
### Aliases: PlotSignature

### ** Examples

data(PhyloExpressionSetExample)

# plot TAI pattern and perform flat line test
PlotSignature(PhyloExpressionSetExample, 
              measure       = "TAI", 
              permutations  = 100,
              TestStatistic = "FlatLineTest",
              ylab = "Transcriptome Age Index")
              



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSignature", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotSignatureTransformed")
### * PlotSignatureTransformed

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotSignatureTransformed
### Title: Plot evolutionary signatures across transcriptomes and RNA-seq
###   transformations
### Aliases: PlotSignatureTransformed

### ** Examples

## Not run: 
##D data(PhyloExpressionSetExample)
##D 
##D # Flat line test
##D PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
##D                     TestStatistic = "FlatLineTest",
##D                     transforms = c("none", "log2", "sqrt", "rank", "squared"))
##D 
##D # Reductive hourglass test
##D PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
##D                      TestStatistic = "ReductiveHourglassTest",
##D                      transforms = c("none", "log2", "sqrt", "rank", "squared"),
##D                      modules = list(early = 1:2, mid = 3:5, late = 6:7))
##D 
##D library(DESeq2)
##D PlotSignatureTransformed(ExpressionSet = PhyloExpressionSetExample,
##D                      TestStatistic = "ReductiveHourglassTest",
##D                      transforms = c("none", "log2", "sqrt", "vst", "rank", "squared"),
##D                      modules = list(early = 1:2, mid = 3:5, late = 6:7))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotSignatureTransformed", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("PlotVars")
### * PlotVars

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: PlotVars
### Title: Plot Variance of Expression Profiles
### Aliases: PlotVars

### ** Examples

### Example using a PhyloExpressionSet
### and DivergenceExpressionSet
# load PhyloExpressionSet
data(PhyloExpressionSetExample)

# load PhyloExpressionSet
data(DivergenceExpressionSetExample)

# plot evolutionary old PS (PS1-3) vs evolutionary young PS (PS4-12)
PlotVars(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = TRUE)

# if users wish to not adjust the y-axis scale when 
# 2 groups are selected they can specify: adjust.range = FALSE
PlotVars(PhyloExpressionSetExample,
          Groups = list(c(1:3), c(4:12)), 
          legendName = "PS",
          adjust.range = FALSE)
          
          
# plot conserved DS (DS1-5) vs divergent DS (PS6-10)
# NOTE: DS are always defined in the range 1, 2, ... , 10.
# Hence, make sure that your groups are within this range!
PlotVars(DivergenceExpressionSetExample,
          Groups = list(c(1:5), c(6:10)), 
          legendName = "DS",
          adjust.range = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("PlotVars", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("RE")
### * RE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: RE
### Title: Transform to Relative Expression Levels
### Aliases: RE

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)

# relative expression profile of PS1 genes
RE(PhyloExpressionSetExample[ which(PhyloExpressionSetExample[ , 1] == 1), 3:9 ])





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("RE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("REMatrix")
### * REMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: REMatrix
### Title: Compute a Relative Expression Matrix
### Aliases: REMatrix

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
REMatrix(PhyloExpressionSetExample)

# example DivergenceExpressionSet
REMatrix(DivergenceExpressionSetExample)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("REMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReductiveHourglassTest")
### * ReductiveHourglassTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReductiveHourglassTest
### Title: Perform the Reductive Hourglass Test
### Aliases: ReductiveHourglassTest

### ** Examples


data(PhyloExpressionSetExample)

# perform the reductive hourglass test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
# stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
# module 3 = late
ReductiveHourglassTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7), 
                       permutations = 1000)


# use your own permutation matrix based on which p-values (ReductiveHourglassTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

ReductiveHourglassTest(PhyloExpressionSetExample,
                     modules = list(early = 1:2, mid = 3:5, late = 6:7),
                     custom.perm.matrix = custom_perm_matrix)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReductiveHourglassTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ReverseHourglassTest")
### * ReverseHourglassTest

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ReverseHourglassTest
### Title: Perform the Reverse Hourglass Test
### Aliases: ReverseHourglassTest

### ** Examples


data(PhyloExpressionSetExample)

# perform the reductive hourglass test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
# stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
# module 3 = late
ReverseHourglassTest(PhyloExpressionSetExample,
                       modules = list(early = 1:2, mid = 3:5, late = 6:7),
                       permutations = 1000)


# use your own permutation matrix based on which p-values (ReverseHourglassTest)
# shall be computed
custom_perm_matrix <- bootMatrix(PhyloExpressionSetExample,100)

ReverseHourglassTest(PhyloExpressionSetExample,
                     modules = list(early = 1:2, mid = 3:5, late = 6:7),
                     custom.perm.matrix = custom_perm_matrix)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ReverseHourglassTest", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("SelectGeneSet")
### * SelectGeneSet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: SelectGeneSet
### Title: Select a Subset of Genes in an ExpressionSet
### Aliases: SelectGeneSet

### ** Examples

data(PhyloExpressionSetExample)

# receive a subset ExpressionSet for the fist 5 genes stored in
# the PhyloExpressionSetExample
SelectGeneSet(ExpressionSet = PhyloExpressionSetExample,
            gene.set      = PhyloExpressionSetExample[1:5, 2])
            
            
# get a gene subset using only a phylostratihraphic map
ExamplePSMap <- PhyloExpressionSetExample[ , 1:2]

SelectGeneSet(ExpressionSet = ExamplePSMap,
              gene.set      = PhyloExpressionSetExample[1:5, 2],
              use.only.map  = TRUE)          
            



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("SelectGeneSet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("TAI")
### * TAI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: TAI
### Title: Compute the Transcriptome Age Index (TAI)
### Aliases: TAI

### ** Examples


# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample)

# computing the TAI profile of a given PhyloExpressionSet object
TAIs <- TAI(PhyloExpressionSetExample)






base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("TAI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("TDI")
### * TDI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: TDI
### Title: Compute the Transcriptome Divergence Index (TDI)
### Aliases: TDI

### ** Examples


# reading a standard DivergenceExpressionSet
data(DivergenceExpressionSetExample)

# computing the TDI profile of a given DivergenceExpressionSet object
TDIs <- TDI(DivergenceExpressionSetExample)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("TDI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("TEI")
### * TEI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: TEI
### Title: Compute the Transcriptome Evolutionary Index (TEI)
### Aliases: TEI

### ** Examples


# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample, package = "myTAI")

# computing the TEI profile of a given PhyloExpressionSet object
TEI <- TEI(PhyloExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("TEI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("TPI")
### * TPI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: TPI
### Title: Compute the Transcriptome Polymorphism Index (TPI)
### Aliases: TPI

### ** Examples

## Not run: 
##D # reading a standard PolymorphismExpressionSet
##D data(PolymorphismExpressionSetExample)
##D 
##D # computing the TPI profile of a given PolymorphismExpressionSet object
##D TPIs <- TPI(PolymorphismExpressionSet)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("TPI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("age.apply")
### * age.apply

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: age.apply
### Title: Age Category Specific apply Function
### Aliases: age.apply

### ** Examples

 
 # source the example dataset
 data(PhyloExpressionSetExample)
 
# Example 1
# get the relative expression profiles for each phylostratum
age.apply(PhyloExpressionSetExample, RE)

# this is analogous to 
REMatrix(PhyloExpressionSetExample)
# Example 2
# compute the mean expression profiles for each phylostratum
age.apply(PhyloExpressionSetExample, colMeans)

# Example 3
# compute the variance profiles for each phylostratum
age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , var))

# Example 4
# compute the range for each phylostratum
# Note: in this case, the range() function returns 2 values for each phylostratum
# and each developmental stage, hence one should use the argument 'as.list = TRUE'
# to make sure that the results are returned properly 
age.apply(PhyloExpressionSetExample, function(x) apply(x , 2 , range), as.list = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("age.apply", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bar.colors")
### * bar.colors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bar.colors
### Title: Color palette for barplots
### Aliases: bar.colors
### Keywords: internal

### ** Examples


# get 5 different colors for 5 different bars
barplot_colors <- bar.colors(5)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bar.colors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bootMatrix")
### * bootMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootMatrix
### Title: Compute a Permutation Matrix for Test Statistics
### Aliases: bootMatrix

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet using 100 permutations
randomTAI.Matrix <- bootMatrix(PhyloExpressionSetExample, permutations = 100)

# example DivergenceExpressionSet using 100 permutations
randomTDI.Matrix <- bootMatrix(DivergenceExpressionSetExample, permutations = 100)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("bootTEI")
### * bootTEI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: bootTEI
### Title: Compute a Permutation Matrix of Transcriptome Evolutionary Index
###   (TEI)
### Aliases: bootTEI

### ** Examples


# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample, package = "myTAI")

# computing partial TEI contribution per gene
bM <- bootTEI(PhyloExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("bootTEI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("ecScore")
### * ecScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: ecScore
### Title: Compute the Hourglass Score for the EarlyConservationTest
### Aliases: ecScore

### ** Examples


 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # Example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the early conservation score for the TAI profile
 ec_score <- ecScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7)


 # Example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the early conservation score for the TDI profile
 ec_score <- ecScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7)
 
 # compute ecScore() vector from bootMatrix()
 apply(bootMatrix(PhyloExpressionSetExample,10),1,ecScore,early = 1:2,mid = 3:5,late = 6:7)
 
 # get warning if the expected pattern isn't followed
 ec_score <- ecScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,profile.warn=TRUE)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("ecScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("geom.mean")
### * geom.mean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: geom.mean
### Title: Geometric Mean
### Aliases: geom.mean

### ** Examples

x <- 1:10

geom.mean(x)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("geom.mean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("harm.mean")
### * harm.mean

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: harm.mean
### Title: Harmonic Mean
### Aliases: harm.mean

### ** Examples

x <- 1:10

harm.mean(x)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("harm.mean", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("is.ExpressionSet")
### * is.ExpressionSet

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: is.ExpressionSet
### Title: Test ExpressionSet Standard
### Aliases: is.ExpressionSet
### Keywords: internal

### ** Examples


# read example PhyloExpressionSet
data(PhyloExpressionSetExample)

is.ExpressionSet(PhyloExpressionSetExample)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("is.ExpressionSet", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lcScore")
### * lcScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lcScore
### Title: Compute the Hourglass Score for the LateConservationTest
### Aliases: lcScore

### ** Examples


 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # Example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the late conservation score for the TAI profile
 lc_score <- lcScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7)


 # Example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the late conservation score for the TDI profile
 lc_score <- lcScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7)
 
 # compute lcScore() vector from bootMatrix()
 apply(bootMatrix(PhyloExpressionSetExample,10),1,lcScore,early = 1:2,mid = 3:5,late = 6:7)
 
 # get warning if the expected pattern isn't followed
 lc_score <- lcScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,profile.warn=TRUE)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lcScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("omitMatrix")
### * omitMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: omitMatrix
### Title: Compute TAI or TDI Profiles Omitting a Given Gene
### Aliases: omitMatrix

### ** Examples


# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
omMatrix_ps <- omitMatrix(PhyloExpressionSetExample)

# example DivergenceExpressionSet
omMatrix_ds <- omitMatrix(DivergenceExpressionSetExample)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("omitMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pMatrix")
### * pMatrix

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pMatrix
### Title: Compute Partial TAI or TDI Values
### Aliases: pMatrix

### ** Examples



# read standard phylotranscriptomics data
data(PhyloExpressionSetExample)
data(DivergenceExpressionSetExample)

# example PhyloExpressionSet
PTM_ps <- pMatrix(PhyloExpressionSetExample)

# example DivergenceExpressionSet
PTM_ds <- pMatrix(DivergenceExpressionSetExample)

# boxplot of the pMatrix
boxplot(pMatrix(PhyloExpressionSetExample),outline = FALSE)

# boxplot of the pMatrix using log2 transformed expression levels
boxplot(pMatrix(tf(PhyloExpressionSetExample,log2)))
 




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pMatrix", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pMatrixTEI")
### * pMatrixTEI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pMatrixTEI
### Title: Compute Partial Transcriptome Evolutionary Index (TEI) Values
### Aliases: pMatrixTEI

### ** Examples


# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample, package = "myTAI")

# computing partial TEI contribution per gene
pMT <- pMatrixTEI(PhyloExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pMatrixTEI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pStrata")
### * pStrata

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pStrata
### Title: Compute Partial Strata Values
### Aliases: pStrata

### ** Examples


data(PhyloExpressionSetExample)

# compute partial TAI values for each Phylostratum
partialStrata <- pStrata(PhyloExpressionSetExample)

# show that colSums of pStrata is equavalent to the TAI values
all.equal(colSums(partialStrata),TAI(PhyloExpressionSetExample))

# show that colSums of pStrata is equavalent to colSums of pMatrix(PhyloExpressionSetExample)
all.equal(colSums(partialStrata),colSums(pMatrix(PhyloExpressionSetExample)))





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pStrata", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pStrataTEI")
### * pStrataTEI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pStrataTEI
### Title: Compute Partial Transcriptome Evolutionary Index (TEI) Strata
###   Values
### Aliases: pStrataTEI

### ** Examples


# reading a standard PhyloExpressionSet
data(PhyloExpressionSetExample, package = "myTAI")

# computing partial TEI contribution per gene
pS <- pStrataTEI(PhyloExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pStrataTEI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pTAI")
### * pTAI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pTAI
### Title: Compute the Phylostratum Contribution to the Global
###   Transcriptome Age Index
### Aliases: pTAI

### ** Examples


data(PhyloExpressionSetExample)

# get the partial contribution of phylostrata to the global
# TAI pattern
pTAI(PhyloExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pTAI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pTDI")
### * pTDI

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pTDI
### Title: Compute the Divergence Stratum Contribution to the Global
###   Transcriptome Divergence Index
### Aliases: pTDI

### ** Examples


data(DivergenceExpressionSetExample)

# get the partial contribution of divergence strata to the global
# TDI pattern
pTAI(DivergenceExpressionSetExample)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pTDI", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("pairScore")
### * pairScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: pairScore
### Title: Compute the Pairwise Difference in TAI Score
### Aliases: pairScore

### ** Examples


 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # Example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the pair score for the first two stages in the TAI profile
 # we test whether TAI in contrast1 is greater than contrast 2.
 pair_score <- pairScore(age_vals = TAIs,contrast1 = 1,contrast2 = 2,
                         altHypothesis="greater")


 # Example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the pair score for the first two stages in the TDI profile
 # we test whether TDI in contrast1 is greater than contrast 2.
 pair_score <- pairScore(age_vals = TDIs,contrast1 = 1,contrast2 = 2,
                         altHypothesis="greater")
 
 # compute pairScore() vector from bootMatrix()
 apply(bootMatrix(PhyloExpressionSetExample,10),1,
       pairScore,contrast1 = 1,contrast2 = 2, altHypothesis="greater")
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("pairScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rcpp_boottei_parallel")
### * rcpp_boottei_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rcpp_boottei_parallel
### Title: rcpp_boottei_parallel
### Aliases: rcpp_boottei_parallel

### ** Examples

## load example PhyloExpressionSetExample

data("PhyloExpressionSetExample", package="myTAI")

## convert into sparseMatrix - rownames GeneID

spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
    "sparseMatrix")
rownames(spmat) <- PhyloExpressionSetExample$GeneID

## create named Phylostratum vector

ps <- setNames(PhyloExpressionSetExample$Phylostratum,
    PhyloExpressionSetExample$GeneID)

## get permutations
rcpp_boottei_parallel(spmat, ps, 100, 1)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rcpp_boottei_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rcpp_pMatrix_parallel")
### * rcpp_pMatrix_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rcpp_pMatrix_parallel
### Title: rcpp_pMatrix_parallel
### Aliases: rcpp_pMatrix_parallel

### ** Examples

## load example PhyloExpressionSetExample

data("PhyloExpressionSetExample", package="myTAI")

## convert into sparseMatrix - rownames GeneID

spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
    "sparseMatrix")
rownames(spmat) <- PhyloExpressionSetExample$GeneID

## create named Phylostratum vector

ps <- setNames(PhyloExpressionSetExample$Phylostratum,
    PhyloExpressionSetExample$GeneID)

## get pMatrix
rcpp_pMatrix_parallel(spmat, ps)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rcpp_pMatrix_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rcpp_pStrata_parallel")
### * rcpp_pStrata_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rcpp_pStrata_parallel
### Title: rcpp_pStrata_parallel
### Aliases: rcpp_pStrata_parallel

### ** Examples

## load example PhyloExpressionSetExample

data("PhyloExpressionSetExample", package="myTAI")

## convert into sparseMatrix - rownames GeneID

spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]),
    "sparseMatrix")
rownames(spmat) <- PhyloExpressionSetExample$GeneID

## create named Phylostratum vector

ps <- setNames(PhyloExpressionSetExample$Phylostratum,
    PhyloExpressionSetExample$GeneID)
psgroup <- sort(unique(ps))

## get pStrata
rcpp_pStrata_parallel(spmat, ps, psgroup)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rcpp_pStrata_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rcpp_tei_parallel")
### * rcpp_tei_parallel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rcpp_tei_parallel
### Title: rcpp_tei_parallel
### Aliases: rcpp_tei_parallel

### ** Examples

## load example sequence data
data("PhyloExpressionSetExample", package="myTAI")
spmat <- as(data.matrix(PhyloExpressionSetExample[,-c(1,2)]), "sparseMatrix")
rownames(spmat) <- PhyloExpressionSetExample$GeneID
ps <- setNames(PhyloExpressionSetExample$Phylostratum, PhyloExpressionSetExample$GeneID)
rcpp_tei_parallel(spmat, ps)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rcpp_tei_parallel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("reversehourglassScore")
### * reversehourglassScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: reversehourglassScore
### Title: Compute the Reverse Hourglass Score for the Reverse Hourglass
###   Test
### Aliases: reversehourglassScore

### ** Examples


 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the global reverse hourglass destruction score 
 # for the TAIs profile using reduction method: mean(mean-mean)
 reversehourglass_score <- reversehourglassScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean")


 # example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the global reverse hourglass destruction score for the TDIs profile 
 # using reduction method: mean(mean-mean)
 reversehourglass_score <- reversehourglassScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean")
                     
 # get warning if the expected pattern isn't followed
 reversehourglass_score <- reversehourglassScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean",profile.warn=TRUE)
 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("reversehourglassScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("rhScore")
### * rhScore

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: rhScore
### Title: Compute the Hourglass Score for the Reductive Hourglass Test
### Aliases: rhScore

### ** Examples


 # read standard phylotranscriptomics data
 data(PhyloExpressionSetExample)
 data(DivergenceExpressionSetExample)

 # example PhyloExpressionSet:

 # compute the TAI profile
 TAIs <- TAI(PhyloExpressionSetExample)

 # compute the global hourglass destruction score 
 # for the TAIs profile using reduction method: mean(mean-mean)
 rh_score <- rhScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean")


 # example DivergenceExpressionSet:

 # compute the TDI profile
 TDIs <- TDI(DivergenceExpressionSetExample)

 # compute the global hourglass destruction score for the TDIs profile 
 # using reduction method: mean(mean-mean)
 rh_score <- rhScore(age_vals = TDIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean")
                     
 # get warning if the expected pattern isn't followed
 rh_score <- rhScore(age_vals = TAIs,early = 1:2,mid = 3:5,late = 6:7,
                     method = "mean",scoringMethod = "mean-mean",profile.warn=TRUE) 



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("rhScore", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("taxid")
### * taxid

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: taxid
### Title: Retrieve taxonomy categories from NCBI Taxonomy
### Aliases: taxid

### ** Examples

## Not run: 
##D # download categories.dmp file to current working directory 
##D # and filter for 'Archea' taxids
##D Archea.taxids <- taxid(db.path = getwd(), filter = "Archea", download = TRUE)
##D 
##D # Once the NCBI Taxonomy 'categories.dmp' file is downloaded to your machine ('download = TRUE')
##D # the 'taxid()' function can be proceed on the local 'categories.dmp' file
##D # e.g. filter for Virus taxids
##D Virus.taxids <- taxid(db.path = getwd(), filter = "Viruses")
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("taxid", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tf")
### * tf

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tf
### Title: Transform Gene Expression Levels
### Aliases: tf

### ** Examples

## Not run: 
##D data(PhyloExpressionSetExample)
##D 
##D # a simple example is to transform the gene expression levels
##D # of a given PhyloExpressionSet using a sqrt or log2 transformation
##D 
##D PES.sqrt <- tf(PhyloExpressionSetExample, sqrt)
##D 
##D PES.log2 <- tf(PhyloExpressionSetExample, log2)
##D 
##D # plotting the TAI using log2 transformed expression levels
##D # and performing the Flat Line Test to obtain the p-value
##D PlotSignature(ExpressionSet = tf(PhyloExpressionSetExample, log2),
##D             permutations = 1000)
##D 
##D 
##D 
##D # in case the expression matrix contains 0s, a pseudocount can be added prior
##D # to certain transformations, e.g. log2(x+1) where 1 is the pseudocount.
##D 
##D PhyloExpressionSetExample[4,3] = 0
##D PES.log2 <- tf(PhyloExpressionSetExample, log2, pseudocount = 0)
##D 
##D # this should return -Inf at PES.log2[4,3] the issue here is that 
##D # -Inf cannot be used to compute the phylotranscriptomic profile.
##D 
##D PES.log2 <- tf(PhyloExpressionSetExample, log2, pseudocount = 1)
##D # log2 transformed expression levels can now be used in downstream analyses.
##D 
##D 
##D # to perform rank transformation
##D 
##D PES.rank <- tf(PhyloExpressionSetExample, FUN = function(x) apply(x, 2, base::rank))
##D 
##D 
##D # rlog and vst transformations are now also possible by loading the DESeq2 package
##D # and transforming the data with the parameter integerise = TRUE.
##D library(DESeq2) # make sure the DESeq2 version >= 1.29.15 for rlog
##D PES.vst <- tf(PhyloExpressionSetExample, vst, integerise = TRUE)
## End(Not run)





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tf", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tfPS")
### * tfPS

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tfPS
### Title: Transform Phylostratum Values
### Aliases: tfPS

### ** Examples

# source the example dataset
data(PhyloExpressionSetExample)
 
# get the relative expression profiles for each phylostratum
tfPES <- tfPS(PhyloExpressionSetExample, transform = "qr")
head(tfPES)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tfPS", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("tfStability")
### * tfStability

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: tfStability
### Title: Perform Permutation Tests Under Different Transformations
### Aliases: tfStability

### ** Examples


data(PhyloExpressionSetExample)

# perform the reductive hourglass test for a PhyloExpressionSet
# here the prior biological knowledge is that stages 1-2 correspond to module 1 = early,
# stages 3-5 to module 2 = mid (phylotypic module), and stages 6-7 correspond to
# module 3 = late
tfStability(ExpressionSet = PhyloExpressionSetExample,
                     TestStatistic = "ReductiveHourglassTest",
                     permutations       = 100,
                     transforms = c("log2", "sqrt", "none"),
                     modules = list(early = 1:2, mid = 3:5, late = 6:7))


# it is also possible to test the phylotranscriptomic pattern using rlog
# and vst transforms from DESeq2

library(DESeq2)
tfStability(ExpressionSet = PhyloExpressionSetExample,
                     TestStatistic = "ReductiveHourglassTest",
                     permutations       = 100,
                     transforms = c("log2", "sqrt", "none", "vst"),
                     modules = list(early = 1:2, mid = 3:5, late = 6:7))





base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("tfStability", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
