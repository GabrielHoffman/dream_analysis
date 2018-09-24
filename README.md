## dream: Powerful differential expression analysis for repeated measures designs

`dream` is available into the Bioconductor package [variancePartition](http://bioconductor.org/packages/release/bioc/html/variancePartition.html).  See [documentation](http://bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/dream.html)


### Simulations
See code in `sims/master_analysis.sh`

### Data analysis

Data and annotation are stored on [Synapse](synapse.org).  In order to run these scripts, you will need the [Synapse client for R](https://docs.synapse.org/articles/getting_started.html) and a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)

- ##### Timothy Syndrome

`rmarkdown::render("src/iPSC_Timothy.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Alzheimer's disease

`rmarkdown::render("src/AMP_AD.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Childhood onset schizophrenia

`rmarkdown::render("src/COS.Rmd", output_dir='./', intermediates_dir='./')`





