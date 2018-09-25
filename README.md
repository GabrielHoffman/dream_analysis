## dream: Powerful differential expression analysis for repeated measures designs

<p align="left">
<img src="https://raw.githubusercontent.com/GabrielHoffman/gabrielhoffman.github.io/master/img/dream_icon.png" width="400">
</p>

`dream` is available into the Bioconductor package [variancePartition](http://bioconductor.org/packages/release/bioc/html/variancePartition.html).  See [documentation](http://bioconductor.org/packages/devel/bioc/vignettes/variancePartition/inst/doc/dream.html)


### Simulations
See code in `sims/master_analysis.sh`

### Data analysis

Data and annotation are stored on [Synapse](https://www.synapse.org).  In order to run these scripts, you will need the [Synapse client for R](https://docs.synapse.org/articles/getting_started.html) and a [Synapse account](https://www.synapse.org/#!RegisterAccount:0)

- ##### Timothy Syndrome | [paper](https://www.nature.com/articles/nm.2576), [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25542), [results](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/7f756df5/results/iPSC_Timothy.html)

`rmarkdown::render("src/iPSC_Timothy.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Alzheimer's disease | [paper](https://www.nature.com/articles/sdata2018185), [data](https://www.synapse.org/#!Synapse:syn3159438), [results](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/7f756df5/results/AMP_AD.html)

`rmarkdown::render("src/AMP_AD.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Childhood onset schizophrenia | [paper](https://www.nature.com/articles/s41467-017-02330-5), [data](www.synapse.org/hiPSC_COS), [results](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/7f756df5/results/COS.html)

`rmarkdown::render("src/COS.Rmd", output_dir='./', intermediates_dir='./')`



https://github.com/GabrielHoffman/gabrielhoffman.github.io/blob/master/img/dream_icon_wide.png

https://raw.githubusercontent.com/GabrielHoffman/gabrielhoffman.github.io/blob/master/img/dream_icon_wide.png