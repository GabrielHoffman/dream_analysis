## dream: Powerful differential expression analysis for repeated measures designs

<p align="left">
<img src="https://raw.githubusercontent.com/GabrielHoffman/gabrielhoffman.github.io/master/img/dream_icon.png" width="400">
</p>

`dream` is available in the Bioconductor package [variancePartition](http://bioconductor.org/packages/release/bioc/html/variancePartition.html) ⩾ v1.10.1 for Bioconductor ⩾ v3.7.  You may need to install [R](https://www.r-project.org) ⩾ v3.5.

See dream [documentation](http://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/dream.html).


### Simulations
See code in [sims/master_analysis.sh](https://github.com/GabrielHoffman/dream_analysis/tree/master/sims/master_analysis.sh)

### Data analysis

Data and annotation are stored on [Synapse](https://www.synapse.org).  In order to run these scripts, you will need the [Synapse client for R](https://docs.synapse.org/articles/getting_started.html) and a [Synapse account](https://www.synapse.org/#!RegisterAccount:0).

Links point to original publication, public data, html files from current analysis, and text files from differential expression analysis. `rmarkdown::render` code snippet will regenerate results and html files. 

- ##### Timothy Syndrome | [paper](https://www.nature.com/articles/nm.2576), [data](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25542), [html](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/master/results/iPSC_Timothy.html), [results](https://github.com/GabrielHoffman/dream_analysis/tree/master/results/files/iPSC_Timothy)

`rmarkdown::render("src/iPSC_Timothy.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Alzheimer's disease | [paper](https://www.nature.com/articles/sdata2018185), [data](https://www.synapse.org/#!Synapse:syn3159438), [html](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/master/results/AMP_AD.html), [results](https://github.com/GabrielHoffman/dream_analysis/tree/master/results/files/AMP_AD)

`rmarkdown::render("src/AMP_AD.Rmd", output_dir='./', intermediates_dir='./')`

- ##### Childhood onset schizophrenia | [paper](https://www.nature.com/articles/s41467-017-02330-5), [data](www.synapse.org/hiPSC_COS), [html](https://cdn.rawgit.com/GabrielHoffman/dream_analysis/master/results/COS.html), [results](https://github.com/GabrielHoffman/dream_analysis/tree/master/results/files/COS)

`rmarkdown::render("src/COS.Rmd", output_dir='./', intermediates_dir='./')`

