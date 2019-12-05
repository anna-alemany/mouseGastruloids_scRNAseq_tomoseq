# Tomo-seq notebooks

### Summary of the data

We start with the following datasets:

* 3 tomo-seq replicates of E14 gastruloids
* 5 tomo-seq replicates of Lfng gastruloids
* 3 tomo-seq replicates of E8.5 mouse embryos
* 3 replicates of a published microarray dataset where the posterior mesoderm (from the tail bud to the newly formed somite) of E9.5 mouse embryos was dissected and sectioned in 7 pieces (accession number GSE39615).

All tomo-seq RNA-seq datasets can be found in the Gene Expression Omnibus (GEO) under accession code [GSE123187](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187) and can be explored   [here](https://avolab.hubrecht.eu/MouseGastruloids2019)

### Organization of notebooks in this repository subfolder

For each dataset, a folder is presented that contains python notebooks used to filter and normalized the pre-processed sequencing (or microarray) data, correct batch effects (when necessary), and identify reproducible genes between replicates. For the notebooks to work, input datasets are required (most of them found [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123187); otherwise please contact me. 

Additionally, python notebooks containing comparisons between different systems (E14 gastruloids vs Lfng gastruloids; E14 gastruloids vs Lfng gastruloids vs E8.5 mouse embryos; etc) are provided. The main idea is always to select interesting genes to compare based on our reproducibility criteria and cluster the according to expression patterns along the anterior/posterior direction using self-organizing maps and hierarchical clustering. 

### Additional note

Even not explicitely stated, all python notebooks contain the following libraries by default:

```{python}
import numpy as np
import pandas as pd
from pandas.io.parsers import read_csv
from collections import Counter
import matplotlib.pyplot as plt
```
