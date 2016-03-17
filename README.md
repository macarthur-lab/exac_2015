# ExAC analysis scripts

## Analysis and figures from ExAC flagship paper, bioRxiv 2015

This package can recreate the figures from the [main ExAC paper](http://biorxiv.org/content/early/2015/10/30/030338).

It requires the following packages:

```R
install.packages(c('binom', 'plyr', 'ggplot2', 'data.table', 'reshape', 'plotrix', 'dplyr', 'Hmisc', 'gdata', 'magrittr', 'vioplot'))
```

## Basic usage

Warning: this process requires a decent amount of memory, ideally on a machine that has at least 16G RAM (newer Macbook Pros are fine, but not Macbook Airs).
On Linux and machines with stricter memory policies, we recommend 32G or even 64G to allow for overhead in loading.

```R
source('exac_constants.R')
exac = load_exac_data()
```

exac is then the data frame of ALL ~10M variants in the dataset with over ~100 annotations per variant. Each variant is now its own entry, as opposed to a typical VCF where multi-allelic variants are combined into a single line. **Note that these are unfiltered variant calls**. Make sure to filter to either PASS sites, or better yet, our criteria for high-quality calls (given by the column `use`, see [this blog post](http://macarthurlab.org/2016/03/17/reproduce-all-the-figures-a-users-guide-to-exac-part-2) for details):

```R
filtered_calls = subset(exac, filter == 'PASS')
```
or:
```R
use_data = subset(exac, use)
```

The code for figures is in `exac_figures.R`. This process is very memory intensive, so depending on your system, it may be easier to open the script and run sections of the code independently.

## Additional analyses

_Note for ExAC analysts_: When writing additional analyses, for consistency, it is suggested to start your script with the following
(including exac, constraint, etc., as needed):

```R
source('exac_constants.R')
if (!("exac" %in% ls(globalenv()))) {
  exac = load_exac_data()
}
if (!("use_data" %in% ls(globalenv()))) {
  use_data = subset(exac, use)
}
if (!("constraint" %in% ls(globalenv()))) {
  constraint = load_constraint_data()
}
```

Then, use `exac` as the full exac data frame, and `use_data` as the high-confidence variant set.

Alternatively, particularly for scripts that will be used on arbitrary subsets of data,
wrap your scripts in functions that take an `exac` variable as an argument.

If you are using `source('script.R')`, assume the user is running the code from the root of `exac_papers`.