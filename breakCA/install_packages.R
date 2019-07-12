install.packages(c('data.table', 'plyr', 'dplyr', 'pbapply', 'readr', 'reshape', 'rmutil', 'lattice', 'stringr', 'mlr', 'randomForest'), repos='https://cloud.r-project.org/')

source("https://bioconductor.org/biocLite.R");
biocLite(c("biovizBase", "rtracklayer", "Rsamtools", "BSgenome.Hsapiens.UCSC.hg19", "GenomicAlignments", "VariantAnnotation"))
