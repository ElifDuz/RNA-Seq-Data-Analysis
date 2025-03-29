This file contains the data analysis of RNA-Seq datasets. 
These datasets can be downloaded from the GEO database. 
There are multiple ways to generate transcriptomic data: from raw count files (by performing Deseq2 normalization), directly from the GEO database, or manually from the GEO web page. 
Additionally, the code performs DEG analysis using limma, visualizes results with a volcano plot, and conducts enrichment analysis with g:Profiler.
