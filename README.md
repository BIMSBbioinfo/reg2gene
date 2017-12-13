<a name="logo"/>
<div align="center">
<img src="https://github.com/IngaPa/reg2gene/blob/master/pkg/inst/hex-reg2gene.png" alt="hex Logo"  ></img>
</a>
</div>


---
title: "reg2gene"
author: Inga Patarcic, Altuna Akalin
output: github_document

---

# reg2gene project
https://github.com/BIMSBbioinfo/reg2gene
[![Build Status](https://travis-ci.org/BIMSBbioinfo/reg2gene.svg?branch=master)](https://travis-ci.org/BIMSBbioinfo/reg2gene)
![codecov.io](https://codecov.io/github/BIMSBbioinfo/reg2gene/coverage.svg?branch=master)






This repo contains the code for predicting the target genes for regulatory elements

It contains the following:


-  regulatory activity quantification per pre-defined region: 
    - the code used to quantify enhancer activties across different markers and studies: *regActivity()*
    - the code used to quantify gene expression across different studies: *bwToGeneExp()*
    - the code used to match enhancers to potential target genes within pre-defined window around gene TSSes: *regActivityAroundTSS()*
- enhancer target prediction (targetPrediction)
    - the code for target prediction methods: *associateReg2Gene()*
- benchmark predicted interactions:
   - the code for benchmarking predicted interactions given some benchmark dataset: *benchmarkGI()*
- performance metrics (performance):
    -  the code for performance metrics via calculating and reporting confusion matrix for tested benchmarked dataset: *confusionMatrix()* + *filterPreBenchGI()*
    
    
    
    # package installation
    
        library(devtools)
        install_github("BIMSBbioinfo/reg2gene")
    

