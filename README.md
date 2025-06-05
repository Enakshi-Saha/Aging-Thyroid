# Sex differences in thyroid aging and their implications in thyroid disorders: insights from gene regulatory networks

This repository contains all the necessary Python code for generating the individual-specific gene regulatory networks and all R code used in analyzing these networks to derive the results in our paper titled "Sex differences in thyroid aging and their implications in thyroid disorders: insights from gene regulatory networks". doi. TBD

## Discovery and Validation Datasets
We used RNA-sequencing data from thyroid tissue samples from the Genotype Tissue Expression (GTEx) Project to identify aging-related gene regulatory changes in males and females. To validate our results on aging we used an RNa-sequencing dataset from the Gene Expression Omnibus (GEO) with accession number GSE165724. 

We used the GTEx data, again to identify gene regulatory patterns associated with Hashimoto's thyroiditis (HT). A microarray dataset from GEO with accession number GSE29315 was used as the validation dataset to validate results on HT-related gene regulatory changes.

We used two other microarray datasets from GEO with accession numbers GSE33630 and GSE27155 as discovery and validation datasets respectively to identify gene regulatory changes associated with anaplastic thyroid carcinoma. 

## Constructing Individual-specific Gene Regulatory Networks
We used gene expression data from the above mentioned discovery and validation datasets to first compute individual-specific gene co-expression networks using BONOBO (https://genome.cshlp.org/content/34/9/1397.short). These co-expression networks were combined with sex-specific TF-gene motif priors (from CIS-BP) and TF protein-protein interaction (PPI) prior network (from StringDB), using PANDA (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0064832) to infer individual-specific gene regulatory networks. PANDA was implemented using python package netZooPy (https://netzoopy.readthedocs.io/en/stable/).

The networks are stored in an Amazon Web Services s3 bucket and will be made available upon reasonable request.

## Code
R code for replicating the analysis are documented in /Code/README.txt

