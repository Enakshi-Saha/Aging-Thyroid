# Compute liger (bonobo-panda) networks for GTEx thyroid samples

import os
#import s3fs
import pandas as pd
import numpy as np
#from psutil import *
import netZooPy
from netZooPy.panda import Panda
from netZooPy.ligress import Ligress

# PANDA on CPU (GSE165724)
expression_file = "/home/esaha/Aging_thyroid/data/validation_data/GSE165724/GSE165724_expression_normalThyroid.txt"
priors_table_file = "/home/esaha/Aging_thyroid/data/validation_data/GSE165724/GSE165724_normalThyroid_motif.csv"
e = pd.read_csv(expression_file, sep = "\t", index_col=0)
ppi_file = "/home/ubuntu/prior_large/ppi_997.txt"
output_liger = "./BonoboPanda_GSE165724/"
ligress_obj = Ligress(e, priors_table_file, ppi_file = ppi_file, output_folder = output_liger, mode_process = 'intersection')
ligress_obj.run_ligress(keep_coexpression=False, delta = 0.3, tune_delta = True, precision = 'single')

# PANDA on CPU (GSE29315)
expression_file = "/home/esaha/Aging_thyroid/data/validation_data/GSE29315/GSE29315_expression.txt"
priors_table_file = "/home/esaha/Aging_thyroid/data/validation_data/GSE29315/GSE29315_motif.csv"
e = pd.read_csv(expression_file, sep = "\t", index_col=0)
ppi_file = "/home/ubuntu/prior_large/ppi_997.txt"
output_liger = "./BonoboPanda_GSE29315/"
ligress_obj = Ligress(e, priors_table_file, ppi_file = ppi_file, output_folder = output_liger, mode_process = 'intersection')
ligress_obj.run_ligress(keep_coexpression=False, delta = 0.3, tune_delta = True, precision = 'single')