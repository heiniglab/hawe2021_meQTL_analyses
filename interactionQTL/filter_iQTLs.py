"""
Created on Tue Aug  4 10:37:58 2020

@author: katharina.schmid
"""

import pandas as pd
#Note: for read_parquet also pyarrow or fastparquet required
import os.path

def filter_iqtl(prefix):
    
    #Get the number of tested iQTLs
    numTests = 0
    for i in range(1, 23):
        filename = "%s_chr%i.cis_qtl_pairs.%i.parquet" % (prefix, i, i)
        if os.path.isfile(filename):
            iqtl = pd.read_parquet(filename)
            numTests += iqtl.shape[0]
    
    pval=0.05/numTests
    print("%i total tests leading to p-value of %s" % (numTests,pval))
    
    all_iqtl = None
    summary = pd.DataFrame()
    for i in range(1, 23):
        filename = "%s_chr%i.cis_qtl_pairs.%i.parquet" % (prefix, i, i)
        if os.path.isfile(filename):
            iqtl = pd.read_parquet(filename)
            filtered = iqtl[iqtl.pval_gi < pval]
            print("chromosome %i: %i tested %i iQTL" %(i, iqtl.shape[0], filtered.shape[0]))
            summary = summary.append(pd.DataFrame(data={"chromosome":[i], "ntested":[iqtl.shape[0]], "iQTL":[filtered.shape[0]]}))
            if all_iqtl is None:
                all_iqtl = filtered
            else:
                all_iqtl = all_iqtl.append(filtered)
        else:
            print("chromosome %i: file missing" %(i))
    summary.to_csv("%s_summary.csv" % prefix)
    all_iqtl.to_csv("%s_filtered.csv" % prefix)
    return None

print("KORA - BMI")
filter_iqtl(prefix="results/current/tensorqtl/interactions_bmi")
print("KORA - CD4T")
filter_iqtl(prefix="results/current/tensorqtl/interactions_CD4T")
print("KORA - CD8T")
filter_iqtl(prefix="results/current/tensorqtl/interactions_CD8T")
print("KORA - Mono")
filter_iqtl(prefix="results/current/tensorqtl/interactions_Mono")
print("KORA - sex")
filter_iqtl(prefix="results/current/tensorqtl/interactions_sex")
print("KORA - smoking")
filter_iqtl(prefix="results/current/tensorqtl/interactions_smoking")

print("LOLIPOP - BMI")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_bmi")
print("LOLIPOP - CD4T")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_CD4T")
print("LOLIPOP - CD8T")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_CD8T")
print("LOLIPOP - Mono")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_Mono")
print("LOLIPOP - sex")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_sex")
print("LOLIPOP - smoking")
filter_iqtl(prefix="results/current/tensorqtl/OmniEE/interactions_smoking")
