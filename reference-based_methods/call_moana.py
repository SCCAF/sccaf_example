import warnings
warnings.filterwarnings("ignore")
import sys

from SCCAF import *

ad_ref = sc.read("../write/%s.h5"%sys.argv[1])
ad = sc.read("../write/%s.h5"%sys.argv[2])

from moana.core import ExpMatrix, CellAnnVector
from moana.classify import CellTypeClassifier
def run_moana(ad, ad_ref):
    # ad_ref is the reference
    ad1 = ad_ref[:,ad_ref.var_names.isin(ad.var_names)]
    matrix1 = ExpMatrix(X=ad1.X.T.todense(), genes=ad1.var_names, cells=ad1.obs_names)
    ds = CellAnnVector(ad1.obs['cell'],cells=ad1.obs_names)

    clf = CellTypeClassifier()
    clf.fit(matrix=matrix1, cell_labels = ds)

    matrix = ExpMatrix(X=ad[:,ad1.var_names].X.T.todense(), genes=ad[:,ad1.var_names].var_names,\
                        cells=ad[:,ad1.var_names].obs_names)
    labs = clf.predict(matrix)
    return(labs.values)

import time
t0 = time.time()
x = run_moana(ad, ad_ref)
t1 = time.time() - t0
print(t1)
pd.Series(x).to_csv("%s_%s_moana.csv"%(sys.argv[1],sys.argv[2]))