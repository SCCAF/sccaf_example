import warnings
warnings.filterwarnings("ignore")
import sys

from SCCAF import *

ad_ref = sc.read("../write/%s.h5"%sys.argv[1])
ad = sc.read("../write/%s.h5"%sys.argv[2])

def run_LR(ad, ad_ref):
    # ad_ref is the reference
    ad1 = ad_ref[:,ad_ref.var_names.isin(ad.var_names)]
    
    y_prob, y_pred, y_test, clf, cvsm, accuracy_test = SCCAF_assessment(ad1.X, \
                    ad1.obs['cell'], n=100)
    return(clf.predict(ad[:,ad1.var_names].X))

import time
t0 = time.time()
x = run_LR(ad, ad_ref)
t1 = time.time() - t0
print(t1)
pd.Series(x).to_csv("%s_%s_LogisticRegression.csv"%(sys.argv[1],sys.argv[2]))
