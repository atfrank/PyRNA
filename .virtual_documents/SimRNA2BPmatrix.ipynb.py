from PyRNA import SimRNAfile2bp, dot2bp, visualize_structure, basepair_matrix2CT
import pandas as pd
from sklearn.metrics import recall_score, precision_score, matthews_corrcoef


# Sensitivity, hit rate, recall, or true positive rate: TPR = TP/(TP+FN)
# Specificity or true negative rate: TNR = TN/(TN+FP) 
# Precision or positive predictive value: PPV = TP/(TP+FP)
# Matthews correlation coefficient (MCC): https://scikit-learn.org/stable/modules/generated/sklearn.metrics.matthews_corrcoef.html#sklearn.metrics.matthews_corrcoef


bp_matrix_ref = SimRNAfile2bp("data/SimRNA/reference.ss")
bp_matrix_pre = SimRNAfile2bp("data/SimRNA/predicted.ss")


ref = pd.DataFrame.from_records(bp_matrix_ref).values.flatten()
pre = pd.DataFrame.from_records(bp_matrix_pre).values.flatten()
tpr, ppv, mcc = recall_score(y_true = ref, y_pred=pre), precision_score(y_true = ref, y_pred=pre), matthews_corrcoef(y_true = ref, y_pred=pre)
print(tpr, ppv, mcc)


ct = basepair_matrix2CT(bp_matrix_ref)
ct.columns = ["get_ipython().run_line_magic("s"%ct.shape[0],", " \"\", \"\", \"\", \"\", \"\"]")
ct.to_csv("data/SimRNA/reference.ct", sep = " ", index=False)


visualize_structure(bp_matrix_ref, label="FSW-reference")
visualize_structure(bp_matrix_pre, label="FSW-predicted")
