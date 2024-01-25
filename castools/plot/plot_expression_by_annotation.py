import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

expression = pd.read_csv("SL4-2_bulk_sc_merged.tsv", delimiter = "\t")
print(expression.head())
mapping = pd.read_csv("mapped_locations_chromHMM_100821.tsv", delimiter = "\t")
print(mapping.head())
expression = expression.merge(mapping, on = ["tBC"])
print(expression.shape)
print(expression.head())

sns.set(rc={"figure.figsize":(6, 14)}) #width=8, height=4
plt.figure(0)
g = sns.boxplot(data = expression, y = "exp_x", x = "annotation")
plt.xticks(rotation=45, size = 8)
plt.savefig("bulkmean_annotationsraw.pdf")
plt.figure(1)
g = sns.boxplot(data = expression, y = "mean", x = "annotation")
plt.xticks(rotation=45)
plt.savefig("mean_annotationsraw.pdf")
plt.figure(2)
g = sns.boxplot(data = expression, y = "var", x = "annotation")
plt.xticks(rotation=45)
plt.savefig("var_annotationsraw.pdf")
plt.figure(3)
g = sns.boxplot(data = expression, y = "mu", x = "annotation")
plt.xticks(rotation=45)
plt.savefig("mu_annotationsraw.pdf")
plt.figure(4)
g = sns.boxplot(data = expression, y = "alpha", x = "annotation")
plt.xticks(rotation=45)
plt.savefig("alpha_annotationsraw.pdf")


annotation_dict = {
     "13_Heterochrom/lo": "repressed",
     "11_Weak_Txn": "rest",
     "12_Repressed": "repressed",
     "10_Txn_Elongation": "rest",
     "7_Weak_Enhancer": "enhancer",
     "4_Strong_Enhancer": "enhancer",
     "9_Txn_Transition": "rest",
     "1_Active_Promoter": "promoter",
      "6_Weak_Enhancer": "enhancer",
      "2_Weak_Promoter": "promoter",
      "5_Strong_Enhancer": "enhancer",
      "8_Insulator": "rest",
      "14_Repetitive/CNV": "rest"
}
expression["annotation_summary"] = expression["annotation"].map(annotation_dict)
print(expression.head(20))

plt.figure(5)
sns.boxplot(data = expression, y = "exp_x", x = "annotation_summary")
plt.savefig("bulkmean_annotationsummary.pdf")
plt.figure(6)
sns.boxplot(data = expression, y = "mean", x = "annotation_summary")
plt.savefig("mean_annotationsummary.pdf")
plt.figure(7)
sns.boxplot(data = expression, y = "var", x = "annotation_summary")
plt.savefig("var_annotationsummary.pdf")
plt.figure(8)
sns.boxplot(data = expression, y = "mu", x = "annotation_summary")
plt.savefig("mu_annotationsummary.pdf")
plt.figure(9)
sns.boxplot(data = expression, y = "alpha", x = "annotation_summary")
plt.savefig("alpha_annotationsummary.pdf")

