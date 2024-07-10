##chapter 1
import numpy as np
from sklearn.model_selection import train_test_split, cross_validate
from sklearn.ensemble import RandomForestClassifier, AdaBoostClassifier, ExtraTreesClassifier
from sklearn.svm import SVC
from sklearn.neural_network import MLPClassifier
from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score, make_scorer
import lightgbm as lgb
from sklearn.naive_bayes import GaussianNB
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import SGDClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.datasets import make_classification
from sklearn.model_selection import cross_val_predict, train_test_split
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.metrics import roc_curve, auc, roc_auc_score
from sklearn.model_selection import GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC
from sklearn.ensemble import AdaBoostClassifier
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neural_network import MLPClassifier
import lightgbm as lgb
import matplotlib.pyplot as plt
import os
import pandas as pd
from sklearn.ensemble import StackingClassifier
from sklearn.model_selection import StratifiedKFold
import numpy as np
from numpy import interp
from sklearn.model_selection import cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.ensemble import RandomForestClassifier
from mlxtend.classifier import StackingCVClassifier
from sklearn.metrics import *
from sklearn.tree import DecisionTreeClassifier
from sklearn.metrics import confusion_matrix
import random

# 设置 Pandas 显示选项，以便显示所有的行和列
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)


##chapter 2
def getBPM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_BPM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)


        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)

    for i in HCC:
        j = i.replace("_BPM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)

        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)

    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getsBPM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_sBPM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_sBPM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)
    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getEDM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_EDM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_EDM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)
    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getOLM(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        j = i.replace("_OLM.txt", '')

        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(0)
        elif j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif/np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif))/(np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(0)
    for i in HCC:
        j = i.replace("_OLM.txt", '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            train_X.append(motif)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['Motif', 'frequency'])
            motif = tmp_data['frequency'].to_list()
            # motif = motif / np.sum(motif)
            motif = (motif - np.mean(motif)) / np.std(motif)
            # motif = (motif - np.min(motif)) / (np.max(motif) - np.min(motif))
            test_X.append(motif)
            test_y.append(1)
    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


def getCNV(HCC_path, Control_path, X_train, X_test):
    Control = os.listdir(Control_path)
    HCC = os.listdir(HCC_path)
    train_X = []
    train_y = []
    test_X = []
    test_y = []
    for i in Control:
        # j = i.replace('.logR.txt', '')
        j = i.replace('.logR.txt', '').replace('_logR.txt', '')
        if j in X_train:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['log'])
            CNV = tmp_data['log'].to_list()
            train_X.append(CNV)
            train_y.append(0)
        if j in X_test:
            tmp_data = pd.read_csv(Control_path + i, sep='\t', header=None, names=['log'])
            CNV = tmp_data['log'].to_list()
            test_X.append(CNV)
            test_y.append(0)
    for i in HCC:
        # j = i.replace('.logR.txt', '')
        j = i.replace('.logR.txt', '').replace('_logR.txt', '')
        if j in X_train:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['log'])
            CNV = tmp_data['log'].to_list()
            train_X.append(CNV)
            train_y.append(1)
        if j in X_test:
            tmp_data = pd.read_csv(HCC_path + i, sep='\t', header=None, names=['log'])
            CNV = tmp_data['log'].to_list()
            test_X.append(CNV)
            test_y.append(1)
    return np.array(train_X), np.array(train_y), np.array(test_X), np.array(test_y)


##chapter 3

def tStacking_BPM():
    clf1 = RandomForestClassifier(n_estimators=300, random_state=20)
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(n_estimators=300, random_state=20)
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(max_iter=5000, random_state=20)
    clf6 = SVC(C=1000, random_state=9)
    clf7 = KNeighborsClassifier(n_neighbors=1)
    clf9 = lgb.LGBMClassifier(random_state=0, n_estimators=1000, learning_rate=0.001)
    clf10 = AdaBoostClassifier(random_state=0)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf4, clf1, clf2, clf5],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5,
                                 random_state=0
                                 )
    return model


def Stacking_BPM():
    clf1 = RandomForestClassifier(random_state=0, )
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(random_state=49, )
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(random_state=9, max_iter=500)
    clf6 = SVC(random_state=9)
    clf7 = KNeighborsClassifier()
    clf8 = LogisticRegression(max_iter=500)
    clf9 = lgb.LGBMClassifier(random_state=18, n_estimators=300)
    clf10 = AdaBoostClassifier(random_state=5878)  # AdaBoostClassifier(random_state=0)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf1, clf3, clf5, clf8],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5,
                                 random_state=5
                                 )
    return model


def Stacking_sBPM():
    clf1 = RandomForestClassifier(random_state=0, )
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(random_state=49, )
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(random_state=9, max_iter=500)
    clf6 = SVC(random_state=9)
    clf7 = KNeighborsClassifier()
    clf8 = LogisticRegression(max_iter=500)
    clf9 = lgb.LGBMClassifier(random_state=18, n_estimators=300)
    clf10 = AdaBoostClassifier(random_state=5878)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf1, clf3, clf5, clf8],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5, random_state=0
                                 )
    return model


def Stacking_EDM():
    clf1 = RandomForestClassifier(random_state=0, )
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(random_state=49, )
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(random_state=9, max_iter=500)
    clf6 = SVC(random_state=9)
    clf7 = KNeighborsClassifier()
    clf8 = LogisticRegression(max_iter=500)
    clf9 = lgb.LGBMClassifier(random_state=18, n_estimators=300)
    clf10 = AdaBoostClassifier(random_state=5878)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf1, clf3, clf5, clf8],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5, random_state=0
                                 )
    return model


def Stacking_OLM():
    clf1 = RandomForestClassifier(random_state=0, )
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(random_state=49, )
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(random_state=9, max_iter=500)
    clf6 = SVC(random_state=9)
    clf7 = KNeighborsClassifier()
    clf8 = LogisticRegression(max_iter=500)
    clf9 = lgb.LGBMClassifier(random_state=18, n_estimators=300)
    clf10 = AdaBoostClassifier(random_state=5878)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf1, clf3, clf5, clf8],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5, random_state=0
                                 )
    return model


def Stacking_CNV():
    clf1 = RandomForestClassifier(random_state=0, )
    clf2 = GaussianNB()
    clf3 = ExtraTreesClassifier(random_state=49, )
    clf4 = SGDClassifier(random_state=0)
    clf5 = MLPClassifier(random_state=9, max_iter=500)
    clf6 = SVC(random_state=9)
    clf7 = KNeighborsClassifier()
    clf8 = LogisticRegression(max_iter=500)
    clf9 = lgb.LGBMClassifier(random_state=18, n_estimators=300)
    clf10 = AdaBoostClassifier(random_state=5878)
    lr = LogisticRegression(random_state=20)
    model = StackingCVClassifier(classifiers=[clf1, clf3, clf5, clf8],  # 第一层分类器
                                 meta_classifier=lr,  # 第二层分类器
                                 cv=5, random_state=0
                                 )
    return model


##chapter 4   获得训练集及验证集的名字
#
Control = os.listdir('path/to/your/data/')
HCC = os.listdir('path/to/your/data')
all_name = []
False_name = []
postive_name = []
n_label = []

for i in Control:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_BPM", "")
    all_name.append(final_name)
    False_name.append(final_name)
    n_label.append(0)

for i in HCC:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_BPM", "")
    all_name.append(final_name)
    postive_name.append(final_name)
    n_label.append(1)
print(all_name)


print(n_label)


def random_sample_with_seed(data, num_samples, seed):
    random.seed(seed)
    return random.sample(data, num_samples)


def random_sample_with_seed_second(data, num_samples, seed):
    random.seed(seed)

    first_sample = random.sample(data, num_samples)

    second_sample = random.sample(list(filter(lambda x: x not in first_sample, data)), num_samples)

    return first_sample, second_sample


false_name, HCC_false_name = random_sample_with_seed_second(False_name, 10,2)
print('false_name:', false_name)
combine_name = false_name + postive_name

# X_train, X_test, y_train, y_test = train_test_split(all_name, n_label, stratify=n_label,test_size=0.3, random_state=2)
##chapter 5
Control_disease = os.listdir('path/to/your/data/')
disease = os.listdir('path/to/your/data/')
all_name_disease = []
HCC_positive = []
n_label_disease = []

for i in Control_disease:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_BPM", "")
    all_name_disease.append(final_name)
    # HCC_false.append(final_name)
    n_label_disease.append(0)

for i in disease:
    final_name = os.path.splitext(os.path.basename(i))[0]
    final_name = final_name.replace("_BPM", "")
    all_name_disease.append(final_name)
    HCC_positive.append(final_name)
    n_label_disease.append(1)

# HCC_false_name = random_sample_with_seed(HCC_false,10,0)
HCC_positive_name = random_sample_with_seed(HCC_positive, 10, 0)
HCC_combine_name = HCC_positive_name + HCC_false_name
print('HCC_false_name', HCC_false_name)

##chapter 6
# 训练集、验证集路径
BPM_Control_path = 'path/to/your/data/'
BPM_HCC_path = 'path/to/your/data/'
sBPM_Control_path = 'path/to/your/data/'
sBPM_HCC_path = 'path/to/your/data/'
EDM_Control_path = 'path/to/your/data/'
EDM_HCC_path = 'path/to/your/data/'
OLM_Control_path = 'path/to/your/data/'
OLM_HCC_path = 'path/to/your/data/'
CNV_Control_path = 'path/to/your/data/'
CNV_HCC_path = 'path/to/your/data/'

#
BPM_disease_Control_path = 'path/to/your/data/'
BPM_disease_Cancer_path = 'path/to/your/data/'
EDM_disease_Control_path = 'path/to/your/data/'
EDM_disease_Cancer_path = 'path/to/your/data/'
sBPM_disease_Control_path = 'path/to/your/data/'
sBPM_disease_Cancer_path = 'path/to/your/data/'
OLM_disease_Control_path = 'path/to/your/data/'
OLM_disease_Cancer_path = 'path/to/your/data/'
CNV_disease_Control_path = 'path/to/your/data/'
CNV_disease_Cancer_path = 'path/to/your/data/'

train_test = []

BPM_X_train, BPM_y_train, BPM_X_test, BPM_y_test = getBPM(BPM_HCC_path, BPM_Control_path, combine_name, train_test)
EDM_X_train, EDM_y_train, EDM_X_test, EDM_y_test = getEDM(EDM_HCC_path, EDM_Control_path, combine_name, train_test)
sBPM_X_train, sBPM_y_train, sBPM_X_test, sBPM_y_test = getsBPM(sBPM_HCC_path, sBPM_Control_path, combine_name, train_test)
OLM_X_train, OLM_y_train, OLM_X_test, OLM_y_test = getOLM(OLM_HCC_path, OLM_Control_path, combine_name, train_test)
CNV_X_train, CNV_y_train, CNV_X_test, CNV_y_test = getCNV(CNV_HCC_path, CNV_Control_path, combine_name, train_test)


#
disease_train = []
#
BPM_disease_X_train, BPM_disease_y_train, BPM_disease_X_test, BPM_disease_y_test = getBPM(BPM_disease_Cancer_path,
                                                                                          BPM_disease_Control_path,
                                                                                          disease_train,
                                                                                          HCC_combine_name)
EDM_disease_X_train, EDM_disease_y_train, EDM_disease_X_test, EDM_disease_y_test = getEDM(EDM_disease_Cancer_path,
                                                                                          EDM_disease_Control_path,
                                                                                          disease_train,
                                                                                          HCC_combine_name)
sBPM_disease_X_train, sBPM_disease_y_train, sBPM_disease_X_test, sBPM_disease_y_test = getsBPM(sBPM_disease_Cancer_path,
                                                                                               sBPM_disease_Control_path,
                                                                                               disease_train,
                                                                                               HCC_combine_name)
OLM_disease_X_train, OLM_disease_y_train, OLM_disease_X_test, OLM_disease_y_test = getOLM(OLM_disease_Cancer_path,
                                                                                          OLM_disease_Control_path,
                                                                                          disease_train,
                                                                                          HCC_combine_name)
CNV_disease_X_train, CNV_disease_y_train, CNV_disease_X_test, CNV_disease_y_test = getCNV(CNV_disease_Cancer_path,
                                                                                          CNV_disease_Control_path,
                                                                                          disease_train,
                                                                                          HCC_combine_name)
# print('EDM_disease_y_test', EDM_disease_y_test)
# print('BPM_disease_y_test', BPM_disease_y_test)
# print('CNV_disease_y_test', CNV_disease_y_test)

from sklearn.model_selection import LeaveOneOut

# 初始化loo对象
loo = LeaveOneOut()

###融合数据
# 二次stacked的数据
BPM = Stacking_BPM()
EDM = Stacking_EDM()
sBPM = Stacking_sBPM()
OLM = Stacking_OLM()
CNV = Stacking_CNV()

##chapter 7
n_folds = 5  # 改折数
skf = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=5)
df_eval = pd.DataFrame(columns=['Accuracy', 'Precision', 'Recall', 'F1_score', 'auc'])


def training_BaseModel(X_train, y_train, disease_test, skf, clf, blend_train, j):
    blend_disease_j = np.zeros((disease_test.shape[0], 5))
    tprs = []
    mean_fpr = np.linspace(0, 1, 100)
    for i, (train_index, cv_index) in enumerate(skf.split(X_train, y_train)):
        print('Fold [%s]' % (i))
        tr_X = X_train[train_index]
        tr_y = y_train[train_index]
        cv_X = X_train[cv_index]
        cv_y = y_train[cv_index]
        clf.fit(tr_X, tr_y)
        blend_train[[cv_index], j] = clf.predict(cv_X)

        blend_disease_j[:, i] = clf.predict(disease_test)

    return blend_train, blend_disease_j


tprs = []
r = []
mean_fpr1 = np.linspace(0, 1, 100)
blend_train = np.zeros((BPM_X_train.shape[0], 5))
blend_test = np.zeros((BPM_X_test.shape[0], 5))
blend_disease = np.zeros((BPM_disease_X_test.shape[0], 5))
j = 0
for clf, label in zip([BPM, EDM, sBPM, OLM, CNV], ['BPM', 'EDM', 'sBPM', 'OLM', 'CNV']):
    print("Training BaseModel [%s]" % (label))

    if label == 'BPM':
        blend_train, blend_disease_j = training_BaseModel(BPM_X_train, BPM_y_train, BPM_disease_X_test, skf, clf,
                                                          blend_train, j)

    if label == 'EDM':
        blend_train, blend_disease_j = training_BaseModel(EDM_X_train, EDM_y_train, EDM_disease_X_test, skf, clf,
                                                          blend_train,
                                                          j)
    if label == 'sBPM':
        blend_train, blend_disease_j = training_BaseModel(sBPM_X_train, sBPM_y_train, sBPM_disease_X_test, skf, clf,
                                                          blend_train,
                                                          j)
    if label == 'OLM':
        blend_train, blend_disease_j = training_BaseModel(OLM_X_train, OLM_y_train, OLM_disease_X_test, skf, clf,
                                                          blend_train,
                                                          j)
    if label == 'CNV':
        blend_train, blend_disease_j = training_BaseModel(CNV_X_train, CNV_y_train, CNV_disease_X_test, skf, clf,
                                                          blend_train,
                                                          j)

    blend_disease[:, j] = blend_disease_j.mean(1)

    j = j + 1

TP_total = 0
TN_total = 0
FP_total = 0
FN_total = 0

clf1 = RandomForestClassifier(random_state=32, min_samples_split=5)
clf3 = ExtraTreesClassifier(random_state=49, n_estimators=600, min_samples_leaf=8, criterion='entropy')
clf5 = MLPClassifier(random_state=0, max_iter=300, learning_rate_init=0.05, learning_rate='invscaling', solver='sgd')
clf8 = LogisticRegression(max_iter=500)


clf_list = {
    'RF': clf1,
    'ExtraTree': clf3,
    'MLP': clf5,
    'Logistic': clf8,
}


import pandas as pd

datasets = {
    'BPM': (BPM_X_train, BPM_y_train, BPM_disease_X_test, BPM_disease_y_test),
    'EDM': (EDM_X_train, EDM_y_train, EDM_disease_X_test, EDM_disease_y_test),
    'CNV': (CNV_X_train, CNV_y_train, CNV_disease_X_test, CNV_disease_y_test),
    'sBPM': (sBPM_X_train, sBPM_y_train, sBPM_disease_X_test, sBPM_disease_y_test),
    'OLM': (OLM_X_train, OLM_y_train, OLM_disease_X_test, OLM_disease_y_test),
    'blend': (blend_train, BPM_y_train, blend_disease, BPM_disease_y_test),
}

results = []
results_disease = []
results_disease_matrix = []
roc = []
for clf_name, clf in clf_list.items():
    for dataset_name, (X_train, y_train, disease_X_test, disease_y_test) in datasets.items():
        print(f", Dataset: {dataset_name}")
        TP_total = 0
        TN_total = 0
        FP_total = 0
        FN_total = 0
        # CLF_TN = 0
        # CLF_TP = 0
        # CLF_FP = 0
        # CLF_FN = 0

        disease_Accuracy = 0
        disease_Specificity = 0
        disease_Recall = 0
        disease_F1 = 0

        AUC = 0
        fpr = 0
        tpr = 0

        for train_index, test_index in loo.split(X_train):
            print(f"Classifier: {clf_name}, Dataset: {dataset_name}")
            print(train_index, test_index)
            #
            X_train_split, X_test_split = X_train[train_index], X_train[test_index]
            y_train_split, y_test_split = y_train[train_index], y_train[test_index]

            clf.fit(X_train_split, y_train_split)
            y_pred = clf.predict(X_test_split)

            #
            TP = np.sum(np.logical_and(y_test_split == 1, y_pred == 1))
            TN = np.sum(np.logical_and(y_test_split == 0, y_pred == 0))
            FP = np.sum(np.logical_and(y_test_split == 0, y_pred == 1))
            FN = np.sum(np.logical_and(y_test_split == 1, y_pred == 0))
            print('tptnfpfn', TP, TN, FP, FN)

            TP_total += TP
            TN_total += TN
            FP_total += FP
            FN_total += FN
            print('TP_total,TN_total,FP_total,FN_total', TP_total, TN_total, FP_total, FN_total)

            HCC_pred = clf.predict(disease_X_test)
            print('disease_y_test,', disease_y_test)
            print('HCC_pred,', HCC_pred)
            clf_tn, clf_fp, clf_fn, clf_tp = confusion_matrix(disease_y_test, HCC_pred).ravel()

            # CLF_TN += clf_tn
            # CLF_FP += clf_fp
            # CLF_FN += clf_fn
            # CLF_TP += clf_tp

            disease_accuracy = accuracy_score(disease_y_test, HCC_pred)
            disease_recall = recall_score(disease_y_test, HCC_pred)
            disease_specificity = clf_tn / (clf_tn + clf_fp)
            disease_f1 = f1_score(disease_y_test, HCC_pred)
            disease_Accuracy += disease_accuracy
            disease_Specificity += disease_specificity
            disease_Recall += disease_recall
            disease_F1 += disease_f1
            print(f'{clf_name}的TP:', clf_tp)
            print(f'{clf_name}的HCC准确率:', disease_accuracy)
            print(f'{clf_name}的HCC召回率:', disease_recall)

            ###计算AUC的值
            ### 画另一癌症的roc曲线
            clf_disease_prob = clf.predict_proba(disease_X_test)[:, 1]

            np.random.seed(0)

            clf_disease_prob_with_noise = clf_disease_prob + np.random.uniform(low=-(1e-1), high=1e-1,
                                                                               size=clf_disease_prob.shape)

            clf_disease_prob_with_noise = np.round(clf_disease_prob_with_noise, decimals=10)
            #
            clf_disease_fpr, clf_disease_tpr, clf_disease_thresholds = roc_curve(disease_y_test,
                                                                                 clf_disease_prob_with_noise)
            clf_disease_mean_fpr1 = np.linspace(0, 1, 100)
            clf_disease_tprs = []
            clf_disease_tprs.append(interp(clf_disease_mean_fpr1, clf_disease_fpr, clf_disease_tpr))
            clf_disease_tprs[-1][0] = 0.0
            clf_disease_mean_tpr1 = np.mean(clf_disease_tprs, axis=0)
            clf_disease_mean_tpr1[-1] = 1.0
            clf_disease_mean_auc = auc(clf_disease_mean_fpr1, clf_disease_mean_tpr1)
            clf_disease_tprs_mean = np.mean(clf_disease_tprs, axis=0)
            clf_disease_roc_auc = auc(clf_disease_mean_fpr1, clf_disease_mean_tpr1)
            AUC += clf_disease_roc_auc
            fpr += clf_disease_mean_fpr1
            tpr += clf_disease_mean_tpr1

        #
        accuracy = (TP_total + TN_total) / (TP_total + TN_total + FP_total + FN_total)
        specificity = TN_total / (TN_total + FP_total)
        recall = TP_total / (TP_total + FN_total)
        f1 = 2 * TP_total / (2 * TP_total + FP_total + FN_total)
        disease_Accuracy_mean = disease_Accuracy / 20
        disease_Specificity_mean = disease_Specificity / 20
        disease_Recall_mean = disease_Recall / 20
        disease_F1_mean = disease_F1 / 20
        disease_AUC_mean = AUC / 20
        disease_fpr = fpr / 20
        disease_tpr = tpr / 20
        # CLF_FN_mean = CLF_FN / 20
        # CLF_TP_mean = CLF_TP / 20
        # CLF_FP_mean = CLF_FP / 20
        # CLF_TN_mean = CLF_TN / 20

        #
        roc.append([clf_name, disease_fpr, disease_tpr, disease_AUC_mean])
        results.append([clf_name, dataset_name, accuracy, specificity, recall, f1])
        results_disease.append(
            [clf_name, dataset_name, disease_Accuracy_mean, disease_Specificity_mean, disease_Recall_mean,
             disease_F1_mean, disease_AUC_mean])
        # results_disease_matrix.append([clf_name, dataset_name, CLF_TP_mean, CLF_TN_mean, CLF_FP_mean, CLF_FN_mean])

df = pd.DataFrame(results, columns=['Classifier', 'Dataset', 'Accuracy', 'Specificity', 'Recall', 'F1'])
df_disease = pd.DataFrame(results_disease,
                          columns=['Classifier', 'Dataset', 'Accuracy', 'Specificity', 'Recall', 'F1', 'AUC'])
# df_disease_matrix = pd.DataFrame(results_disease_matrix,columns=['Classifier','Dataset','TP_mean','TN_mean','FP_mean','FN_mean'])

# 打印表格
print(df)
print(df_disease)


# roc是一个大列表，里面包含了24个小列表，每个小列表有4个元素：分类器名称，fpr，tpr，auc，每个小列表正好用来画一条roc曲线

def plot_roc_curves(roc, n, clfname):
    plt.figure(figsize=(8, 6))
    plt.rcParams['font.size'] = 12


    plt.plot(roc[n][1], roc[n][2], label='BPM (AUC = %0.4f)' % roc[n][3], color='red')

    plt.plot(roc[n + 1][1], roc[n + 1][2], label='EDM (AUC = %0.4f)' % roc[n + 1][3], color='blue')

    plt.plot(roc[n + 2][1], roc[n + 2][2], label='CNV (AUC = %0.4f)' % roc[n + 2][3], color='green')

    plt.plot(roc[n + 4][1], roc[n + 4][2], label='nOLM (AUC = %0.4f)' % roc[n + 4][3], color='yellow')

    plt.plot(roc[n + 3][1], roc[n + 3][2], label='sBPM (AUC = %0.4f)' % roc[n + 3][3], color='pink')

    plt.plot(roc[n + 5][1], roc[n + 5][2], label='Stacked Feature (AUC = %0.4f)' % roc[n + 5][3], color='purple')

    plt.plot([0, 1], [0, 1], color='navy', lw=2, linestyle='--')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic of {}'.format(clfname))
    plt.legend(loc='lower right')
    # 保存为 JPEG 和 EPS 格式
    plt.savefig("path/to/your/{}.jpg".format(clfname), format="jpeg")
    plt.savefig("path/to/your/{}.eps".format(clfname), format="eps")
    # 显示图像
    plt.show()


plot_roc_curves(roc, 0, 'RF')
plot_roc_curves(roc, 6, 'ExtraTree')
plot_roc_curves(roc, 12, 'MLP')
plot_roc_curves(roc, 18, 'Logistic')










