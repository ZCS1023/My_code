import pandas as pd
from sklearn.model_selection import train_test_split
from sksurv.ensemble import RandomSurvivalForest
from sksurv.metrics import concordance_index_censored
from sksurv.util import Surv
import numpy as np

# 加载数据
data = pd.read_excel('python用.xlsx')

# 假设第一列至第十三列是因素，第十四列是生存时间，第十五列是最终结局
factors = data.iloc[:, :13]  # 取前13列作为因素
time = data.iloc[:, 13]  # 第14列是生存时间
event = data.iloc[:, 14]  # 第15列是生存结局，1为死亡，0为生存

# 将数据分成训练集和测试集
X_train, X_test, y_train, y_test = train_test_split(factors, Surv.from_arrays(event == 1, time), test_size=0.2)

# 构建随机生存森林模型
rsf = RandomSurvivalForest(n_estimators=100, min_samples_split=10, random_state=42)
rsf.fit(X_train, y_train)

# 获取特征重要性
importances = rsf.feature_importances_
feature_names = X_train.columns

# 输出特征重要性
importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': importances})
print("Feature Importances:")
print(importance_df)

# 设置阈值筛选重要变量（例如，重要性 > 0.01）
threshold = 0.01
important_features = importance_df[importance_df['Importance'] > threshold]['Feature'].tolist()

# 根据筛选的特征重新构建训练集和测试集
X_train_filtered = X_train[important_features]
X_test_filtered = X_test[important_features]

# 使用筛选后的特征重新训练随机生存森林模型
rsf_filtered = RandomSurvivalForest(n_estimators=100, min_samples_split=10, random_state=42)
rsf_filtered.fit(X_train_filtered, y_train)

# 对测试集进行预测
pred_surv = rsf_filtered.predict_survival_function(X_test_filtered)

# 计算C指数（衡量模型的预测能力）
c_index = concordance_index_censored(y_test['event'], y_test['time'], rsf_filtered.predict(X_test_filtered))[0]
print("C-index:", c_index)

# 预测1年、3年、5年的生存概率
times = [365, 1095, 1825]  # 1年、3年、5年对应的天数
survival_probs = [rsf_filtered.predict_survival_function(X_test_filtered, return_array=True)[:, time_idx] for time_idx in times]
