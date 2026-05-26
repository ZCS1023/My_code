import pandas as pd
import xgboost as xgb
import shap
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import numpy as np
# 设置字体为黑体，以支持中文
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False  # 解决负号显示问题

# 读取CSV文件
data = pd.read_csv("../Card Fraud.csv", encoding="GBK")

# 2. 提取特征和标签
X = data.drop('Fraud', axis=1)
y = data['Fraud']
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 3. 训练模型
model = xgb.XGBRegressor()
model.fit(X_train, y_train)

# 4. 计算 SHAP 值
explainer = shap.TreeExplainer(model)
shap_values = explainer(X_test)

# 5. 可视化 SHAP 值
# 总体特征重要性
shap.summary_plot(shap_values, X_test)

#特征重要性柱状图
shap.plots.bar(shap_values)

# 输出每个特征的 SHAP 值均值
mean_shap_values = np.mean(np.abs(shap_values.values), axis=0)
feature_importance = pd.DataFrame({
    'Feature': X_test.columns,
    'Mean SHAP Value': mean_shap_values
}).sort_values(by='Mean SHAP Value', ascending=False)

print(feature_importance)