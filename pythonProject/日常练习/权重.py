import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
import numpy as np

# 读取CSV文件
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("train1.csv",encoding="GBK")

# 假设第1到第10列是指标，第11列是事件发生概率
X = data.iloc[:, :10]
y = data.iloc[:, 10]

# 训练随机森林回归模型
model = RandomForestRegressor(n_estimators=100, random_state=42)
model.fit(X, y)

# 提取特征重要性
feature_importances = model.feature_importances_
features = X.columns

# 归一化权重系数
normalized_importances = feature_importances / np.sum(feature_importances)

# 创建一个DataFrame来存储特征及其对应的权重
importance_df = pd.DataFrame({
    'Feature': features,
    'Importance': normalized_importances
})

# 按权重从大到小排序
importance_df = importance_df.sort_values(by='Importance', ascending=False)

# 输出归一化的特征重要性
print("Normalized Feature Importances from Random Forest (sorted):")
print(importance_df)

# 绘制归一化的特征重要性图表
plt.figure(figsize=(10, 6))
plt.bar(importance_df['Feature'], importance_df['Importance'], color='skyblue')
plt.xlabel('特征')
plt.ylabel('重要性')
plt.title('随机森林法特征重要性')
plt.xticks(rotation=45)
plt.show()