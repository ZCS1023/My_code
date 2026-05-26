import pandas as pd
from sklearn.ensemble import RandomForestRegressor
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("train.csv",encoding="GBK")

# 提取指标得分和概率列
X = data.iloc[:, 1:21]  # 第2列到第21列是20个指标
y = data.iloc[:, 21]    # 第22列是发生的概率

# 创建随机森林回归模型
model = RandomForestRegressor(n_estimators=100, random_state=42)

# 拟合模型
model.fit(X, y)

# 获取特征重要性
importances = model.feature_importances_

# 创建一个数据框，包含指标编号和特征重要性
features = data.columns[1:21]  # 获取指标的列名
importance_df = pd.DataFrame({
    'Feature': features,
    'Importance': importances
})

# 按重要性从大到小排序
importance_df = importance_df.sort_values(by='Importance', ascending=False)

# 输出结果
print(importance_df)

# 可视化特征重要性
plt.figure(figsize=(10, 6))
plt.barh(importance_df['Feature'], importance_df['Importance'], color='skyblue')
plt.xlabel('Feature Importance')
plt.ylabel('Feature')
plt.title('Feature Importance Analysis using Random Forest')
plt.gca().invert_yaxis()
plt.show()