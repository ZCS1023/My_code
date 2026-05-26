import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

# 加载数据
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("train.csv",encoding="GBK")

# 假设这些是你的中文列名
chinese_names = ['id','季风强度', '地形排水', '河流管理', '森林砍伐', '城市化',
                 '气候变化', '大坝质量', '淤积', '农业实践', '侵蚀',
                 '无效防灾', '排水系统', '海岸脆弱性', '滑坡', '流域',
                 '基础设施恶化', '人口得分', '湿地损失', '规划不足', '政策因素',
                 '洪水概率']

# 确保列名匹配数据
data.columns = chinese_names[:len(data.columns)]

# 提取洪水发生概率
probability = data['洪水概率'].values.reshape(-1, 1)

# 进行K-means聚类
kmeans = KMeans(n_clusters=3, random_state=42).fit(probability)
data['风险等级'] = kmeans.labels_

# 计算不同风险组的概率均值
mean_probabilities = data.groupby('风险等级')['洪水概率'].mean()

# 根据概率均值排序风险组
sorted_risk_levels = mean_probabilities.sort_values().index
sorted_mean_probabilities = mean_probabilities[sorted_risk_levels]

# 绘制概率均值曲线图
plt.figure(figsize=(10, 6))
plt.plot(range(3), sorted_mean_probabilities.values, marker='o', linestyle='-')
plt.xticks(range(3), ['Low', 'Medium', 'High'])
plt.xlabel('Risk Level')
plt.ylabel('Average Probability')
plt.title('Average Probability for Different Risk Levels')
plt.grid(True)
plt.show()

from sklearn.tree import DecisionTreeClassifier
import numpy as np
import matplotlib.pyplot as plt

# 提取特征和目标变量
X = data.iloc[:, 1:21].values
y = data['洪水概率'].values

# 训练决策树模型
clf = DecisionTreeClassifier(random_state=42)
clf.fit(X, y)

# 计算特征重要度并进行归一化
importances = clf.feature_importances_
importances_normalized = importances / np.sum(importances)

# 绘制特征重要度横向条形图
features = chinese_names[1:21]
indices = np.argsort(importances_normalized)[::-1]

plt.figure(figsize=(12, 8))
plt.barh(range(X.shape[1]), importances_normalized[indices], align='center')
plt.yticks(range(X.shape[1]), [features[i] for i in indices])
plt.xlabel('Feature Importance')
plt.title('Normalized Feature Importances')
for i, v in enumerate(importances_normalized[indices]):
    plt.text(v, i, f"{v:.2f}", va='center', ha='left')
plt.gca().invert_yaxis()  # 倒序排列
plt.show()

import numpy as np
import matplotlib.pyplot as plt


# 定义一个函数来计算灵敏度分析
def sensitivity_analysis(model, feature_index, X_sample, step=0.1, num_steps=20):
    original_value = X_sample[feature_index]
    probabilities = []

    for i in range(-num_steps, num_steps + 1):
        X_sample[feature_index] = original_value + i * step
        probabilities.append(model.predict([X_sample])[0])

    X_sample[feature_index] = original_value
    return np.array(probabilities)


# 选择一个样本进行灵敏度分析
X_sample = X[0].copy()

# 绘制灵敏度分析图
plt.figure(figsize=(20, 16))
for i, feature_index in enumerate(range(X.shape[1])):
    sensitivities = sensitivity_analysis(clf, feature_index, X_sample)
    plt.subplot(4, 5, i + 1)
    plt.plot(range(-20, 21), sensitivities, label=features[feature_index])
    plt.title(features[feature_index])
    plt.xlabel('Change in Feature Value')
    plt.ylabel('Probability of Flood')
    plt.grid(True)

plt.tight_layout()
plt.show()