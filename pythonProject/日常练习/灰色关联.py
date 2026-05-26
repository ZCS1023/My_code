import pandas as pd
import numpy as np

# 读取CSV文件
data = pd.read_csv("train.csv")

# 提取指标得分和概率列
X = data.iloc[:, 1:21]  # 第2列到第21列是20个指标
y = data.iloc[:, 21]    # 第22列是发生的概率

# 归一化函数
def normalize(x):
    return (x - np.min(x)) / (np.max(x) - np.min(x))

# 归一化数据
X_normalized = X.apply(normalize, axis=0)
y_normalized = normalize(y)

# 计算灰色关联度
def calculate_grey_rel_degree(x, y, rho=0.5):
    diff_seq = np.abs(x - y)
    max_diff = np.max(diff_seq)
    min_diff = np.min(diff_seq)
    grey_rel_degree = (min_diff + rho * max_diff) / (diff_seq + rho * max_diff)
    return np.mean(grey_rel_degree)

# 初始化灰色关联度列表
grey_relation_degrees = []

# 计算每个指标与概率之间的灰色关联度
for column in X_normalized.columns:
    xi = X_normalized[column]
    grey_rel_degree = calculate_grey_rel_degree(xi, y_normalized)
    grey_relation_degrees.append((column, grey_rel_degree))

# 将结果转换为数据框
importance_df = pd.DataFrame(grey_relation_degrees, columns=['Feature', 'Grey Relation Degree'])

# 按关联度从大到小排序
importance_df = importance_df.sort_values(by='Grey Relation Degree', ascending=False)

# 输出结果
print(importance_df)

# 可视化灰色关联度
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.barh(importance_df['Feature'], importance_df['Grey Relation Degree'], color='skyblue')
plt.xlabel('Grey Relation Degree')
plt.ylabel('Feature')
plt.title('Grey Relation Degree Analysis')
plt.gca().invert_yaxis()
plt.show()
