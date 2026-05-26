import pandas as pd
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
df = pd.read_csv("train.csv",encoding="GBK")

# 打印所有列名来检查是否正确
print("数据列名:", df.columns)

# 提取20个指标列和事件发生的概率列
indicator_cols = df.columns[1:21]  # 替换为你的20个指标列的实际列名
probability_col = '洪水概率'  # 替换为实际的列名

# 打印指标列名来检查是否正确
print("指标列名:", indicator_cols)

# 提取发生概率列
probabilities = df[[probability_col]]

# K-means 聚类
kmeans = KMeans(n_clusters=3, random_state=42)  # 3个类别：高风险、中风险、低风险
df['风险等级'] = kmeans.fit_predict(probabilities)

# 获取每个聚类的中心并排序
cluster_centers = kmeans.cluster_centers_.flatten()
sorted_indices = np.argsort(cluster_centers)

# 映射聚类结果为高风险、中风险、低风险
risk_mapping = {sorted_indices[0]: '低风险', sorted_indices[1]: '中风险', sorted_indices[2]: '高风险'}
df['风险等级'] = df['风险等级'].map(risk_mapping)

# 标准化数据
scaler = StandardScaler()
indicator_data = scaler.fit_transform(df[indicator_cols])

# PCA主成分分析
pca = PCA(n_components=1)
pca_results = []

for risk_level in ['低风险', '中风险', '高风险']:
    group_data = indicator_data[df['风险等级'] == risk_level]
    if group_data.shape[0] > 0:  # 检查样本数量是否大于0
        pca.fit(group_data)
        pca_results.append(pca.components_[0])
    else:
        # 如果某个风险等级没有数据，添加全零向量
        pca_results.append(np.zeros(len(indicator_cols)))

# 将PCA结果转换为DataFrame
pca_df = pd.DataFrame(pca_results, columns=indicator_cols, index=['低风险', '中风险', '高风险']).T

# 打印PCA结果来检查是否正确
print("PCA结果:\n", pca_df)

# 雷达图绘制
def plot_radar(data, labels, title, colors):
    num_vars = len(labels)
    angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
    angles += angles[:1]

    fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)

    plt.xticks(angles[:-1], labels, color='black', size=10)

    # 隐藏刻度标签
    ax.set_yticklabels([])

    for idx, (risk_level, color) in enumerate(zip(data.columns, colors)):
        values = data[risk_level].values.flatten().tolist()
        values += values[:1]
        ax.plot(angles, values, color=color, linewidth=2, label=risk_level)
        ax.fill(angles, values, color=color, alpha=0.25)

    plt.title(title, size=15, color='black', y=1.1)
    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.show()

# 定义指标标签和颜色
labels = indicator_cols
colors = ['blue', 'orange', 'green']

# 绘制雷达图
plot_radar(pca_df, labels, '不同风险等级的特征重要度雷达图', colors)