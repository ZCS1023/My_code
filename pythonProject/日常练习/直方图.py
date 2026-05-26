import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# 加载数据
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("submit.csv",encoding="GBK")

# 获取概率列
probabilities = data.iloc[:, 1]

# 绘制频数直方图
plt.figure(figsize=(10, 6))
count, bins, ignored = plt.hist(probabilities, bins=30, density=False, alpha=0.6, color='skyblue', edgecolor='black')

# 计算概率列的均值和标准差
mean = np.mean(probabilities)
std = np.std(probabilities)

# 绘制正态分布曲线
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = stats.norm.pdf(x, mean, std) * len(probabilities) * (bins[1] - bins[0])  # 调整正态分布曲线以匹配频数
plt.plot(x, p, 'k', linewidth=2)

# 在每个区间添加竖线
for bin_edge in bins:
    plt.axvline(bin_edge, color='black', linewidth=0.5, linestyle='--')

# 添加标题和标签
title = "拟合结果: 均值 = %.2f,  标准差 = %.2f" % (mean, std)
plt.title(title)
plt.xlabel('概率')
plt.ylabel('频数')

# 显示图形
plt.show()