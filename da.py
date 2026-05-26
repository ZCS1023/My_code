import numpy as np
from statsmodels.stats.power import NormalIndPower
from statsmodels.stats.proportion import proportion_effectsize, proportions_ztest, proportion_confint

# 计算所需样本数量
effect_size = proportion_effectsize(0.08, 0.09)
analysis = NormalIndPower()
n_per_group = analysis.solve_power(effect_size, power=0.8, alpha=0.05, ratio=1)

# 实验结果分析
count = np.array([1201, 1424])
nobs = np.array([15010, 14990])

# 进行 Z 检验
stat, pval = proportions_ztest(count, nobs, alternative='two-sided')

# 计算两组比例的置信区间
conf_int_1 = proportion_confint(count[0], nobs[0], alpha=0.05, method='normal')
conf_int_2 = proportion_confint(count[1], nobs[1], alpha=0.05, method='normal')

# 转化率
conversion_rate_1 = count[0] / nobs[0] * 100
conversion_rate_2 = count[1] / nobs[1] * 100

# 是否显著
is_significant = "显著" if pval < 0.05 else "不显著"

# 输出结果
print("| 组      | 用户数    | 注册成功  | 转化率   | 95% 置信区间        | p-value | 是否显著 |")
print("| ------ | ------ | ----- | ----- | ------------- | ------- | ---- |")
print(f"| A方案 | {nobs[0]:,} | {count[0]} | {conversion_rate_1:.2f}% | [{conf_int_1[0]*100:.2f}, {conf_int_1[1]*100:.2f}] | —       | —    |")
print(f"| B方案 | {nobs[1]:,} | {count[1]} | {conversion_rate_2:.2f}% | [{conf_int_2[0]*100:.2f}, {conf_int_2[1]*100:.2f}] | {pval:.3f}   | {is_significant}   |")