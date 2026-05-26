import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
import numpy as np
import matplotlib.pyplot as plt

# 读取CSV文件
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("train.csv",encoding="GBK")

# 将数据集划分为训练集(70%)和测试集(30%)
X = data.iloc[:, :8]  # 前8列作为特征
y = data.iloc[:, 8]   # 第9列作为目标变量

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 构建随机森林模型
rf_model = RandomForestRegressor(n_estimators=500, random_state=42)
rf_model.fit(X_train, y_train)

# 对测试集进行预测
predictions = rf_model.predict(X_test)

# 计算MSE, RMSE和MAE
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_test, predictions)

# 打印MSE, RMSE和MAE
print(f"MSE: {mse}")
print(f"RMSE: {rmse}")
print(f"MAE: {mae}")

# 可视化真实值和预测值之间的比较
plt.figure(figsize=(10, 6))
plt.scatter(y_test, predictions, edgecolors=(0, 0, 0))
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=4)
plt.xlabel('实际值')
plt.ylabel('预测值')
plt.title('实际值和预测值对比图')
plt.show()

import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import SGDRegressor
from sklearn.metrics import mean_squared_error, mean_absolute_error
from sklearn.preprocessing import StandardScaler
import numpy as np
import matplotlib.pyplot as plt

# 读取CSV文件
plt.rcParams['font.sans-serif'] = ['SimHei']

# 读取CSV文件
data = pd.read_csv("train.csv",encoding="GBK")

# 将数据集划分为训练集(70%)和测试集(30%)
X = data.iloc[:, :8]  # 前8列作为特征
y = data.iloc[:, 8]   # 第9列作为目标变量

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# 对数据进行标准化
scaler = StandardScaler()
X_train = scaler.fit_transform(X_train)
X_test = scaler.transform(X_test)

# 构建基于梯度下降的线性回归模型
sgd_model = SGDRegressor(max_iter=1000, tol=1e-3, random_state=42)
sgd_model.fit(X_train, y_train)

# 对测试集进行预测
predictions = sgd_model.predict(X_test)

# 计算MSE, RMSE和MAE
mse = mean_squared_error(y_test, predictions)
rmse = np.sqrt(mse)
mae = mean_absolute_error(y_test, predictions)

# 打印MSE, RMSE和MAE
print(f"MSE: {mse}")
print(f"RMSE: {rmse}")
print(f"MAE: {mae}")

# 可视化真实值和预测值之间的比较
plt.figure(figsize=(10, 6))
plt.scatter(y_test, predictions, color='red', edgecolors='w')  # 改为红色
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], color='blue', linestyle='--', linewidth=2)  # 参考线
plt.xlabel('Actual')
plt.ylabel('Predicted')
plt.title('Actual vs Predicted')
plt.show()

