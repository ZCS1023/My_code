import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import SGDRegressor
from sklearn.preprocessing import StandardScaler
import joblib

# 读取训练数据
train_data = pd.read_csv("train.csv", encoding="GBK")

# 分离特征和标签
X_train = train_data.drop('洪水概率', axis=1)
y_train = train_data['洪水概率']

# 标准化特征
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)

# 训练基于梯度下降的线性回归模型
model = SGDRegressor()
model.fit(X_train_scaled, y_train)

# 保存模型和标准化器
joblib.dump(model, 'linear_regression_model.pkl')
joblib.dump(scaler, 'scaler.pkl')

import pandas as pd
from sklearn.linear_model import SGDRegressor
from sklearn.preprocessing import StandardScaler
import joblib

# 读取新数据文件
new_data = pd.read_csv("test.csv", encoding="GBK")

# 分离特征
X_new = new_data.drop('洪水概率', axis=1)

# 加载已训练的模型和标准化器
model = joblib.load('linear_regression_model.pkl')
scaler = joblib.load('scaler.pkl')

# 标准化新数据
X_new_scaled = scaler.transform(X_new)

# 使用模型进行预测
predicted_probabilities = model.predict(X_new_scaled)

# 将预测的概率填充到原始数据的最后一列
new_data['洪水概率'] = predicted_probabilities

# 保存填充了预测概率的新数据到新的 CSV 文件
new_data.to_csv("new_data_with_predictions.csv", index=False, encoding="GBK")

print("预测完成，结果已保存到 new_data_with_predictions.csv 文件中")
