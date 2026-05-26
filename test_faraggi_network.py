# ==============================================
# 1. 导入必要库
# ==============================================
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

# 设置随机种子，保证结果可复现
np.random.seed(42)
torch.manual_seed(42)

# ==============================================
# 2. 模拟生存分析数据（适配Faraggi-Simon网络）
# ==============================================
def simulate_survival_data(n_samples=1000, n_covariates=5):
    """
    模拟带有删失的生存数据：
    - 协变量x：正态分布
    - 真实风险指数：包含非线性项（验证Faraggi-Simon的非线性建模能力）
    - 生存时间t：Weibull分布
    - 删失标识delta：0=删失，1=未删失
    """
    # 1. 生成协变量
    x = np.random.normal(0, 1, size=(n_samples, n_covariates))
    
    # 2. 真实风险指数（非线性组合，模拟真实场景）
    # 对比经典Cox的线性βᵀx，这里加入平方项和交互项
    true_risk = 0.5 * x[:, 0] + 0.3 * x[:, 1]**2 + 0.4 * x[:, 2] * x[:, 3] - 0.2 * x[:, 4]
    
    # 3. 生成生存时间（Weibull分布，生存分析常用分布）
    shape = 1.5  # Weibull形状参数
    scale = np.exp(-true_risk / shape)  # 尺度参数与风险指数关联
    survival_time = scale * (-np.log(np.random.uniform(0, 1, n_samples))) ** (1 / shape)
    
    # 4. 生成删失时间（均匀分布），构造删失标识
    censoring_time = np.random.uniform(0, np.percentile(survival_time, 80), n_samples)
    # 观测时间=生存时间和删失时间的最小值
    observed_time = np.minimum(survival_time, censoring_time)
    # 删失标识：1=未删失（生存时间<=删失时间），0=删失
    delta = (survival_time <= censoring_time).astype(int)
    
    return x, observed_time, delta

# 执行数据模拟
n_samples = 1000
n_covariates = 5
x, time, delta = simulate_survival_data(n_samples, n_covariates)

# 数据标准化（神经网络训练必备，提升收敛速度）
scaler = StandardScaler()
x_scaled = scaler.fit_transform(x)

# 划分训练集和测试集（修正：赋值顺序与train_test_split返回规则一致）
x_train, x_test, time_train, time_test, delta_train, delta_test = train_test_split(
    x_scaled, time, delta, test_size=0.2, random_state=42
)

# 转换为PyTorch张量（GPU/CPU自适应）
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
x_train_tensor = torch.tensor(x_train, dtype=torch.float32).to(device)
x_test_tensor = torch.tensor(x_test, dtype=torch.float32).to(device)
time_train_tensor = torch.tensor(time_train, dtype=torch.float32).to(device)
delta_train_tensor = torch.tensor(delta_train, dtype=torch.float32).to(device)

# ==============================================
# 3. 定义Faraggi-Simon网络模型
# ==============================================
class FaraggiSimonNet(nn.Module):
    """
    Faraggi-Simon网络核心结构：
    - 输入层：接收协变量（维度=n_covariates）
    - 隐藏层：1层（经典设置2-3个节点，这里用3个），sigmoid激活（非线性变换）
    - 输出层：单节点，无激活函数（输出风险指数θ(x;w,b)）
    完全遵循Faraggi和Simon 1995年提出的结构
    """
    def __init__(self, n_covariates, hidden_nodes=3):
        super(FaraggiSimonNet, self).__init__()
        # 输入层 -> 隐藏层
        self.fc1 = nn.Linear(n_covariates, hidden_nodes)
        # 隐藏层激活函数：sigmoid（经典选择）
        self.sigmoid = nn.Sigmoid()
        # 隐藏层 -> 输出层（单节点，输出风险指数）
        self.fc2 = nn.Linear(hidden_nodes, 1)
    
    def forward(self, x):
        """前向传播：计算风险指数θ(x)"""
        out = self.fc1(x)
        out = self.sigmoid(out)
        risk_score = self.fc2(out)  # 输出θ(x)，对应替换Cox模型的βᵀx
        return risk_score.squeeze()  # 压缩维度，便于后续计算

# 初始化模型
model = FaraggiSimonNet(n_covariates=n_covariates, hidden_nodes=3).to(device)

# ==============================================
# 4. 定义修正的Cox偏似然损失函数（Faraggi-Simon核心训练目标）
# ==============================================
def cox_partial_likelihood_loss(risk_scores, time, delta):
    """
    计算修正的Cox偏似然损失（最大化偏似然 = 最小化负的对数偏似然）
    适配PyTorch训练，处理删失数据
    参数：
    - risk_scores: 模型输出的风险指数θ(x)
    - time: 观测生存时间
    - delta: 删失标识（1=未删失，0=删失）
    """
    # 步骤1：按观测时间降序排序（Cox偏似然计算要求）
    _, idx = torch.sort(time, descending=True)
    risk_scores_sorted = risk_scores[idx]
    delta_sorted = delta[idx]
    
    # 步骤2：计算对数风险和的累积和（用torch.cumsum实现）
    exp_risk = torch.exp(risk_scores_sorted)
    cum_exp_risk = torch.cumsum(exp_risk, dim=0)
    log_cum_exp_risk = torch.log(cum_exp_risk)
    
    # 步骤3：计算负的对数偏似然（作为损失函数，越小越好）
    # 仅对未删失样本（delta=1）计算损失
    partial_likelihood = torch.sum(delta_sorted * (risk_scores_sorted - log_cum_exp_risk))
    # 取负数，转为最小化问题（PyTorch优化器默认最小化损失）
    neg_log_partial_likelihood = -partial_likelihood / torch.sum(delta_sorted)  # 归一化
    
    return neg_log_partial_likelihood

# ==============================================
# 5. 模型训练
# ==============================================
def train_model(model, x_train, time_train, delta_train, epochs=500, lr=0.01):
    """
    训练Faraggi-Simon网络
    """
    # 优化器：Adam（收敛稳定，优于SGD）
    optimizer = optim.Adam(model.parameters(), lr=lr)
    # 记录训练损失
    train_loss_history = []
    
    model.train()  # 切换训练模式
    for epoch in range(epochs):
        optimizer.zero_grad()  # 清空梯度
        
        # 前向传播：计算风险分数
        risk_scores_train = model(x_train)
        
        # 计算损失
        loss = cox_partial_likelihood_loss(risk_scores_train, time_train, delta_train)
        
        # 反向传播+参数更新
        loss.backward()
        optimizer.step()
        
        # 记录损失
        train_loss_history.append(loss.item())
        
        # 每50轮打印一次训练状态
        if (epoch + 1) % 50 == 0:
            print(f"Epoch [{epoch+1}/{epochs}], Train Loss: {loss.item():.4f}")
    
    return model, train_loss_history

# 执行训练
epochs = 500
lr = 0.01
trained_model, loss_history = train_model(
    model, x_train_tensor, time_train_tensor, delta_train_tensor, epochs, lr
)

# ==============================================
# 6. 模型预测与结果可视化
# ==============================================
# 6.1 测试集预测
trained_model.eval()  # 切换评估模式
with torch.no_grad():  # 关闭梯度计算，提升速度
    risk_scores_test = trained_model(x_test_tensor).cpu().numpy()
    risk_scores_train = trained_model(x_train_tensor).cpu().numpy()

# 6.2 可视化训练损失变化
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.plot(range(1, epochs+1), loss_history, color='blue', linewidth=2)
plt.xlabel("Epoch", fontsize=10)
plt.ylabel("Training Loss (Neg Log Partial Likelihood)", fontsize=10)
plt.title("Faraggi-Simon Network Training Loss", fontsize=12)
plt.grid(alpha=0.3)

# 6.3 可视化测试集风险分数分布
plt.subplot(1, 2, 2)
plt.hist(risk_scores_test, bins=20, color='orange', alpha=0.7, edgecolor='black')
plt.xlabel("Risk Score (θ(x))", fontsize=10)
plt.ylabel("Number of Samples", fontsize=10)
plt.title("Test Set Risk Score Distribution", fontsize=12)
plt.grid(alpha=0.3)

plt.tight_layout()
plt.show()

# 6.4 输出关键结果
print("\n==================== 训练结果汇总 ====================")
print(f"训练集平均风险分数：{np.mean(risk_scores_train):.4f}")
print(f"测试集平均风险分数：{np.mean(risk_scores_test):.4f}")
print(f"最终训练损失：{loss_history[-1]:.4f}")
print("======================================================")