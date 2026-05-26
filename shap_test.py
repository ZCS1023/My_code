import pandas as pd
import numpy as np
import seaborn as sns
import torch
import matplotlib.pyplot as plt
import torchtuples as tt
from sklearn.preprocessing import LabelEncoder, StandardScaler
from pycox.models import CoxPH
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from pycox.evaluation import EvalSurv
from sklearn.calibration import calibration_curve
from sklearn.metrics import roc_auc_score, roc_curve, auc

'''一、数据预处理部分'''

# 读取训练集和测试集
train_df = pd.read_csv(r'C:\Users\zheng\Documents\code\PycharmProjects\train_from_R.csv')
test_df = pd.read_csv(r'C:\Users\zheng\Documents\code\PycharmProjects\test_from_R.csv')

# 定义分类变量映射
category_mappings = {
    'age': {1: "<50", 2: "50-59", 3: "60-69", 4: "70-79", 5: ">80"},
    'sex': {1: "Female", 2: "Male"},
    'race': {1: "White", 2: "Black", 3: "Other"},
    'stage_T': {1: "T0", 2: "T3", 3: "T4", 4: "TX"},
    'stage_N': {1: "N0", 2: "N1", 3: "N2", 4: "N3", 5: "NX"},
    'stage_M': {1: "M0", 2: "M1a", 3: "M1b"},
    'grade': {1: "I", 2: "II", 3: "III", 4: "IV"},
    'laterality': {1: "Left", 2: "Right", 3: "Paire site"},
    'chemotherapy': {1: "Yes", 2: "No/Unknown"},
    'marital': {1: "Married", 2: "Unknown"},
    'income': {
        1: "<40000", 2: "40000-59999", 3: "60000-79999",
        4: "80000-99999", 5: "100000-119999", 6: ">120000"
    },
    'rural.urban': {1: "Metropolitan", 2: "Nonmetropolitan"}
}

# 预处理函数（用于训练集和测试集）
def preprocess_data(df, fitted_encoders=None, fitted_scaler=None):
    df = df.copy()

    # 分类变量映射
    for col, mapping in category_mappings.items():
        df[col] = df[col].map(mapping).astype('category')

    # 独热编码
    categorical_cols = list(category_mappings.keys())
    df = pd.get_dummies(df, columns=categorical_cols, drop_first=True)

    # 标准化连续变量
    if fitted_scaler is None:
        scaler = StandardScaler()
        df['tumor.size'] = scaler.fit_transform(df[['tumor.size']])
    else:
        df['tumor.size'] = fitted_scaler.transform(df[['tumor.size']])
        scaler = fitted_scaler

    return df, scaler

# 分别预处理训练集和测试集
df_train, scaler = preprocess_data(train_df)
df_test, _ = preprocess_data(test_df, fitted_scaler=scaler)

# 准备PyCox格式数据
time_col = 'time'
event_col = 'status'

def get_target(df):
    return (df[time_col].values, df[event_col].values)

y_train = get_target(df_train)
y_test = get_target(df_test)

cols_to_drop = [time_col, event_col]
feature_names = df_train.drop(cols_to_drop, axis=1).columns.tolist()
x_train = df_train.drop(cols_to_drop, axis=1).values.astype('float32')
x_test = df_test.drop(cols_to_drop, axis=1).values.astype('float32')

'''二、模型的定义和训练'''

# 定义和训练DeepSurv模型
in_features = x_train.shape[1]
num_nodes = [32, 32]
batch_norm = True
dropout = 0.1
output_bias = False

net = tt.practical.MLPVanilla(
    in_features, num_nodes, 1, batch_norm, dropout, output_bias=output_bias)

model = CoxPH(net, tt.optim.Adam)
batch_size = 256
epochs = 100
callbacks = [tt.callbacks.EarlyStopping()]

log = model.fit(x_train, y_train, batch_size, epochs, callbacks, val_data=(x_test, y_test))
_ = model.compute_baseline_hazards()

import shap
import warnings
warnings.filterwarnings('ignore')

'''三、变量重要度排序 - SHAP版本'''

def calculate_shap_importance(model, x_train, x_test, feature_names, n_samples=100):
    """
    计算DeepSurv模型的SHAP值
    参数：
        model: 训练好的CoxPH模型
        x_train: 训练数据（用于背景分布）
        x_test: 测试数据（用于解释）
        feature_names: 特征名称列表
        n_samples: 背景样本数（SHAP需要背景分布）
    """
    
    print("正在计算SHAP值...")
    
    # 1. 准备背景数据（随机抽样以减少计算时间）
    if n_samples < len(x_train):
        background_idx = np.random.choice(len(x_train), n_samples, replace=False)
        background = x_train[background_idx]
    else:
        background = x_train
    
    print(f"背景数据形状: {background.shape}")
    print(f"解释数据形状: {x_test.shape}")
    
    # 2. 创建SHAP解释器
    # 注意：DeepSurv模型需要包装预测函数
    def model_predict(x):
        """包装模型预测函数，返回风险得分"""
        # CoxPH模型预测风险得分
        return model.predict(x).flatten()
    
    # 创建KernelExplainer（适用于任何模型）
    explainer = shap.KernelExplainer(
        model_predict,
        background,
        link="identity"
    )
    
    # 3. 计算SHAP值（可以分批计算以避免内存问题）
    print("计算SHAP值（这可能需要一些时间）...")
    
    # 限制样本数以加快计算
    if len(x_test) > 500:
        print("测试集较大，随机抽取500个样本进行SHAP分析...")
        test_idx = np.random.choice(len(x_test), min(500, len(x_test)), replace=False)
        x_test_sampled = x_test[test_idx]
    else:
        x_test_sampled = x_test
    
    # 计算SHAP值
    shap_values = explainer.shap_values(
        x_test_sampled,
        nsamples=100,  # 减少计算量
        silent=True
    )
    
    print(f"SHAP值计算完成，形状: {shap_values.shape}")
    
    return explainer, shap_values, x_test_sampled

def plot_shap_summary(shap_values, features, feature_names):
    """绘制SHAP摘要图"""
    
    plt.figure(figsize=(12, 8))
    
    # 创建SHAP摘要图
    shap.summary_plot(
        shap_values, 
        features, 
        feature_names=feature_names,
        max_display=20,  # 显示最重要的20个特征
        show=False
    )
    
    plt.title("SHAP Feature Importance Summary", fontsize=14, pad=20)
    plt.tight_layout()
    plt.show()

def plot_shap_summary_bar(shap_values, feature_names):
    """绘制SHAP重要性条形图"""
    
    # 计算平均绝对SHAP值作为重要性
    shap_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': np.abs(shap_values).mean(axis=0)
    }).sort_values('importance', ascending=False)
    
    plt.figure(figsize=(10, 12))
    
    # 绘制条形图
    colors = plt.cm.viridis(np.linspace(0.3, 0.9, len(shap_importance)))
    
    bars = plt.barh(
        shap_importance['feature'][:15],  # 显示前15个
        shap_importance['importance'][:15],
        color=colors,
        edgecolor='black',
        alpha=0.8
    )
    
    plt.xlabel('Mean |SHAP value| (average impact on model output)', fontsize=12)
    plt.title('Top 15 Feature Importance based on SHAP Values', fontsize=14, pad=20)
    plt.gca().invert_yaxis()  # 最重要的在顶部
    
    # 添加数值标签
    for i, (bar, imp) in enumerate(zip(bars, shap_importance['importance'][:15])):
        plt.text(imp + 0.001, bar.get_y() + bar.get_height()/2,
                f'{imp:.4f}', va='center', fontsize=10)
    
    plt.grid(axis='x', alpha=0.3, linestyle='--')
    plt.tight_layout()
    plt.show()
    
    return shap_importance

def plot_shap_dependence_plots(shap_values, features, feature_names, top_n=5):
    """绘制Top N特征的依赖图"""
    
    # 获取最重要的特征
    shap_importance = pd.DataFrame({
        'feature': feature_names,
        'importance': np.abs(shap_values).mean(axis=0)
    }).sort_values('importance', ascending=False)
    
    top_features = shap_importance['feature'].head(top_n).tolist()
    
    print(f"\n绘制最重要的 {top_n} 个特征的SHAP依赖图:")
    
    # 创建子图
    fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    axes = axes.flatten()
    
    for i, (idx, row) in enumerate(shap_importance.head(6).iterrows()):
        feat_name = row['feature']
        feat_idx = feature_names.index(feat_name)
        
        # 获取特征值
        feature_values = features[:, feat_idx]
        
        # 绘制依赖图
        shap.dependence_plot(
            feat_idx,
            shap_values,
            features,
            feature_names=feature_names,
            ax=axes[i],
            show=False,
            alpha=0.5,
            dot_size=12
        )
        
        axes[i].set_title(f"SHAP Dependence: {feat_name}", fontsize=11)
        axes[i].set_xlabel(feat_name, fontsize=10)
        axes[i].set_ylabel('SHAP Value (impact on risk score)', fontsize=10)
        
        # 添加回归线
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            feature_values, shap_values[:, feat_idx]
        )
        x_range = np.linspace(feature_values.min(), feature_values.max(), 100)
        axes[i].plot(x_range, slope * x_range + intercept, 
                    'r--', alpha=0.8, linewidth=2,
                    label=f'slope={slope:.3f}, r²={r_value**2:.3f}')
        axes[i].legend(fontsize=9)
        
        # 只在最后一个子图显示图例
        if i < 5:
            axes[i].get_legend().remove()
    
    # 隐藏多余的子图
    for i in range(len(top_features), len(axes)):
        axes[i].set_visible(False)
    
    plt.suptitle('SHAP Dependence Plots for Top Features', fontsize=14, y=1.02)
    plt.tight_layout()
    plt.show()

def aggregate_shap_by_original_feature_simple(shap_values, feature_names, category_mappings):
    """
    简化版的SHAP值聚合函数
    """
    # 获取原始分类变量名列表
    original_categorical_vars = list(category_mappings.keys())
    
    # 创建从编码后特征名到原始特征名的映射
    feature_to_original = {}
    
    for feat in feature_names:
        # 检查是否为分类变量的编码特征
        found = False
        for orig_var in original_categorical_vars:
            if feat.startswith(orig_var + '_'):
                feature_to_original[feat] = orig_var
                found = True
                break
        
        # 如果不是分类变量的编码特征，保持原样
        if not found:
            feature_to_original[feat] = feat
    
    # 根据原始特征名聚合SHAP值
    aggregated = {}
    for feat_idx, feat in enumerate(feature_names):
        orig_feat = feature_to_original[feat]
        shap_vals = shap_values[:, feat_idx]
        
        if orig_feat not in aggregated:
            aggregated[orig_feat] = {
                'abs_shap_values': [],
                'shap_values': []
            }
        
        aggregated[orig_feat]['abs_shap_values'].append(np.abs(shap_vals))
        aggregated[orig_feat]['shap_values'].append(shap_vals)
    
    # 计算统计量
    result = []
    for orig_feat, vals in aggregated.items():
        # 合并所有子特征的SHAP值
        all_abs_shap = np.concatenate(vals['abs_shap_values'])
        all_shap = np.concatenate(vals['shap_values'])
        
        mean_abs_shap = np.mean(all_abs_shap)
        std_abs_shap = np.std(all_abs_shap)
        mean_shap = np.mean(all_shap)
        
        # 确定方向
        if mean_shap > 0.01:  # 加入阈值避免浮点误差
            direction = 'Positive'
        elif mean_shap < -0.01:
            direction = 'Negative'
        else:
            direction = 'Neutral'
        
        result.append({
            'feature': orig_feat,
            'shap_importance_mean': mean_abs_shap,
            'shap_importance_std': std_abs_shap,
            'mean_shap_value': mean_shap,
            'direction': direction,
            'n_subfeatures': len(vals['abs_shap_values'])
        })
    
    # 转为 DataFrame 并排序
    agg_df = pd.DataFrame(result).sort_values(by='shap_importance_mean', ascending=False).reset_index(drop=True)
    return agg_df

# 主调用部分
print("\n========== SHAP特征重要性分析 ==========")

# 计算SHAP值
explainer, shap_values, x_test_sampled = calculate_shap_importance(
    model, x_train, x_test, feature_names, n_samples=200
)

# 1. 绘制SHAP摘要图
print("\n1. 绘制SHAP摘要图...")
plot_shap_summary(shap_values, x_test_sampled, feature_names)

# 2. 绘制SHAP重要性条形图
print("\n2. 绘制SHAP重要性条形图...")
shap_importance_df = plot_shap_summary_bar(shap_values, feature_names)

print("\nTop 10 特征重要性（基于SHAP）:")
print(shap_importance_df.head(10))

# 3. 绘制依赖图
print("\n3. 绘制SHAP依赖图...")
plot_shap_dependence_plots(shap_values, x_test_sampled, feature_names, top_n=6)

# 4. 聚合到原始变量级别
def aggregate_shap_by_original_feature_simple(shap_values, feature_names, category_mappings):
    """
    简化版的SHAP值聚合函数
    """
    # 获取原始分类变量名列表
    original_categorical_vars = list(category_mappings.keys())
    
    # 创建从编码后特征名到原始特征名的映射
    feature_to_original = {}
    
    for feat in feature_names:
        # 检查是否为分类变量的编码特征
        found = False
        for orig_var in original_categorical_vars:
            if feat.startswith(orig_var + '_'):
                feature_to_original[feat] = orig_var
                found = True
                break
        
        # 如果不是分类变量的编码特征，保持原样
        if not found:
            feature_to_original[feat] = feat
    
    # 根据原始特征名聚合SHAP值
    aggregated = {}
    for feat_idx, feat in enumerate(feature_names):
        orig_feat = feature_to_original[feat]
        shap_vals = shap_values[:, feat_idx]
        
        if orig_feat not in aggregated:
            aggregated[orig_feat] = {
                'abs_shap_values': [],
                'shap_values': []
            }
        
        aggregated[orig_feat]['abs_shap_values'].append(np.abs(shap_vals))
        aggregated[orig_feat]['shap_values'].append(shap_vals)
    
    # 计算统计量
    result = []
    for orig_feat, vals in aggregated.items():
        # 合并所有子特征的SHAP值
        all_abs_shap = np.concatenate(vals['abs_shap_values'])
        all_shap = np.concatenate(vals['shap_values'])
        
        mean_abs_shap = np.mean(all_abs_shap)
        std_abs_shap = np.std(all_abs_shap)
        mean_shap = np.mean(all_shap)
        
        # 确定方向
        if mean_shap > 0.01:  # 加入阈值避免浮点误差
            direction = 'Positive'
        elif mean_shap < -0.01:
            direction = 'Negative'
        else:
            direction = 'Neutral'
        
        result.append({
            'feature': orig_feat,
            'shap_importance_mean': mean_abs_shap,
            'shap_importance_std': std_abs_shap,
            'mean_shap_value': mean_shap,
            'direction': direction,
            'n_subfeatures': len(vals['abs_shap_values'])
        })
    
    # 转为 DataFrame 并排序
    agg_df = pd.DataFrame(result).sort_values(by='shap_importance_mean', ascending=False).reset_index(drop=True)
    return agg_df

# 主调用部分 - 使用简化版本
print("\n4. 聚合到原始变量级别...")
shap_agg_df = aggregate_shap_by_original_feature_simple(shap_values, feature_names, category_mappings)

print("\nTop 10 特征重要性（按原始变量聚合）:")
print(shap_agg_df[['feature', 'shap_importance_mean', 'direction', 'mean_shap_value']].head(10))

print("\nTop 10 特征重要性（按原始变量聚合）:")
print(shap_agg_df[['feature', 'shap_importance_mean', 'direction']].head(10))

# 5. 绘制原始变量级别的SHAP重要性
plt.figure(figsize=(12, 10))

top_n_original = min(15, len(shap_agg_df))
plot_data = shap_agg_df.head(top_n_original).copy()
plot_data['feature'] = plot_data['feature'].astype('category')
plot_data['feature'] = plot_data['feature'].cat.reorder_categories(
    plot_data['feature'].tolist()[::-1]  # 让最重要的在顶部
)

# 根据影响方向设置颜色
colors = ['#FF6B6B' if d == 'Positive' else '#4ECDC4' if d == 'Negative' else '#C7C7C7' 
          for d in plot_data['direction']]

bars = plt.barh(
    plot_data['feature'],
    plot_data['shap_importance_mean'],
    xerr=plot_data['shap_importance_std'],
    color=colors,
    capsize=5,
    alpha=0.8,
    edgecolor='black'
)

plt.xlabel('Mean |SHAP value| (Average Impact)', fontsize=12)
plt.ylabel('Covariate', fontsize=12)
plt.title(f'Top {top_n_original} Covariate Importance based on SHAP Values\n(Color indicates direction of impact)', 
          fontsize=14, pad=20)

# 添加图例
from matplotlib.patches import Patch
legend_elements = [
    Patch(facecolor='#FF6B6B', label='Positive Impact (increases risk)'),
    Patch(facecolor='#4ECDC4', label='Negative Impact (decreases risk)'),
    Patch(facecolor='#C7C7C7', label='Mixed Impact')
]
plt.legend(handles=legend_elements, loc='lower right', fontsize=10)

plt.grid(axis='x', alpha=0.3, linestyle='--')
plt.tight_layout()
plt.show()

# 6. 输出详细分析报告
print("\n" + "="*60)
print("SHAP分析报告总结")
print("="*60)

print("\n1. 最重要特征分析:")
for i, (_, row) in enumerate(shap_agg_df.head(5).iterrows()):
    feat = row['feature']
    importance = row['shap_importance_mean']
    direction = row['direction']
    
    if direction == 'Positive':
        effect = "增加患者的死亡风险"
    elif direction == 'Negative':
        effect = "降低患者的死亡风险"
    else:
        effect = "影响方向不确定"
    
    print(f"  {i+1}. {feat}")
    print(f"     重要性得分: {importance:.4f}")
    print(f"     影响方向: {direction} - {effect}")
    print(f"     置信区间: ±{row['shap_importance_std']:.4f}")

print("\n2. 临床意义:")
print("   - 红色特征：值越高，患者风险越高")
print("   - 蓝色特征：值越高，患者风险越低")
print("   - SHAP值的大小表示特征对预测的影响力大小")