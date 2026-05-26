import pandas as pd
import numpy as np
import seaborn as sns
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

'''三、变量重要度排序'''

# 变量重要度排序
def calculate_permutation_importance(model, x, y, feature_names, n_repeats=5, random_state=42):
    
    # 计算排列重要度：通过打乱单个特征的值，观察模型C-index的下降幅度（幅度越大，特征越重要）

    np.random.seed(random_state)
    # 先计算原始模型的C-index
    surv_original = model.predict_surv_df(x)
    durations, events = y
    ev_original = EvalSurv(surv_original, durations, events, censor_surv='km')
    baseline_cindex = ev_original.concordance_td()

    # 初始化重要度数组
    imp_scores = np.zeros((len(feature_names), n_repeats))

    # 对每个特征进行排列验证
    for i, feat in enumerate(feature_names):
        for j in range(n_repeats):
            # 复制数据并打乱当前特征
            x_perm = x.copy()
            x_perm[:, i] = np.random.permutation(x_perm[:, i])

            # 预测并计算打乱后的C-index
            surv_perm = model.predict_surv_df(x_perm)
            ev_perm = EvalSurv(surv_perm, durations, events, censor_surv='km')
            perm_cindex = ev_perm.concordance_td()

            # 重要度 = 原始C-index - 打乱后的C-index（下降幅度）
            imp_scores[i, j] = baseline_cindex - perm_cindex

    # 计算每个特征的平均重要度和标准差
    imp_mean = np.mean(imp_scores, axis=1)
    imp_std = np.std(imp_scores, axis=1)

    # 封装结果
    imp_df = pd.DataFrame({
        'feature': feature_names,
        'perm_importance_mean': imp_mean,
        'perm_importance_std': imp_std
    }).sort_values(by='perm_importance_mean', ascending=False)

    return imp_df

# 合并变量
def aggregate_importance_by_original_feature(imp_df, category_mappings):
    
    # 获取原始分类变量名列表
    original_categorical_vars = set(category_mappings.keys())
    original_continuous_vars = {'tumor.size'}  # 连续变量
    
    # 初始化聚合字典
    aggregated = {}
    
    for _, row in imp_df.iterrows():
        feat = row['feature']
        imp_mean = row['perm_importance_mean']
        imp_std = row['perm_importance_std']
        
        # 判断是否属于某个分类变量（通过前缀匹配）
        matched = False
        for orig_var in original_categorical_vars:
            if feat.startswith(orig_var + '_'):
                if orig_var not in aggregated:
                    aggregated[orig_var] = {'imp_means': [], 'imp_stds': []}
                aggregated[orig_var]['imp_means'].append(imp_mean)
                aggregated[orig_var]['imp_stds'].append(imp_std)
                matched = True
                break
        
        # 如果没匹配到，说明是连续变量或未知变量
        if not matched:
            # 直接保留原名（如 'tumor.size'）
            aggregated[feat] = {
                'imp_means': [imp_mean],
                'imp_stds': [imp_std]
            }
    
    # 聚合策略：对每个原始变量，取 importance 的最大值
    result = []
    for orig_feat, vals in aggregated.items():
        means = np.array(vals['imp_means'])
        stds = np.array(vals['imp_stds'])
        
        # 这里采用最大重要度作为该变量的代表
        idx_max = np.argmax(means)
        agg_mean = means[idx_max]
        agg_std = stds[idx_max]  # 对应最大值的那个 std
        
        result.append({
            'feature': orig_feat,
            'perm_importance_mean': agg_mean,
            'perm_importance_std': agg_std
        })
    
    # 转为 DataFrame 并排序
    agg_df = pd.DataFrame(result).sort_values(by='perm_importance_mean', ascending=False).reset_index(drop=True)
    return agg_df

# 绘制特征重要度图
def plot_covariate_importance_seaborn(imp_df, top_n=15, title="Top Covariate Permutation Importance (DeepSurv)"):
    
    # 使用 seaborn 绘制美观的横向条形图（带误差棒）
    
    # 取 Top N
    plot_df = imp_df.head(top_n).copy()
    plot_df['feature'] = plot_df['feature'].astype('category')
    plot_df['feature'] = plot_df['feature'].cat.reorder_categories(
        plot_df['feature'].tolist()  # 让最重要的在顶部
    )

    plt.figure(figsize=(10, max(6, top_n * 0.4)))  # 自适应高度
    ax = sns.barplot(
        data=plot_df,
        x='perm_importance_mean',
        y='feature',
        xerr=plot_df['perm_importance_std'],
        color='steelblue',
        capsize=0.15,
        errwidth=1.2,
        alpha=0.85
    )
    
    # 美化
    ax.set_xlabel('Mean C-index Drop (Permutation Importance)', fontsize=12)
    ax.set_ylabel('Covariate', fontsize=12)
    ax.set_title(title, fontsize=14, pad=20)
    ax.grid(axis='x', linestyle='--', alpha=0.6)
    sns.despine(left=True, bottom=False)
    
    # 调整布局
    plt.tight_layout()
    plt.show()

print("\n========== 计算协变量排列重要度 ==========")
# 先计算原始（dummy-level）重要度
perm_imp_df_raw = calculate_permutation_importance(
    model, x_train, y_train, feature_names, n_repeats=3, random_state=42
)

# 聚合到原始变量级别
perm_imp_df_agg = aggregate_importance_by_original_feature(perm_imp_df_raw, category_mappings)

print("\nTop 10 排列重要度（按原始变量聚合）：")
print(perm_imp_df_agg.head(10))

# 绘图
plot_covariate_importance_seaborn(perm_imp_df_agg, top_n=15)

'''四、模型的应用'''

# 在训练集和测试集上预测生存函数
surv_train = model.predict_surv_df(x_train)
surv_test = model.predict_surv_df(x_test)

'''五、模型的性能评估'''

# 获取真实事件和生存时间
durations_train, events_train = y_train
durations_test, events_test = y_test

# 定义评估函数
def evaluate_model(surv, durations, events, label=""):
    ev = EvalSurv(surv, durations, events, censor_surv='km')

    # C-index
    cindex = ev.concordance_td()
    print(f"\n[{label}] C-index: {cindex:.3f}")

    # 显示1、3、5年的AUC值
    time_points = [12, 36, 60]
    auc_values = []
    for t in time_points:
        idx = np.searchsorted(surv.index, t)
        if idx >= len(surv):
            idx = len(surv) - 1
        closest_time = surv.index[idx]
        y_true = (durations <= t) & (events == 1)
        y_score = 1 - surv.loc[closest_time].values
        try:
            auc_val = roc_auc_score(y_true, y_score)
            auc_values.append(auc_val)
        except ValueError:
            print(f"[{label}] No valid AUC at {t} months (all 0 or 1 labels)")
            auc_values.append(np.nan)

    print(f"[{label}] Time-dependent AUC: 1-year: {auc_values[0]:.3f}, 3-year: {auc_values[1]:.3f}, 5-year: {auc_values[2]:.3f}")

    # 返回用于绘图的数据
    return {
        'cindex': cindex,
        'auc_values': auc_values,
        'surv': surv,
        'durations': durations,
        'events': events
    }

# 分别评估训练集和测试集
results_train = evaluate_model(surv_train, durations_train, events_train, label="Train")
results_test = evaluate_model(surv_test, durations_test, events_test, label="Test")

# 绘制校准曲线（1年、3年、5年）
# def plot_calibration_curves(results, label=""):
#     durations = results['durations']
#     events = results['events']
#     surv = results['surv']

#     time_points = [12, 36, 60]
#     titles = ['1-year Calibration', '3-year Calibration', '5-year Calibration']

#     plt.figure(figsize=(15, 5))
#     for i, t in enumerate(time_points):
#         plt.subplot(1, 3, i + 1)
#         idx = np.searchsorted(surv.index, t)
#         if idx >= len(surv):
#             idx = len(surv) - 1
#         closest_time = surv.index[idx]

#         pred_prob = 1 - surv.loc[closest_time].values
#         actual = (durations <= t) & (events == 1)

#         prob_true, prob_pred = calibration_curve(actual, pred_prob, n_bins=10)

#         plt.plot(prob_pred, prob_true, 's-', label='Model')
#         plt.plot([0, 1], [0, 1], '--', label='Perfectly calibrated')
#         plt.title(titles[i])
#         plt.xlabel("Predicted probability")
#         plt.ylabel("Actual probability")
#         plt.legend()
#     plt.suptitle(f"Calibration Curves - {label}")
#     plt.tight_layout(rect=[0, 0, 1, 0.95])
#     plt.show()

def get_calibration_data(results, time_points):

    # 获取单个数据集（训练集/测试集）在指定时间点的校准曲线数据
    durations = results['durations']
    events = results['events']
    surv = results['surv']
    calib_data = []
    for t in time_points:
        idx = np.searchsorted(surv.index, t)
        if idx >= len(surv):
            idx = len(surv) - 1
        closest_time = surv.index[idx]
        pred_prob = 1 - surv.loc[closest_time].values
        actual = (durations <= t) & (events == 1)
        prob_true, prob_pred = calibration_curve(actual, pred_prob, n_bins=10)
        calib_data.append((prob_true, prob_pred))
    return calib_data

# 合并函数：绘制2行3列的校准曲线大图（上：训练集，下：测试集；左到右：1/3/5年）
def plot_combined_calibration_curves(results_train, results_test, time_points=[12, 36, 60]):
    calib_train = get_calibration_data(results_train, time_points)
    calib_test = get_calibration_data(results_test, time_points)
    titles = ['1-year', '3-year', '5-year']
    row_labels = ['Training Group', 'Validation Group']

    # 增大图尺寸 + 调整行列间距
    fig, axes = plt.subplots(2, 3, figsize=(20, 14))  # 增大画布尺寸
    fig.suptitle('Calibration Curves (1, 3, 5 Years)', fontsize=18, y=0.99)  # 调大总标题
    
    # 调整子图之间的间距（hspace=行间距，wspace=列间距）
    plt.subplots_adjust(
        left=0.05,  # 左边界
        right=0.95, # 右边界
        top=0.92,   # 上边界
        bottom=0.08,# 下边界
        hspace=0.25,# 行间距
        wspace=0.2  # 列间距
    )

    for row in range(2):
        for col in range(3):
            ax = axes[row, col]
            t = time_points[col]
            if row == 0:
                prob_true, prob_pred = calib_train[col]
                label = row_labels[0]
            else:
                prob_true, prob_pred = calib_test[col]
                label = row_labels[1]
            
            # 调小每个子图的标题和标签字体（避免拥挤）
            ax.plot(prob_pred, prob_true, 's-', label=label, color='steelblue', lw=2)
            ax.plot([0, 1], [0, 1], '--', label='Perfectly Calibrated', color='red', lw=1.5)
            ax.set_title(f'{titles[col]} Calibration ({label})', fontsize=13)  # 子图标题
            ax.set_xlabel("Predicted Probability", fontsize=12)  # X轴标签
            ax.set_ylabel("Actual Probability", fontsize=12)      # Y轴标签
            ax.legend(fontsize=12)  # 图例字体
            ax.grid(alpha=0.3)

    plt.show()

# 画校准曲线
# plot_calibration_curves(results_train, "Train")
# plot_calibration_curves(results_test, "Test")
plot_combined_calibration_curves(results_train, results_test)

# 绘制ROC曲线（1年、3年、5年）
def plot_roc_curves(results, label=""):
    durations = results['durations']
    events = results['events']
    surv = results['surv']

    time_points = [12, 36, 60]
    labels = ['1-year ROC', '3-year ROC', '5-year ROC']

    plt.figure(figsize=(10, 8))
    for t, label_roc in zip(time_points, labels):
        idx = np.searchsorted(surv.index, t)
        if idx >= len(surv):
            idx = len(surv) - 1
        closest_time = surv.index[idx]

        y_score = 1 - surv.loc[closest_time].values
        y_true = (durations <= t) & (events == 1)

        if len(np.unique(y_true)) < 2:
            print(f"[{label}] No valid ROC at {t} months (all 0 or 1 labels)")
            continue

        fpr, tpr, _ = roc_curve(y_true, y_score)
        roc_auc = auc(fpr, tpr)

        plt.plot(fpr, tpr, lw=2, label=f'{label_roc} (AUC = {roc_auc:.3f})')

    plt.plot([0, 1], [0, 1], 'k--', lw=1)
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title(f'ROC Curves at 1, 3, and 5 Years - {label}', fontsize = 13)
    plt.legend(loc="lower right", fontsize = 13)
    plt.grid()
    plt.tight_layout()
    plt.show()

# 画roc曲线
plot_roc_curves(results_train, "Training")
plot_roc_curves(results_test, "Validation")

# 画log-rank检验图
def plot_km_with_logrank(durations, events, group, title=""):
    kmf = KaplanMeierFitter()
    
    plt.figure(figsize=(8, 6))
    for group_val, label in [(0, 'Low Risk Group'), (1, 'High Risk Group')]:
        mask = (group == group_val)
        if mask.sum() > 0:
            kmf.fit(durations[mask], events[mask], label=label)
            kmf.plot_survival_function(ax=plt.gca(), ci_show=False)

    # Log-rank test
    mask_low = (group == 0)
    mask_high = (group == 1)
    results = logrank_test(durations[mask_low], durations[mask_high],
                           event_observed_A=events[mask_low], event_observed_B=events[mask_high])
    p_value = results.p_value

    plt.title(f"Kaplan-Meier Curve by Risk Group\n{title}" + 
              (f" (Log-rank p = {p_value:.1e})" if p_value < 0.001 else f" (Log-rank p = {p_value:.3f})"), fontsize = 13)
    plt.xlabel("Time", fontsize = 13)
    plt.ylabel("Survival Probability", fontsize = 13)
    plt.legend(title='Risk Group', fontsize = 13)
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# 获取训练集和测试集的风险得分
risk_train = model.predict(x_train).flatten()
risk_test = model.predict(x_test).flatten()

# 合并训练集+测试集的风险得分，计算统一的整体中位数
all_risk_scores = np.concatenate([risk_train, risk_test])  # 合并所有风险得分
unified_median_threshold = np.median(all_risk_scores)      # 计算整体中位数作为统一阈值
print(f"整体风险得分中位数（统一分层阈值）：{unified_median_threshold:.3f}")

# 统一阈值应用于训练集和测试集
# 训练集风险分组：> 统一中位数 = 高风险（1），≤ 统一中位数 = 低风险（0）
group_train = (risk_train > unified_median_threshold).astype(int)
# 测试集风险分组：使用与训练集完全相同的统一阈值
group_test = (risk_test > unified_median_threshold).astype(int)

# 打印分组样本数量，便于验证
print(f"训练集 - 低风险样本数：{np.sum(group_train == 0)}, 高风险样本数：{np.sum(group_train == 1)}")
print(f"测试集 - 低风险样本数：{np.sum(group_test == 0)}, 高风险样本数：{np.sum(group_test == 1)}")

# 绘制训练集和测试集的KM曲线（使用统一阈值分层）
plot_km_with_logrank(durations_train, events_train, group_train, title="Training Group")
plot_km_with_logrank(durations_test, events_test, group_test, title="Validation Group")