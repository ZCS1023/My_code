import pandas as pd
import requests
import json
from concurrent.futures import ThreadPoolExecutor
import re

# API 信息
API_URL = "https://api.coze.cn/open_api/v2/chat"
BOT_ID = "7479805290738122761"
API_KEY = "pat_auvBUJI72c18DxJdhz1YUEvr09BRfIVljfXU8hAruJ5Qd5codhUFmYWqhORuactT"

# 读取 Excel
df = pd.read_excel("dialog.xlsx", sheet_name="Sheet1")

# 清理并解析JSON字符串
def clean_and_parse_json(json_str):
    try:
        # 处理可能的单引号问题，将单引号替换成双引号
        json_str = re.sub(r"'", '"', json_str)
        # 处理可能的多余空格
        json_str = json_str.strip()
        data = json.loads(json_str)
        return data
    except json.JSONDecodeError as e:
        print(f"JSON解析错误: {e}，原字符串: {json_str}")
        return {}

# API 请求函数
def extract_labels_from_coze(dialog_text):
    headers = {
        "Authorization": f"Bearer {API_KEY}",
        "Content-Type": "application/json"
    }
    payload = {
        "bot_id": BOT_ID,
        "user": "test_user",  # 用户 ID，可以自定义
        "query": dialog_text,  # 用户输入的对话文本
        "stream": False  # 是否流式返回，False 表示一次性返回完整结果
    }

    try:
        response = requests.post(API_URL, headers=headers, json=payload, timeout=10)  # 设置超时时间
        response_data = response.json()
        print("完整的API响应:", response_data)  # 打印完整响应用于调试

        # 提取返回的 JSON 结果
        if "messages" in response_data:
            return response_data["messages"][0]["content"]  # 返回消息内容
        else:
            return "{}"
    except Exception as e:
        print(f"请求 Coze API 失败: {e}")
        return "{}"

# 使用多线程加速处理
with ThreadPoolExecutor(max_workers=10) as executor:  # 10个线程并发请求
    results = list(executor.map(extract_labels_from_coze, df["对话"]))

# 解析JSON字符串并填充数据
extracted_results = []
for result in results:
    extracted_result = clean_and_parse_json(result)
    extracted_results.append(extracted_result)

# 定义函数将字典转换为标准格式字符串
def format_dict(d):
    if isinstance(d, str):
        try:
            d = eval(d)
        except:
            return d
    lines = [f"{key}:{value}" for key, value in d.items()]
    return '\n'.join(lines)

# 对抽取结果应用格式化函数
formatted_results = [format_dict(result) for result in extracted_results]

df["标签抽取结果"] = formatted_results

# 保存处理后的 Excel
df.to_excel("抽取结果.xlsx", sheet_name="Sheet1", index=False)
print(f"处理完成，结果已保存至抽取结果.xlsx")
