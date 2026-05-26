# import json
# import csv
#
# # 文件路径
# json_file = "350道题相似值.json"
# csv_file = "processed_output.csv"
#
# # 读取 JSON 文件
# with open(json_file, "r", encoding="utf-8") as file:
#     data = json.load(file)
#
# # 提取数据并处理
# processed_data = []
# for main_id, related_items in data["result"].items():
#     for related_id, xiangsizhi in related_items.items():
#         # 数据清洗：替换逗号为下划线，移除换行符
#         clean_main_id = str(main_id).replace(",", "_").replace("\n", "")
#         clean_related_id = str(related_id).replace(",", "_").replace("\n", "")
#         processed_data.append([clean_main_id, clean_related_id, float(xiangsizhi)])
#
# # 按照 xiangsizhi 倒序排列
# processed_data.sort(key=lambda x: x[2], reverse=True)
#
# # 写入 CSV 文件
# with open(csv_file, "w", newline="", encoding="utf-8") as file:
#     writer = csv.writer(file)
#     # 写入表头
#     writer.writerow(["id", "related_id", "xiangsizhi"])
#     # 写入数据
#     writer.writerows(processed_data)
#
# print(f"数据处理完成，结果已保存到 {csv_file}")

# import pandas as pd
# from collections import defaultdict
#
# # 假设您有以下文件路径（替换为实际文件路径）
# files = ["北京版.xlsx", "沪教版.xlsx", "冀教版.xlsx", "北师大版.xlsx"]
# versions = ["北京版", "沪教版", "冀教版", "北师大版"]
#
# # 合并数据
# combined_data = pd.DataFrame()
# for file, version in zip(files, versions):
#     # 读取文件
#     data = pd.read_excel(file)
#     data['version'] = version  # 添加版本信息
#     combined_data = pd.concat([combined_data, data], ignore_index=True)
#
# # 去重并统计重复项
# id_counts = defaultdict(list)
# for _, row in combined_data.iterrows():
#     id_counts[row['题目ID']].append(row['version'])
#
# # 构建去重后的数据和重复统计
# unique_ids = []
# repeated_stats = []
# for question_id, versions_list in id_counts.items():
#     unique_ids.append(question_id)
#     repeated_stats.append({"题目ID": question_id,
#                            "重复次数": len(versions_list),
#                            "重复版本": ",".join(set(versions_list))})
#
# # 去重后的数据
# unique_data = combined_data.drop_duplicates(subset=['题目ID']).reset_index(drop=True)
#
# # 重复统计表
# repeated_stats_df = pd.DataFrame(repeated_stats)
#
# # 保存到一个 Excel 文件中的不同工作表
# with pd.ExcelWriter("处理结果.xlsx", engine="openpyxl") as writer:
#     unique_data.to_excel(writer, sheet_name="去重后的数据", index=False)
#     repeated_stats_df.to_excel(writer, sheet_name="重复统计", index=False)
#
# print("数据处理完成，结果已保存到 '处理结果.xlsx'")

import pandas as pd

# 读取 Excel 文件
df = pd.read_excel('测试题.xlsx')

# 准备一个空的 DataFrame 来存储结果
result_df = pd.DataFrame(columns=['章节名称', '知识点名称'])

# 遍历每一行数据
for index, row in df.iterrows():
    chapter_name = row['章节名称']

    # 拆分知识点名称，假设知识点是用换行符 \n 分隔的
    knowledge_points = row['知识点名称'].split('\n')

    # 将拆分后的知识点名称与章节名称对应
    new_rows = []
    for point in knowledge_points:
        if point:  # 确保不添加空字符串
            new_rows.append({'章节名称': chapter_name, '知识点名称': point.strip()})  # 去除知识点两侧空格

    # 使用 pd.concat 来添加新行
    new_df = pd.DataFrame(new_rows)
    result_df = pd.concat([result_df, new_df], ignore_index=True)

# 保存到新的 Excel 文件
result_df.to_excel('处理后的测试题.xlsx', index=False)

print("处理完成，新的文件已保存为 '处理后的测试题.xlsx'")
