import pandas as pd
import numpy as np

df = pd.read_csv(r'C:\Users\zheng\Documents\code\PycharmProjects\jd_user_behavior.csv')
df['datetime'] = pd.to_datetime(df['timestamp'],unit='s')


df['date'] = df['datetime'].dt.date
df['hour'] = df['datetime'].dt.hour

print(df.head())

# 按date列分组统计个数
daily_count = df.groupby('date').size()

action_count = df.groupby('action').size()
