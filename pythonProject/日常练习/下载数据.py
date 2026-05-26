#导入csv数据并读取文件头
'''import pandas as pd
df=pd.read_csv(r'sitka_weather_07-2014.csv')
print("读取某行某列指定单元格:",df.iloc[0,0])
for index,column_header in enumerate(header_row):
        print(index,column_header)'''
import pandas as pd

#提取并读取数据
'''import csv

filename='sitka_weather_07-2014.csv'
with open(filename) as f:
    reader=csv.reader(f)
    header_row=next(reader)

    highs=[]
    for row in reader:
        high=int(row[1])
        highs.append(high)
    print(highs)'''

'''a=["张三","李四","王二"]
sr=pd.Series(a)
index=sr[sr.values=="李四"].index
print(index)'''

import pandas as pd
score=pd.read_excel(r'C11.xlsx')
print(score[['卷面成绩','作业成绩','实验成绩']])
score['总成绩']=score['卷面成绩']*0.6+score['作业成绩']*0.2+score['实验成绩']*0.2
score.to_excel('C11.xlsx',index=False)
print(score[['卷面成绩','作业成绩','实验成绩','总成绩']])

print("\n")
print("全班学生的总成绩均值为：")
full_average=0
for i in range(0,110):
    full_average=sum(score['总成绩'].values)/110
print(full_average)

score['离差']=score['总成绩']-full_average
score.to_excel('C11.xlsx',index=False)
print(score[['卷面成绩','作业成绩','实验成绩','总成绩','离差']])

from matplotlib import pyplot as plt
plt.hist(score['总成绩'].values,10,edgecolor='black')
plt.show()