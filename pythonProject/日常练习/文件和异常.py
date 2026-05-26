#简单的文件操作
'''with open('digits.txt') as file_object:
    lines=file_object.readlines()
    for line in lines:
        print(line)'''
import pandas

#使用文件的内容
'''with open('digits.txt') as file:
    lines=file.readlines()
pi_string=''
for line in lines:
    pi_string+=line.strip()

print(pi_string)
print(len(pi_string))'''

#文件练习
'''filename='study experience.txt'
with open(filename) as file:
    contens=file.read()
    print(contens)
    lines=file.readlines()
study_string=""
study_string1=""
for line in lines:
    study_string+=line
print(study_string)
study_string1=study_string.replace('python','c')
print(study_string1)'''

#写入多行内容并尝试追加内容
'''filename='programme'
with open(filename,'a') as file:
    file.write(str(22)+"\n")
    file.write('I also like playing games.\n')'''

#关于编程的调查
'''def reason(y):
    file.write(y+"\n")
    return file
filename='reason'
with open(filename,'a') as file:
    while True:
        print("您愿意接受调查吗？（请输入是或否）")
        x=input()
        if x=='否':
            break
        else:
            print("请输入原因：")
            y=input()
            reason(y)'''

#try-except初尝试
'''try:
    print(5/0)
except ZeroDivisionError:
    print("You can't do that!")'''

#try-excep实战
'''print("请输入两个数字：（输入q来退出程序）")
while True:
    x=input("请输入第一个数字:\n")
    if x=='q':
        break
    y=input("请输入第二个数字：\n")
    if y=='q':
        break
    try:
        answer=int(x)/int(y)
    except ZeroDivisionError:
        print("出错！")
    else:
        print("计算结果为："+str(answer))'''

#json.dump函数存入文件
'''import json
x=input("What's your name?\n")
filename='number'
with open(filename,'w') as file:
    json.dump(x,file)'''

#json.load函数读取文件
'''import json
filename='number'
with open(filename) as file:
    name=json.load(file)
    print("你好，"+name)'''

#json练习
'''import json
filename='users_name'
try:
    with open(filename) as file:
        username=json.load(file)
        print("欢迎回来"+username)
except FileNotFoundError:
    username=input("What's your name?\n")
    with open(filename,'w') as file:
        json.dump(username,file)
        print("我记住你了"+username)'''

#在指定行位置追加文件内容
filename="poem.txt"
with open(filename,'r+',encoding='utf-8') as file:
    content=file.read()
    print(content)

with open('record.txt','w') as f:
    f.write('今天是2023年11月17日，星期五，天气晴\n')
    f.write("现在开始上课！\n")

with open('record.txt','r+') as fp:
    lines=[]
    for line in fp:
        line=line.strip()
        lines.append(line)
    fp.close()
lines.insert(1,content)
s="\n".join(lines)
fp=open('record.txt','w')
fp.write(s)
fp.close()

#读取文件
'''import pandas as pd
df=pd.read_excel('aa.xlsx',sheet_name='Sheet1',index_col=0)
print(df.head(2))
print(df.describe())'''