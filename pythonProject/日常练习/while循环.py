#while选择何时退出
'''prompt="I will repeat it to you:\n"
prompt+="Enter 'quit' to end the programe.\n"
message=""
while message!='quit':
    message=input(prompt)
    if message!='quit':
        print(message)'''

#while使用标志
'''prompt="I will repeat it to you:\n"
prompt+="Enter 'quit' to end the programe.\n"
message=""
active=True
while active:
    message=input(prompt)
    if message=='quit':
        active=False
    else:
        print(message)'''

#列表之间的元素移动
'''unconfirmed_users=['Ming','Hong','Ding']
confirmed_users=[]
while unconfirmed_users:
    current_user=unconfirmed_users.pop()
    print("Verifying user:"+current_user)
    confirmed_users.append(current_user)
print("\n以下人员已确认：")
for i in range(len(confirmed_users)):
    print(confirmed_users[i])'''

#while删除特定元素
'''pets=['dog','cat','dragon','bird','cat','rabbit']
while 'cat' in pets:
    pets.remove('cat')
print(pets)'''

#使用用户输入填充字典
responses={}
active=True
while active:
    name=input("你的名字是？\n")
    response=input("你想去爬哪座山？\n")

    responses[name]=response

    repeat=input("你还想让其他人回答吗？请回答yes或no\n")
    if repeat=='no':
        active=False
print("------结果------")
for name,response in responses.items():
    print(name+"想去爬"+response)