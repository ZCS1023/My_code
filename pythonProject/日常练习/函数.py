#等效函数的调用
'''def city(city_name,city_country='China'):
    print("City "+city_name+" belongs to country "+city_country)

city(city_name='Beijing')
city(city_name='Chengdu')
city(city_name='Los Angels',city_country='USA')'''

#返回字典（动态选择是否返回年龄）
'''def person_name(first_name,last_name,age=''):
    person={
        'first':first_name,
        'last':last_name,
    }
    if age:
        person['age']=age
    return person
musician=person_name('James','Lebron')
print(musician)'''

#while与函数的结合
'''def tell_name(first_name,last_name):
    full_name=first_name+last_name
    return full_name
active=True
while active:
    print("你好，请告诉我你的名字：")
    print("如果想退出该程序，请输入quit")
    f_name=input("你的姓是：")
    if f_name=='quit':
        break
    l_name=input("你的名是：")
    if l_name=='quit':
        break
    full_name=tell_name(f_name,l_name)
    print("你的名字是："+full_name+"!\n")'''

#函数输出所输入的内容
'''def greet_user(user_name):
    for user in user_name:
        print("Hello, "+user.title()+"!")

print("How many users?")
x=int(input())
print("What are these users' names?")
user_names=[]
for i in range(x):
    user_names.append(input())
greet_user(user_names)'''

#魔术师传奇
'''def show_magicians(magician_name):
    for name in magician_name:
        print("你好, "+name.title()+"魔术师!")

magicians_names=['郑策书','王飞','丁超龙','王一凡']
show_magicians(magicians_names)

def make_great():
    for i in range(4):
        magicians_names[i]="The Great"+magicians_names[i]

make_great()
show_magicians(magicians_names)'''

'''names=['黄焖鸡','牛肉面','麻辣烫']
prices=[12,11,8.5]
dictionary=dict(zip(names,prices))
print(dictionary)'''

#使用任意数量的形参和实参
'''def make_pizza(size,*toppings):
    print("Making a "+str(size)+"-inch pizza with the following toppings:")
    for topping in toppings:
        print("-"+topping)'''
