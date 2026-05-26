'''favorite_flower={
    'James':'yueji',
    'Irving':'meigui',
    'Curry':'mudan',
}
friends=['Irving']
for name,value in favorite_flower.items():
    print(name)
    if name in friends:
        print("Hi "+name+",I……is "+value+"!")
if 'durant' not in favorite_flower.keys():
    print("durant,please……")'''

#字典嵌套
'''alien_0={'color':'red','points':5}
alien_1={'color':'green','points':2}
alien_2={'color':'blue','points':3}
aliens=[alien_0,alien_1,alien_2]
for alien in aliens:
    print(alien)'''

#字典列表
'''aliens=[]
for alien_number in range(30):
    new_alien={'color':'red','speed':'slow','points':5}
    aliens.append(new_alien)
for alien in aliens[1:3]:
    if alien['color']=='red':
        alien['color']='green'
        alien['speed']='fast'
        alien['points']=10
for i in range(5):
    print(aliens[i])
print("一共创建了"+str(len(aliens))+"个外星人")'''

#字典中存列表
'''pizza={
    'crust':'thick',
    'topping':['mushrooms','extra cheese'],
}
print("You ordered a "+pizza['crust']+"-crust pizza "+"with toppings:")
for topping in pizza['topping']:
    print(topping)'''

#字典中存储字典
'''users={
    'james':{
        'first':'Irving',
        'last':'Kyrie',
        'location':'USA'
    },
    'lebron':{
      'first':'wang',
        'last':'yi',
        'location':'China',
    },
}
for user,user_info in users.items():
    print("Username:"+user)
    full_name=user_info['first']+user_info['last']
    print(full_name)
    location=user_info['location']
    print("Location:"+user_info['location']+"\n")'''

#城市信息
'''cities={
    'Chengdu':{
        'country':'China',
        'population':'1600w',
        'fact':'bashi',
    },
    'London':{
        'country':'England',
        'population':'500w',
        'fact':'hitorical',
    },
}
for name,country_info in cities.items():
    print("\n")
    print("Name:"+name.title())
    print("所属国家："+country_info['country'])
    print("人口："+country_info['population'])
    print("事实："+country_info['fact'])'''
