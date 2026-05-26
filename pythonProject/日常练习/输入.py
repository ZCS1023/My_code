#基本输入
'''age=input("What's your age?\n")
print("Your age is "+age+"!")
if int(age)>=20:
    print("YES")
else:
    print("NO")'''

#判断输入数字奇偶性
'''number=input("请输入所要判断的数字：\n")
print("您输入的数字为:"+number+"!")
number=int(number)
if number%2==0:
    print("这是一个偶数！")
else:
    print("这是一个奇数！")'''

#判断闰年
'''year=input("请输入一个年份：\n")
print("你输入的年份是："+str(year)+"!")
year=int(year)
if year%4==0 and year%100!=0 or year%400==0:
    print("这是一个闰年！")
else:
    print("这不是一个闰年！")'''

#交税问题
'''salary=input("请输入薪水：\n")
salary=int(salary)
y=0
sum=0
print("你输入的薪水为："+str(salary)+"!")
if salary<5000:
    print("恭喜您，无需交税！")
elif 5000<=salary<8000:
    y=(salary-5000)*0.03
elif 8000<=salary<12000:
    y=(salary-8000)*0.05+90
elif 12000<=salary<20000:
    y=(salary-12000)*0.1+290
elif 20000<=salary<100000:
    y=(salary-20000)*0.15+1090
else:
    y=(salary-100000)*0.25+13090
print("您一共需要交税："+'%.2f'%y+"元！")
sum=salary-y
print("您交完税后所得个人工资为："+str(sum)+"元!")'''
