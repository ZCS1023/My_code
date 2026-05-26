'''class Dog():
    def __init__(self,name,age):
        self.name=name
        self.age=age
    def sit(self):
        print(self.name+" is sitting!")
    def roll_over(self):
        print(self.name+" is rolling over!")
x=Dog('Pussy',3)

print(x.name)
print(x.age)'''
import numpy as np

#描述餐馆信息
'''class restaurant():
    def __init__(self,restaurant_name,cuisine_type):
        self.name=restaurant_name
        self.type=cuisine_type
    def describe_restaurant(self):
        print("该餐馆名字为："+self.name)
        print("该餐馆菜品种类为："+self.type)
    def open_restaurant(self):
        print("该餐馆正在营业！")

x=restaurant("牛牛餐馆","川菜")
x.describe_restaurant()
x.open_restaurant()'''

#描述汽车信息
'''class Car():
    def __init__(self,name,year,model):
        self.name=name
        self.year=year
        self.model=model
        self.odometer_reading=0
    def describe_car(self):
        long_name=self.name+' '+self.model+' '+str(self.year)
        print("汽车的详细信息为："+long_name)
    def upgrade_odometer(self,mile):
        if mile>=self.odometer_reading:
            self.odometer_reading=mile
        else:
            print("禁止将里程数回调！")
    def read_odometer(self):
        print("这辆车已经跑了"+str(self.odometer_reading)+"公里了")

#继承
class Electric(Car):
    def __init__(self,name,year,model):
        super().__init__(name,year,model)
        self.battery=100
    def describe_battery(self):
        print("目前电池容量为："+str(self.battery))

my_car=Electric('BMW','2022','X5')
my_car.describe_car()
my_car.describe_battery()'''

#将实例用作属性
'''class Car():
    def __init__(self,name,year,model):
        self.name=name
        self.year=year
        self.model=model
        self.odometer_reading=0
    def describe_car(self):
        long_name=self.name+' '+self.model+' '+str(self.year)
        print("汽车的详细信息为："+long_name)
    def upgrade_odometer(self,mile):
        if mile>=self.odometer_reading:
            self.odometer_reading=mile
        else:
            print("禁止将里程数回调！")
    def read_odometer(self):
        print("这辆车已经跑了"+str(self.odometer_reading)+"公里了")

class Battery():
    def __init__(self,size=75):
        self.size=size
    def describe_size(self):
        print("这辆车的电容量是："+str(self.size))
    def get_range(self):
        if self.size==70:
            range=240
            print("这辆车续航里程大约为："+str(range))
        else:
            range=280
            print("这辆车满电续航里程大约为："+str(range))
#继承
class Electric(Car):
    def __init__(self,name,year,model):
        super().__init__(name,year,model)
        self.battery=Battery()
        
my_car=Electric('BMW','2022','X5')
my_car.describe_car()
my_car.battery.describe_size()
my_car.battery.get_range()'''

'''class Foods():
    def __init__(self,name,fds):
        self.fds=fds
        self.name=name
        self.number=len(fds)

    def show_foods(self):
        print("小店"+self.name+"现有美食：")
        for i in self.fds:
            print(i)
        print("欢迎前来品尝！")

cuisine=Foods("食来食往",["蛋烘糕","片片鱼","钵钵鸡","翘脚牛肉"])
cuisine.show_foods()'''

import numpy as np
class Matrix():
    def __init__(self,row,column,elements1,elements2):
        self.row=row
        self.column=column
        self.elements1=elements1
        self.elements2=elements2

    def Matrix_sum(self):
        result = np.zeros((self.row,self.column))

        for i in range(self.row):
            for j in range(self.column):
                result[i][j]=self.elements1[i][j]+self.elements2[i][j]
        for r in result:
            print(r)

    def Marrix_minus(self):
        result = np.zeros((self.row, self.column))

        for i in range(self.row):
            for j in range(self.column):
                result[i][j] = self.elements1[i][j] - self.elements2[i][j]
        for r in result:
            print(r)

    def Matrix_multiply(self):
        result = np.zeros((self.row, self.column))
        result = self.elements1 * self.elements2
        print(result)

    def Submatrix(self,i,j):
        m = self.row
        n = self.column
        C = [[self.elements1[x][y] for y in range(n) if y != j] for x in range(m) if x != i]
        return C

    def Det(self):
        m = self.row
        n = self.column
        value = 0
        for j in range(n):
            value += (((-1) ** (j + 2)) * self.elements1[0][j] *
                      self.Det(self.submatrix(self.elements1, 0, j)))
            return value

    def inverseofmatrix(self):
        m = self.row
        n = self.column
        C = np.zeros((self.row, self.column))
        d = self.Det(self.elements1)
        for i in range(m):
            for j in range(n):
                C[i][j] = ((-1) ** (i + j + 2)) * self.Det(self.Submatrix(self.elements1, j, i))
                C[i][j] = C[i][j] / d
        return C

matrix1=Matrix(2,3,[[5,6,7],[1,2,3]],[[7,8,5],[2,3,4]])
print("矩阵A+矩阵B得到：")
matrix1.Matrix_sum()
print("\n")
print("矩阵A-矩阵B得到：")
matrix1.Marrix_minus()
matrix2=Matrix(2,2,np.array([[3,2],[2,4]]),np.array([[5,3],[1,3]]))
print("\n")
print("矩阵C×矩阵D得到：")
matrix2.Matrix_multiply()
print(matrix2.inverseofmatrix())