#matplotlib函数初体验
'''import matplotlib.pyplot as plt
input_values=[1,2,3,4,5]
squares=[]
for i in input_values:
    square=i*2+1
    squares.append(square)
plt.plot(input_values,squares,linewidth=2)
plt.title("S-N",fontsize=24)
plt.xlabel("Value",fontsize=14)
plt.ylabel("S-of-V",fontsize=14)
plt.tick_params(axis='both',labelsize=10)

plt.show()'''

#颜色映射
'''import matplotlib.pyplot as plt
x_values=list(range(1,1001))
y_values=[x**2 for x in x_values]

plt.scatter(x_values,y_values,c=y_values,cmap=plt.cm.Blues,edgecolor='none',s=40)
plt.title("S-N",fontsize=24)
plt.xlabel("Value",fontsize=14)
plt.ylabel("S-of-V",fontsize=14)
plt.axes([0,1100,0,1100000])
plt.show()'''

#立方数
'''import matplotlib.pyplot as plt
x_values=list(range(1,5001))
y_values=[x**3 for x in x_values]
plt.title("Function Figure")
plt.scatter(x_values,y_values,c=y_values,cmap=plt.cm.Blues,edgecolor='none',s=40)
plt.axes([0,5100,0,130000000000])
plt.show()'''

#随机漫步(离散点版)
'''from random import choice
class RandomWalk():
    def __init__(self,num_points):
        #设置初始值
        self.num_points=num_points
        self.x_points=[0]
        self.y_points=[0]

    def fill_walk(self):
        while len(self.x_points)<self.num_points:
            #x控制左右方向
            x_direction=choice([-1,1])
            x_distance=choice([1,2,3,4,5])
            x_step=x_direction*x_distance

            #y控制上下方向
            y_direction=choice([-1,1])
            y_distance=choice([1,2,3,4,5])
            y_step=y_direction*y_distance

            #计算下一个点的坐标位置
            next_x=self.x_points[-1]+x_step
            next_y=self.y_points[-1]+y_step
            self.x_points.append(next_x)
            self.y_points.append(next_y)

import matplotlib.pyplot as plt
while True:
    x=RandomWalk(50000)
    x.fill_walk()
    number=list(range(x.num_points))
    plt.scatter(x.x_points,x.y_points,c=number,cmap=plt.cm.Blues,edgecolors='none',s=5)

    #突出起点和终点
    plt.scatter(0,0,c='green',edgecolors='none',s=50)
    plt.scatter(x.x_points[-1],x.y_points[-1],c='red',edgecolors='none',s=50)
    plt.show()

    print("想再画一次吗？(输入yes或no)")
    y=input()
    if y=='no':
        break'''

#模拟花粉运动
'''from random import choice
class RandomWalk():
    def __init__(self,num_points):
        #设置初始值
        self.num_points=num_points
        self.x_points=[0]
        self.y_points=[0]

    def fill_walk(self):
        while len(self.x_points)<self.num_points:
            #x控制左右方向
            x_direction=choice([-1,1])
            x_distance=choice([1,2,3,4,5])
            x_step=x_direction*x_distance

            #y控制上下方向
            y_direction=choice([-1,1])
            y_distance=choice([1,2,3,4,5])
            y_step=y_direction*y_distance

            #计算下一个点的坐标位置
            next_x=self.x_points[-1]+x_step
            next_y=self.y_points[-1]+y_step
            self.x_points.append(next_x)
            self.y_points.append(next_y)

import matplotlib.pyplot as plt
while True:
    x=RandomWalk(5000)
    x.fill_walk()
    number=list(range(x.num_points))
    plt.plot(x.x_points,x.y_points,linewidth=2)

    #突出起点和终点
    plt.scatter(0,0,c='green',edgecolors='none',s=50)
    plt.scatter(x.x_points[-1],x.y_points[-1],c='red',edgecolors='none',s=50)
    plt.show()
print("想再画一次吗？(输入yes或no)")
y = input()
if y == 'no':
    break
    '''