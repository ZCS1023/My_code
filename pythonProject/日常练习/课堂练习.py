#输入一个数判断是否为素数且判断孪生素数
'''import math
def judge(number):
    k=int(math.sqrt(number))
    for i in range(2,k+1):
        if number%i==0:
            return ("不是一个素数")
            break
    else:
            return ("是一个素数")
while True:
    print("请输入一个正整数，输入1来退出程序:")
    x=int(input())
    if x==1:
        break
    else:
        result=judge(x)
    print(result)
    if result=='是一个素数':
        result1=judge(x+2)
        if result1=='不是一个素数':
            print(str(x)+"和"+str(x+2)+"不是孪生素数\n")
        else:
            print(str(x)+"和"+str(x+2)+"是孪生素数\n")'''

#找出2-1000的孪生素数
'''import math
def judge(number):
    k=int(math.sqrt(number))
    for i in range(2,k+1):
        if number%i==0:
            return ("不是一个素数")
            break
    else:
            return ("是一个素数")

for i in range(2,1001):
    if judge(i) and judge(i+2)=='是一个素数':
        print(str(i)+"和"+str(i+2)+"是卵生素数")
    else:
        continue'''

'''import time
def check(names):
    attendances={}
    for name in names:
        print(name)
        att=int(input("到了吗？"))
        time.sleep(0.5)
        if att==1:
            attendances[name]='到课'
        elif att==2:
            attendances[name]='请假'
        else:
            attendances[name]='旷课'
    return attendances

print("How many students?")
x=int(input())
print("What are these students' names?")
students_names=[]
for i in range(x):
    students_names.append(input())
print("考勤结果为：(1到课，2请假，其余数字旷课)")
result=check(students_names)
print(result)'''

'''from random import randint
class Die():
    def __init__(self,num_sides):
        self.num_sides=num_sides

    def roll(self):
        self.result=randint(1,self.num_sides)

die=Die(6)
results=[]
for roll_result in range(1000):
    die.roll()
    results.append(die.result)

frequencies=[]
for value in range(1,die.num_sides+1):
    frequency=results.count(value)
    frequencies.append(frequency)

print(frequencies)'''

'''import numpy as np
class Matrix():
    def __init__(self, row, column, elements1, elements2):
        self.row = row
        self.column = column
        self.elements1 = elements1
        self.elements2 = elements2

    def Matrix_sum(self):
        result = np.zeros((self.row, self.column))

        result=self.elements1+self.elements2
        print(result)


matrix = Matrix(2, 3, np.array([[5, 6, 7], [1, 2, 3]]), np.array([[7, 8, 5], [2, 3, 4]]))
matrix.Matrix_sum()'''

import copy
# 本函数求矩阵的逆
def matrix_ni(matrix):
    extend_matrix = copy.deepcopy(matrix)
    l = len(matrix)

    for i in range(0, l):    #在矩阵右边补充一个单位矩阵，使用初等变换求逆矩阵
        extend_matrix[i].extend([0]*i)
        extend_matrix[i].extend([1])
        extend_matrix[i].extend([0]*(l-i-1))
    for i in range(0, len(extend_matrix)):    #判断矩阵对角线上是否有0，有0则置换，如置换不了，则没有逆矩阵
        if extend_matrix[i][i] == 0:
            for j in range(i, len(extend_matrix)):
                if extend_matrix[j][i] != 0:
                    extend_matrix[i], extend_matrix[j] = extend_matrix[j], extend_matrix[i]
                    break
            if j >= len(extend_matrix):
                print('没有逆矩阵')
                return 0
            break
    for i in range(0, len(extend_matrix)):    #开始计算逆矩阵
        f = extend_matrix[i][i]
        for j in range(0, len(extend_matrix[i])):    #先把行换为0
            extend_matrix[i][j] /= f
        for m in range(0, len(extend_matrix)):
            if m == i:
                continue
            b = extend_matrix[m][i]
            for n in range(0, len(extend_matrix[i])):    #在把列换为0
                extend_matrix[m][n] -= extend_matrix[i][n] * b
    for i in range(0, len(extend_matrix)):
        extend_matrix[i] = extend_matrix[i][l:]
    return extend_matrix

A=[[3,2],[1,4]]
print(matrix_ni(A))

