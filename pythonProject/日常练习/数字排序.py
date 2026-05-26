print("请输入数字1：")
num1 = int(input())
print("请输入数字2：")
num2 = int(input())
print("请输入数字3：")
num3 = int(input())
if num1 > num2:
    temp = num1
    num1 = num2
    num2 = temp
if num1 > num3:
    temp = num1
    num1 = num3
    num3 = temp
if num2 > num3:
    temp = num2
    num2 = num3
    num3 = temp
print(num1, num2, num3)
