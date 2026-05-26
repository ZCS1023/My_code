print("请输入积分下限：")
a=int(input())
print("请输入积分上限：")
b=int(input())
N=1000
value=0
deltax=(b-a)/N
for i in range(0,N):
    ia=a+i*deltax
    ib=ia+deltax
    area=(ia**2+ib**2)/2*deltax
    value+=area
print("积分区域面积为：\n"+str(value))
