num = [0,1,3,1,0,4]

non_zero = [a for a in num if a != 0]
zero_num = len(num) - len(non_zero)
finnal = non_zero + [0] * zero_num

print(finnal)