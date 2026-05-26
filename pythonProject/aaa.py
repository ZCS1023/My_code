for a in range(1,17):
    for b in range(1,17):
        for c in range(1,17):
            if a+b+c==16 and (a-1)*(c-1)==(b-1)**2:
                print(a,b,c)