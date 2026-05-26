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