def judge_square(stick, target, index, sides):
    if index == len(stick):
        return all(side == target for side in sides)
    for i in range(4):
        if sides[i] + stick[index] > target:
            continue
        sides[i]  += stick[index]
        if judge_square(stick, target, index + 1, sides):
            return True
        sides[i] -= stick[index]
    return False
num = int(input())
for _ in range(num):
    parts = list(map(int, input().split()))
    n = parts[0]
    stick = parts[1:]
    lengh = sum(stick)
    if n < 4 or lengh % 4 != 0:
        print('no')
        continue
    target = lengh // 4
    stick.sort(reverse = True)
    if stick[0] > target:
        print('no')
        continue
    sides = [0] *4
    if judge_square(stick, target, 0, sides):
        print('yes')
    else:
        print('no')