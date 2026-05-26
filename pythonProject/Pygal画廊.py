from random import randint
class Die():
    def __init__(self, num_sides):
        self.num_sides = num_sides
    def roll(self):
        return randint(1, self.num_sides)

die1 = Die(8)
die2 = Die(10)
results = []
for roll_num in range(5000):
    result = die1.roll()+die2.roll()
    results.append(result)

frequencies=[]
for value in range(2,die1.num_sides+die2.num_sides+1):
    frequency=results.count(value)
    frequencies.append(frequency)

import pygal
hist=pygal.Bar()
hist.title='5000次结果图'
sum=die1.num_sides+die2.num_sides
hist.x_title='结果'
hist.x_labels=[]
for i in range(2,sum+1):
    hist.x_labels.append(i)
hist._y_title='Frequncy'

hist.add('D6',frequencies)
hist.render_to_file('die.svg')