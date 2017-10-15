#!/usr/bin/python3

lines1 = open('lr-train','r').readlines()
model1 = open('lr-model','r').readlines()

lines2 = open('lr-train-strict','r').readlines()
model2 = open('lr-model-strict','r').readlines()

means1 = [
sum([float(line.split()[i]) for line in lines1]) / len(lines1)
for i in range(len(lines1[0].split()))
]
means2 = [
sum([float(line.split()[i]) for line in lines2]) / len(lines2)
for i in range(len(lines2[0].split()))
]

coef1 = list(map(float, model1[1].split()))
inter1 = float(model1[2])
k1 = coef1 + [inter1]

coef2 = list(map(float, model2[1].split()))
inter2 = float(model2[2])
k2 = coef2 + [inter2]

for i in range(len(k1)):
    print('%3d %10.7lf %10.7lf      %10.7lf' % (i+1, k1[i], means1[i], k1[i]*means1[i]))
for i in range(len(k2)):
    print('%3d %10.7lf %10.7lf      %10.7lf' % (i+1, k2[i], means2[i], k2[i]*means2[i]))


for i in range(len(k2)):
    print('$k_{%d}$ & %10.7lf & %10.7lf & %10.7lf & %10.7lf  \\\\' % (i+1, k1[i], k2[i], k1[i]*means1[i], k2[i]*means2[i]))
