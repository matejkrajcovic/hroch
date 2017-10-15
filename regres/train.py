#!/usr/bin/python

from random import randint, shuffle
from math import exp
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import PolynomialFeatures
from sys import stdin, argv

orig,cross1,cross2,op = 'lr-train','lr-train-strict','lr-train-special','lr-model'
if len(argv) > 1 and argv[1] == 'strict':
    cross1,orig,cross2,op = 'lr-train','lr-train-strict','lr-train-special','lr-model-strict'
if len(argv) > 1 and argv[1] == 'special':
    cross2,cross1,orig,op = 'lr-train','lr-train-strict','lr-train-special','lr-model-special'


lines = open(orig,'r').readlines()
linesc1 = open(cross1,'r').readlines()
linesc2 = open(cross2,'r').readlines()
output = open(op,'w')
shuffle(lines)
shuffle(linesc1)
shuffle(linesc2)
linesc0 = lines[:5000]
linesc1 = linesc1[:5000]
linesc2 = linesc2[:5000]
lines = lines[5000:]

X = []
y = []
A, B = 0, -1

for line in lines:
    data = map(float, line.split())
    X.append(data[A:B])
    y.append(data[-1])

model = LogisticRegression(tol=0.0001)
model = model.fit(X,y)

coef = model.coef_[0]

print >>output, len(coef)
print >>output, ' '.join(map(str,coef))
print >>output, model.intercept_[0]

print(model.score(X,y))
print(sum(y)/float(len(y)))
print(coef)
print(model.intercept_)

def predict_on(linesc, name):
    X = []
    y = []
    for line in linesc:
        data = map(float, line.split())
        X.append(data[A:B])
        y.append(data[-1])
    predict = model.predict_proba(X)
    samew = sum([(y[i]-predict[i][1])**2 for i in range(len(y))])
    same = sum([(y[i]==(predict[i][1]>0.5)) for i in range(len(y))])
    print name, 1.*samew/len(y), 1.*same/len(y) 

predict_on(lines, orig)
predict_on(linesc0, orig)
predict_on(linesc1, cross1)
predict_on(linesc2, cross2)
