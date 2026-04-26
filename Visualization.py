import numpy as np
import matplotlib.pyplot as plt

with open("data.txt") as file:
    M, M_vis, N, K = map(int, file.readline().split())
    x, y, x_segment, y_segment = [], [], [], []
    k = 0
    for i in range(M):
        xi, yi = map(float, file.readline().split())
        x.append(xi)
        y.append(yi)
        if (i == (N-1)*k) & (k <= K):
            k += 1
            x_segment.append(xi)
            y_segment.append(yi)
    dot, func, lagr = [], [], []
    for i in range(M_vis):
        di, fi, li = map(float, file.readline().split())
        dot.append(di)
        func.append(fi)
        lagr.append(li)


plt.figure(figsize=(7, 7))
plt.plot(dot, func, 'b-', label='f(x)')
plt.plot(dot, lagr, 'r--', label='L(x)')
plt.scatter(x, y, c = 'g', marker = 'o')
plt.scatter(x_segment, y_segment, c = 'b', marker = 'o')
ax = plt.gca()  
ax.spines['left'].set_position('zero')   
ax.spines['bottom'].set_position('zero') 
ax.spines['right'].set_color('none')     
ax.spines['top'].set_color('none')

plt.grid(True)
plt.legend()
plt.show()
