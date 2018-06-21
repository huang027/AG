import os
os.getcwd()
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
data=pd.read_excel('C:\\Users\\guojia\\Desktop\\distant.xlsx',index_col=None,header=None)
coordinates=[]
for i in range(1,len(data)+1):
    coordinates.append(i)
distmat=np.array(data)
numant = 4  # 蚂蚁个数
numcity = len(coordinates)  # 城市个数
alpha = 1  # 信息素重要程度因子
beta = 5  # 启发函数重要程度因子
rho = 0.05  # 信息素的挥发速度
Q = 1
iter = 0
itermax = 250
etatable = 1.0 / (distmat + np.diag([1e10] * numcity))  # 启发函数矩阵，表示蚂蚁从城市i转移到矩阵j的期望程度
pheromonetable = np.ones((numcity, numcity))  # 信息素矩阵
pathtable = np.zeros((numant, numcity)).astype(int)  # 路径记录表
lengthaver = np.zeros(itermax)  # 各代路径的平均长度
lengthbest = np.zeros(itermax)  # 各代及其之前遇到的最佳路径长度
pathbest = np.zeros((itermax, numcity))  # 各代及其之前遇到的最佳路径长度
while iter < itermax:
    # 随机产生各个蚂蚁的起点城市
    if numant <= numcity:  # 城市数比蚂蚁数多
        pathtable[:, 0] = np.random.permutation(range(0, numcity))[:numant]
    else:  # 蚂蚁数比城市数多，需要补足
        pathtable[:numcity, 0] = np.random.permutation(range(0, numcity))[:]
        pathtable[numcity:, 0] = np.random.permutation(range(0, numcity))[:numant - numcity]
    length = np.zeros(numant)  # 计算各个蚂蚁的路径距离
    for i in range(numant):
        # i=0
        visiting = pathtable[i, 0]  # 当前所在的城市
        # visited = set() #已访问过的城市，防止重复
        # visited.add(visiting) #增加元素
        unvisited = set(range(numcity))  # 未访问的城市
        unvisited.remove(visiting)  # 删除元素
        for j in range(1, numcity):  # 循环numcity-1次，访问剩余的numcity-1个城市
            # j=1
            # 每次用轮盘法选择下一个要访问的城市
            listunvisited = list(unvisited)
            probtrans = np.zeros(len(listunvisited))
            for k in range(len(listunvisited)):
                probtrans[k] = np.power(pheromonetable[visiting][listunvisited[k]], alpha) \
                               * np.power(etatable[visiting][listunvisited[k]], alpha)
            cumsumprobtrans = (probtrans / sum(probtrans)).cumsum()
            cumsumprobtrans -= np.random.rand()
            k = listunvisited[np.where(cumsumprobtrans>0)[0][0]]  # 下一个要访问的城市
            pathtable[i, j] = k
            unvisited.remove(k)
            # visited.add(k)
            length[i] += distmat[visiting][k]
            visiting = k
        length[i] += distmat[visiting][pathtable[i, 0]]  # 蚂蚁的路径距离包括最后一个城市和第一个城市的距离
    # print length
    # 包含所有蚂蚁的一个迭代结束后，统计本次迭代的若干统计参数
    lengthaver[iter] = length.mean()
    if iter == 0:
        lengthbest[iter] = length.min()
        pathbest[iter] = pathtable[length.argmin()].copy()
    else:
        if length.min() > lengthbest[iter - 1]:
            lengthbest[iter] = lengthbest[iter - 1]
            pathbest[iter] = pathbest[iter - 1].copy()
        else:
            lengthbest[iter] = length.min()
            pathbest[iter] = pathtable[length.argmin()].copy()
         # 更新信息素
    changepheromonetable = np.zeros((numcity, numcity))
    for i in range(numant):
        for j in range(numcity - 1):
            changepheromonetable[pathtable[i, j]][pathtable[i, j + 1]] += Q / distmat[pathtable[i, j]][
                pathtable[i, j + 1]]
        changepheromonetable[pathtable[i, j + 1]][pathtable[i, 0]] += Q / distmat[pathtable[i, j + 1]][pathtable[i, 0]]
    pheromonetable = (1 - rho) * pheromonetable + changepheromonetable
    iter += 1  # 迭代次数指示器+1
    # 观察程序执行进度，该功能是非必须的

    if (iter - 1) % 20 == 0:
        print(iter - 1)
        # 做出平均路径长度和最优路径长度
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12, 10))
axes[0].plot(lengthaver, 'k', marker=u'')
axes[0].set_title('Average Length')
axes[0].set_xlabel(u'iteration')
axes[1].plot(lengthbest, 'k', marker=u'')
axes[1].set_title('Best Length')
axes[1].set_xlabel(u'iteration')
fig.savefig('Average_Best.png', dpi=500, bbox_inches='tight')
plt.close()
# 作出找到的最优路径图

bestpath = pathbest[-1]
'''
plt.plot(coordinates[:, 0], coordinates[:, 1], 'r.', marker=u'$\cdot$')
plt.xlim([-100, 2000])
plt.ylim([-100, 1500])
'''
for i in range(numcity - 1):  #
    m, n = bestpath[i], bestpath[i + 1]
    print(int(m), int(n))
print(pathbest)
'''
plt.plot([coordinates[int(m)][0], coordinates[int(n)][0]], [coordinates[int(m)][1], coordinates[int(n)][1]], 'k')
plt.plot([coordinates[int(bestpath[0])][0], coordinates[int(n)][0]], [coordinates[int(bestpath[0])][1], coordinates[int(n)][1]], 'b')
ax = plt.gca()
ax.set_title("Best Path")
ax.set_xlabel('X axis')
ax.set_ylabel('Y_axis')
plt.savefig('Best Path.png', dpi=500, bbox_inches='tight')
plt.close()
'''


