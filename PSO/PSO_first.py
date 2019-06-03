# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 15:05:23 2019

@author: david
"""

import numpy as np
import matplotlib.pyplot as plt
import time

'''
改变问题规模时需要改动的参数：工件参数workpiece，工序数目process，机器数目machine，
机器选择部分的范围上限，迭代次数（规模小的时候可能降低迭代次数）
'''

#读取原始数据
contents = []
with open("data_first.txt") as f:
    string = f.readlines()
    for item in string:
        contents.append(item.strip().split(" "))
        
#对工序部分进行处理
def handle(x):
    #输入：粒子的位置，输出：对工序部分处理后的列表
    piece_mark = np.zeros(workpiece) #统计工序的标志
    array = [] #经过处理后的工序列表
    for i in range(total_process):
        piece_mark[int(x[i]-1)] += 1
        array.append((int(x[i]), int(piece_mark[int(x[i]-1)])))
    return array
    
def initpopvfit():
    pop = np.zeros((popsize,total_process*2))
    v = np.zeros((popsize,total_process*2))
    fitness = np.zeros(popsize)

    for i in range(popsize):
        #初始化工序部分
        for j in range(workpiece):
            for p in range(process):
                pop[i][j*process + p] = j+1
        np.random.shuffle(pop[i][:total_process])
        
        #初始化机器部分
        for j in range(total_process):
            index = np.random.randint(0, machine)
            while contents[j][index] == "-":
                index = np.random.randint(0, machine)
            pop[i][j+total_process] = index+1
                
        #计算各粒子初始的适应度 
        fitness[i] = calculate(pop[i])
    return pop,v,fitness

def calculate(x):
    # 输入:粒子位置，输出:粒子适应度值
    Tm = np.zeros(machine) #每个机器上的完工时间
    Te = np.zeros((workpiece, process)) #每个工序的完成时间
    array = handle(x) #经过处理后的工序部分

    for i in range(total_process):
        machine_index = int(x[total_process+(array[i][0]-1)*process+(array[i][1]-1)])-1 #contents数组中的纵坐标
        process_index = (array[i][0]-1)*process + (array[i][1]-1) #contents数组中的横坐标
        process_time = int(contents[process_index][machine_index])
        if array[i][1] == 1:
            Tm[machine_index] += process_time
            Te[array[i][0]-1][array[i][1]-1] = Tm[machine_index]
        else:
            Tm[machine_index] = max(Te[array[i][0]-1][array[i][1]-2], Tm[machine_index]) + process_time
            Te[array[i][0]-1][array[i][1]-1] = Tm[machine_index]
    return max(Tm)

def getinitbest(fitness,pop):
    # 群体最优的粒子位置及其适应度值
    gbestpop,gbestfitness = pop[fitness.argmin()].copy(),fitness.min()
    #个体最优的粒子位置及其适应度值,使用copy()使得对pop的改变不影响pbestpop，pbestfitness类似
    pbestpop,pbestfitness = pop.copy(),fitness.copy()
    return gbestpop,gbestfitness,pbestpop,pbestfitness

if __name__ == "__main__":  
    workpiece = 10 #工件数目
    process = 5 #每个工件的工序数目
    total_process = workpiece * process #工序的总数
    machine = 6 #机器数目
    maxgen = 500 #最大迭代次数
    w = 0.9 #惯性权重
    lr = (2,2) #加速因子 
    popsize = 50 #种群规模
    rangepop = (1,6) #粒子编码中机器选择部分的范围
    clean_contents = []
    for i in range(total_process):
        clean_contents.append([[int(contents[i][j]), j+1] for j in range(machine) if contents[i][j] != "-"])
        temp_sum = 0
        for j in range(len(clean_contents[i])):
            temp_sum += 1/clean_contents[i][j][0]
        for j in range(len(clean_contents[i])):
            clean_contents[i][j][0] = (1/clean_contents[i][j][0])/temp_sum
        clean_contents[i].sort()
        cumulation = 0
        for j in range(len(clean_contents[i])):
            cumulation += clean_contents[i][j][0]
            clean_contents[i][j][0] = cumulation
    
    pop,v,fitness = initpopvfit()
    gbestpop,gbestfitness,pbestpop,pbestfitness = getinitbest(fitness,pop)
    
    iter_process = np.zeros(maxgen)
    pso_base = np.zeros(maxgen)
    
    begin = time.time()
    for i in range(maxgen):
        #t=0.5
        #速度更新
        for j in range(popsize):
            v[j] = w*v[j]+lr[0]*np.random.rand()*(pbestpop[j]-pop[j])+lr[1]*np.random.rand()*(gbestpop-pop[j])

        #粒子位置更新
        #工序部分
        for j in range(popsize):
            store = []
            before = pop[j][:total_process].copy()
            pop[j] += v[j]
            reference = v[j][:total_process].copy()
            for p in range(total_process):
                store.append((reference[p], before[p]))
            store.sort()
            for p in range(total_process):
                pop[j][p] = store[p][1]
        
        pop = np.ceil(pop) 
        #机器部分
        for j in range(popsize):
            array = handle(pop[j])
            for p in range(total_process):
                if (pop[j][total_process+(array[p][0]-1)*process+(array[p][1]-1)] < rangepop[0] or pop[j][total_process+(array[p][0]-1)*process+(array[p][1]-1)] > rangepop[1]) \
                or (contents[(array[p][0]-1)*process+(array[p][1]-1)][int(pop[j][total_process+(array[p][0]-1)*process+(array[p][1]-1)]-1)] == "-"):
                    row = (array[p][0]-1)*process+(array[p][1]-1)
                    pop[j][total_process+(array[p][0]-1)*process+(array[p][1]-1)] = clean_contents[row][len(clean_contents[row])-1][1]

        iter_process[i] = fitness.min()
        pso_base[i] = gbestfitness
        #适应度更新
        for j in range(popsize):
            fitness[j] = calculate(pop[j])

        for j in range(popsize):
            if fitness[j] < pbestfitness[j]:
                pbestfitness[j] = fitness[j]
                pbestpop[j] = pop[j].copy()
                
        if pbestfitness.min() < gbestfitness :
            gbestfitness = pbestfitness.min()
            gbestpop = pop[pbestfitness.argmin()].copy()

    print("按照完全随机初始化的pso算法求得的最好的最大完工时间：",min(pso_base))
    print("按照完全随机初始化的pso算法求得的最好的工艺方案：",gbestpop)
    end = time.time()
    print("整个迭代过程所耗用的时间：{:.2f}s".format(end-begin))
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121)
    ax1.set_title("全局最优解的变化情况")
    ax1.plot(pso_base)
    ax2 = fig.add_subplot(122)
    ax2.set_title("每次迭代后种群适应度最小值的变化情况")
    ax2.plot(iter_process)
    plt.show()    
        