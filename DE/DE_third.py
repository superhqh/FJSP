# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 10:21:26 2019

@author: 大卫
"""
import numpy as np
import matplotlib.pyplot as plt
import time
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
    
def initpopulation():
    pop = np.zeros((popsize,total_process*2))
    v = np.zeros((popsize,total_process*2))
    fitness = np.zeros(popsize)
    time_and_machindex = []
    for i in range(total_process):
        time_and_machindex.append([[int(contents[i][j]), j+1] for j in range(machine) if contents[i][j] != "-"])
    for i in range(global_size):
        global_init(pop[i], time_and_machindex)
        np.random.shuffle(pop[i][:total_process])
        fitness[i] = calculate(pop[i])
    for i in range(global_size, global_size+local_size):
        local_init(pop[i], time_and_machindex)
        np.random.shuffle(pop[i][:total_process])
        fitness[i] = calculate(pop[i])
    return pop, v, fitness

def global_init(x, time_and_machindex):
    machload = np.zeros(machine)
    job_list = [i for i in range(workpiece)]
    cur_begin = 0
    while True:
        if not job_list:
            break
        job_index = np.random.randint(0, len(job_list))
        for i in range(process):
            x[cur_begin+i] = job_list[job_index]+1
        process_list = [i for i in range(process)]
        while True:
            if not process_list:
                break
            process_index = np.random.randint(0, len(process_list))
            cur = time_and_machindex[job_list[job_index]*process+process_list[process_index]]
            temp_machload = np.zeros(machine)
            for i in cur:
                temp_machload[i[1]-1] = i[0] + machload[i[1]-1]
            min_load = np.inf
            for i in range(0, machine):
                if temp_machload[i] != 0 and temp_machload[i] < min_load:
                    min_load = temp_machload[i]
            min_load_list = [i+1 for i in range(len(temp_machload)) if temp_machload[i] == min_load]
            min_time = np.inf
            for i in min_load_list:
                for j in cur:
                    if j[1] == i and j[0] < min_time:
                        min_time = j[0]
                        break
            min_time_list = [i[1] for i in cur if i[0] == min_time]
            x[total_process+job_list[job_index]*process+process_list[process_index]] = min_time_list[np.random.randint(0, len(min_time_list))]
            machload[int(x[total_process+job_list[job_index]*process+process_list[process_index]]-1)] += min_time
            del process_list[process_index]
        cur_begin += process
        del job_list[job_index]

def local_init(x, time_and_machindex):
    machload = np.zeros(machine)
    cur_begin = 0
    for p in range(workpiece):
        for i in range(process):
            x[cur_begin+i] = p+1
        process_list = [i for i in range(process)]
        while True:
            if not process_list:
                break
            process_index = np.random.randint(0, len(process_list))
            cur = time_and_machindex[p*process+process_list[process_index]]
            temp_machload = np.zeros(machine)
            for i in cur:
                temp_machload[i[1]-1] = i[0] + machload[i[1]-1]
            min_load = np.inf
            for i in range(0, machine):
                if temp_machload[i] != 0 and temp_machload[i] < min_load:
                    min_load = temp_machload[i]
            min_load_list = [i+1 for i in range(len(temp_machload)) if temp_machload[i] == min_load]
            min_time = np.inf
            for i in min_load_list:
                for j in cur:
                    if j[1] == i and j[0] < min_time:
                        min_time = j[0]
                        break
            min_time_list = [i[1] for i in cur if i[0] == min_time]
            x[total_process+p*process+process_list[process_index]] = min_time_list[np.random.randint(0, len(min_time_list))]
            machload[int(x[total_process+p*process+process_list[process_index]]-1)] += min_time
            del process_list[process_index]
        cur_begin += process
        for i in range(len(machload)):
            machload[i] = 0

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

def func(x1, x2):
    #输入：两条父染色体，输出：一条子染色体
    #POX交叉
    seq = [i+1 for i in range(workpiece)]
    random_length1 = np.random.randint(2, len(seq)-1)
    set1 = set()
    for i in range(random_length1):
        index = np.random.randint(0, len(seq))
        set1.add(seq[index])
        seq.pop(index)
    set2 = set(seq)
    child1 = np.copy(x1)
    child2 = np.copy(x2)
    remain1 = [i for i in range(total_process) if x1[i] in set2]
    remain2 = [i for i in range(total_process) if x2[i] in set2]
    cursor1, cursor2 = 0, 0
    for i in range(total_process):
        if x2[i] in set2:
            child1[remain1[cursor1]] = x2[i]
            cursor1 += 1
        if x1[i] in set2:
            child2[remain2[cursor2]] = x1[i]
            cursor2 += 1
    if calculate(child1) < calculate(child2):
        child = child1
    else:
        child = child2
    return child

def mutation():
    #对种群中的每条染色体进行变异操作
    mid_population = np.zeros((popsize,total_process*2))
    for i in range(popsize):
        inside, outside = 0, 0
        if np.random.rand() > F:
            inside = np.copy(pop[i])
        else:
            inside = func(pop[i], gbestpop)
        index = np.random.randint(0, popsize)
        while index == i or not (pop[index]-pbestpop).any():
            index = np.random.randint(0, popsize)
        if np.random.rand() > F:
            outside = inside
        else:
            outside = func(inside, pop[index])
        mid_population[i] = outside
    return mid_population

def cross_and_select(mid_population):
    for i in range(popsize):
        individual = np.copy(pop[i])
        for j in range(total_process):
            if np.random.rand() <= Cr:
                individual[j+total_process] = mid_population[i][j+total_process]
            else:
                individual[j+total_process] = pop[i][j+total_process]
        array = handle(individual)
        for j in range(total_process):
            row = (array[j][0]-1)*process+array[j][1]-1
            if contents[row][int(individual[total_process+row])-1] == "-":
                individual[total_process+row] = clean_contents[row][len(clean_contents[row])-1][1]
        if calculate(individual) < calculate(pop[i]):
            pop[i] = individual

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
    popsize = 50 #种群规模
    global_size = popsize//2 #全局初始化的数目
    local_size = popsize - global_size #局部初始化的数目
    maxgen = 500 #最大迭代次数
    rangepop = (1,6) #粒子编码中机器选择部分的范围
    F = 0.1 #变异率
    Cr = 0.1 #交叉率
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
    
    pop,v,fitness = initpopulation()
    gbestpop,gbestfitness,pbestpop,pbestfitness = getinitbest(fitness,pop)
    
    DE_base = []
    iter_process = []
    begin = time.time()
    for i in range(maxgen):
        mid_population = mutation()
        cross_and_select(mid_population)
        
        iter_process.append(fitness.min())
        DE_base.append(gbestfitness)
        for j in range(popsize):
            fitness[j] = calculate(pop[j])
            if fitness[j] < pbestfitness[j]:
                pbestfitness[j] = fitness[j]
                pbestpop[j] = pop[j].copy()
            
        if fitness.min() < gbestfitness:
            gbestfitness = fitness.min()
            gbestpop = pop[fitness.argmin()].copy()

    print("按照基于全局与局部初始化策略的DE算法求得的最好的最大完工时间：",min(DE_base))
    print("按照基于全局与局部初始化策略的DE算法求得的最好的工艺方案：",gbestpop)
    end = time.time()
    print("整个迭代过程所用时间：{:.2f}s".format(end-begin))
    fig = plt.figure(figsize=(10,5))
    ax1 = fig.add_subplot(121)
    ax1.set_title("全局最优解的变化情况")
    ax1.plot(DE_base)
    ax2 = fig.add_subplot(122)
    ax2.set_title("每次迭代后种群适应度最小值的变化情况")
    ax2.plot(iter_process)
    plt.show()
