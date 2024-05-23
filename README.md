# Road-MultiObjective-Optimization-GeneticImmuneAlgorithm

## 1 问题描述

在工程建设中，投资者往往希望用最短的时间，最少的 投入达到最佳的工程质量、安全等要求。如何通过均衡各个目标使总目标最优这 就是工程项目多目标优化问题。

该项目为一简单的办公楼，包括 14 道工序，需要建立进度（工期T）-成本C-安全水平S-质量Q模型。

本仓库内含同一模型、两份数据，使用遗传免疫算法在MATLAB下完成。

**模型**见文件 优化模型及参数.pdf，**第一份数据**见 实例数据1.xlsx，**第二份数据**见 实例数据2.xlsx。

此外还进行了遗传免疫算法和基本遗传算法对比，使用两种测试函数绘制迭代曲线。注意：基本遗传算法只相对于遗传免疫算法修改了 交叉 部分；迭代曲线波动范围大。



## 2 文件说明

main1.m 为求解road 10道工序实例模型（数据1）

main2.m 为求解road 14道工序实例模型（数据2）

draw.m 为绘制四种数据对比图（未用）

test.m 为对比算法-免疫遗传算法，绘制迭代曲线并标注出最大值点

drawfunction.m 为绘制两种测试函数并标注出最大值点



### 3 参考结果

模型求解结果：



![解1 (1)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%A7%A31%20(1).png?raw=true)

![解1 (2)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%A7%A31%20(2).png?raw=true)

![解1 (3)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%A7%A31%20(3).png?raw=true)

![解1 (4)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%A7%A31%20(4).png?raw=true)



算法对比结果：

![测试函数 (1)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E6%B5%8B%E8%AF%95%E5%87%BD%E6%95%B0%20(1).png?raw=true)

![测试函数 (2)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E6%B5%8B%E8%AF%95%E5%87%BD%E6%95%B0%20(2).png?raw=true)

![迭代曲线 (1)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%BF%AD%E4%BB%A3%E6%9B%B2%E7%BA%BF%20(1).png?raw=true)

![迭代曲线 (2)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%BF%AD%E4%BB%A3%E6%9B%B2%E7%BA%BF%20(2).png?raw=true)

![迭代曲线 (3)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%BF%AD%E4%BB%A3%E6%9B%B2%E7%BA%BF%20(3).png?raw=true)

![迭代曲线 (4)](https://github.com/kunzi051/Road-MultiObjective-Optimization-GeneticImmuneAlgorithm/blob/main/img/%E8%BF%AD%E4%BB%A3%E6%9B%B2%E7%BA%BF%20(4).png?raw=true)
