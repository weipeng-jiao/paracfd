

## 二维流体模拟demo

## 编译

```
g++ -O2 -fopenmp -o test test.cpp
```

## SPH结果

![sph](C:\Users\HiWin10\Desktop\LAB\Development Library\cfd\paracfd\paracfd\paracfd\sph.png)

## 线程测试

![test](C:\Users\HiWin10\Desktop\LAB\Development Library\cfd\paracfd\paracfd\paracfd\test.png)

**测试平台Windows Intel(R)core(TM)i5-10500CPI@3.10GHz**

**物理线程6 逻辑线程12**

### 加速比分析

**加速比维持在6左右，在达到逻辑线程达到6（物理线程数量）时已经趋于稳定**

- **数据依赖未处理：内存访问次数太多**

- **分支预测未处理：过多的分支，影响缓存，使得超线程处理没有发挥性能**

## VOF结果

![VOF](C:\Users\HiWin10\Desktop\LAB\Development Library\cfd\paracfd\paracfd\paracfd\VOF.png)

**VOF程序未添加表面张力求解，所以形成上图的液面问题**:joy:





- [ ] **后期添加完整版VOF C++并行版本（debug中...）**:broken_heart:
- [ ] **sph程序性能调优中，目前初步改了一部分代码加速比上升到9左右，后期加入SSE指令上升会快一些**:sleeping:

