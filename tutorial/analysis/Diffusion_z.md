# 计算扩散系数沿某方向的变化

该脚本首先对每个粒子计算msd, 随后根据其z位置进行统计, 最终每个bin内的msd信息会输出到本地的msd.log中, 扩散系数通过对msd进行线性拟合得到, 最终每个bin的扩散系数会输出到filename.dz中.
使用前首先打开Diffusion_z.py脚本设置好msd计算方式:
```python
#----input-------------------------------------------------------------------------
unit_t=0.1 #设置轨迹中每两帧之间的时间间隔, 单位:ps
T=5        #统计时长,即最终msd-t图中的t的范围, 单位:ps
dt=0.1     #统计间隔,每个dt更新一次原点, 单位:ps
ave_time=5  #统计范围, 每次更新后对ave_time进行平均, 单位:ps
```
运行方式: python3 Diffusion_z.py filename nbins direction
```bash
python3 Diffusion_z.py md 50 2 #0,1,2分别对应x,y,z
```

# fortran code
