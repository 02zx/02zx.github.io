# 如何在本地浏览器上打开远程服务器上的jupyter-notebook
1. 在服务器上输入jupyter-notebook
2. 在本地terminal上输入：
   ssh -N -f -L localhost:8888:localhost:8888 username@IP
3. 使用本地浏览器打开远程服务器上jupyter-notebook的链接即可



# 如何在VMD中添加字幕


trace variable vmd_frame(1) w sdf
mol new
proc sdf {args} {global vmd_frame
graphics 0 delete all
graphics 0 color 16
graphics 0 text {0 0 0} "[expr $vmd_frame(1)*1]ps" size 6 thickness 4}

#color 16是黑色
