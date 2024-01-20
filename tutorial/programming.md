# 如何在本地浏览器上打开远程服务器上的jupyter-notebook
1. 在服务器上输入jupyter-notebook
2. 在本地terminal上输入：
   ssh -N -f -L localhost:8888:localhost:8888 username@IP
3. 使用本地浏览器打开远程服务器上jupyter-notebook的链接即可
