# 如何在本地浏览器上打开远程服务器上的jupyter-notebook
1. 在服务器上输入jupyter-notebook
2. 在本地terminal上输入：
   ssh -N -f -L localhost:8888:localhost:8888 username@IP
3. 使用本地浏览器打开远程服务器上jupyter-notebook的链接即可



# GROMACS编译
加载环境：
module purge
source /public/software/compiler/intel/2024.0.1/setvars.sh
module add apps/python/3.8.15/python-3.8.15 compiler/gcc/11.3.0 compiler/cmake/3.25.1 compiler/nvhpc/23.7
编译安装：
mkdir build && cd build
<!-- for CPU -->
<!-- cmake .. -DCMAKE_INSTALL_PREFIX=/public/software/apps/gromacs-2023.5-ultrasound -DCMAKE_C_COMPILER=icx -DCMAKE_CXX_COMPILER=icpx -DGMX_DOUBLE=off   -DGMX_MPI=off  -DGMX_SIMD=AVX_512  -DGMX_USE_LMFIT=none -DGMX_FFT_LIBRARY=fftpack -DGMX_OPENMP=on -->

<!-- for GPU -->
cmake .. -DCMAKE_INSTALL_PREFIX=/public/software/apps/gromacs-2023.5-ultrasound-30   \
-DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_DOUBLE=off   \
-DGMX_FFT_LIBRARY=fftpack -DGMX_MPI=off -DGMX_OPENMP=on   -DGMX_GPU=on \
-DGMX_SIMD=AVX_512 -DGMX_GPU=CUDA \
-DCUDA_TOOLKIT_ROOT_DIR=/public/software/compiler/nvhpc/Linux_x86_64/23.7/cuda  \
-DGMX_USE_LMFIT=none -DGMX_CUDA_TARGET_SM=86 -DGMX_CUDA_TARGET_COMPUTE=86 \
-DCUDA_cufft_LIBRARY=/public/software/compiler/nvhpc/Linux_x86_64/23.7/math_libs/lib64/libcufft.so

make -j 16 && make install
