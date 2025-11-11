总体说明：

1.TDL文件夹下为基于张量字典学习的能谱CT重建的程序。

2.运行前配置：在运行程序前，需要进行一些配置。首先，在本文件夹的Depandency子文件夹中有poblano_toolbox和tensor_toolbox两个工具箱，是Sandia国家实验室发布的张量工具箱。按顺序安装这两个工具箱，在各自文件夹中有安装说明。本程序还会调用到上层路径下的DL_func文件夹，该文件夹是在Ron Rubinstein编写的字典学习工具箱的基础上修改得到的，包含了向量字典和张量字典学习的诸多函数，由于涉及到了与C语言的混编，因此需要在Matlab中通过mex -setup命令来设置混编。

3.ReconTDL_main.m是本方法的主函数,ReconTDL_core.m是方法的程序核心，重建（中间）结果将保存在上层路径下的\Results文件夹下创建的新文件夹中。DL_CPD_ratio.m用于训练张量字典的主函数，kcpd20140705.m是训练张量字典的程序核心，对MOBY仿真数据和小鼠实际数据训练得到的张量字典已经保存在上层路径下的\Results文件夹下。imratioMC.m对多通道图像的各个通道模值进行归一化。

4.程序执行：一般情况下，事先对多通道数据进行FBP重建，然后使用DL_CPD_ratio.m根据重建图像训练得到张量字典并保存。在ReconTDL_main.m中加载投影数据和训练得到的字典，则可重建图像。在本程序中，已经训练并保存了MOBY仿真数据和实际小老鼠数据的张量字典，因此直接运行ReconTDL_main.m即可，大概需要运行8个小时。


General Description:

1.TDL Folder: This folder contains programs for spectral CT reconstruction based on Tensor Dictionary Learning (TDL).

2.Pre-run Configuration:
Before running the program, some configuration steps are required.
In the Dependency subfolder, there are two toolboxes — poblano_toolbox and tensor_toolbox, both released by Sandia National Laboratories. Install these two toolboxes in order; installation instructions can be found in their respective folders.
The program also calls functions from the DL_func folder in the parent directory. This folder is a modified version of Ron Rubinstein’s dictionary learning toolbox and includes many functions for vector and tensor dictionary learning.
Since it involves mixed programming with C language, you need to configure the MEX compiler in MATLAB using the command:
mex -setup

3.ReconTDL_main.m — the main function of this method.
ReconTDL_core.m — the core implementation of the reconstruction algorithm.
The (intermediate) reconstruction results will be saved in a newly created folder under the parent directory’s \Results folder.
DL_CPD_ratio.m — the main function for training tensor dictionaries.
kcpd20140705.m — the core program for tensor dictionary training.
Tensor dictionaries trained on MOBY simulation data and real mouse data have already been saved under the \Results folder in the parent directory.
imratioMC.m — used for normalizing the magnitude values across multiple image channels.

4.Program Execution:
Typically, multi-channel data is first reconstructed using FBP (Filtered Back Projection). Then, DL_CPD_ratio.m is used to train the tensor dictionary based on the reconstructed images, which is subsequently saved.
By loading the projection data and the trained dictionary into ReconTDL_main.m, image reconstruction can be performed.
In this program, tensor dictionaries trained on both MOBY simulation data and real mouse data have already been provided, so you can directly run ReconTDL_main.m. The total runtime is approximately 8 hours.

