# rjMCMC_MT_1D_Inversion

通过可逆跳跃马尔科夫链蒙特卡洛方法实现一维大地电磁反演

测试模型生成：test_model.m

主执行文件：main_TranD.m

主函数：TransD.m

* Bostick反演：bostick_func.m

* 一维正演：forward_func.m

* 生成协方差矩阵：generate_cov.m

* 似然函数：likelihood_func.m

* 建议函数：proposal_func.m

* 使扰动数据点吸附到数据网格：mesh_func.m

变量迭代过程变化曲线绘制：plot_chain.m

绘制结果图：plot_mesh.m

