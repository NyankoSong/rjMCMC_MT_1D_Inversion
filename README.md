# rjMCMC_MT_1D_Inversion

通过可逆跳跃马尔科夫链蒙特卡洛方法实现一维大地电磁反演

| 文件名             | 类型     | 功能             | 文件名            | 类型     | 功能                 |
| ------------------ | -------- | ---------------- | ----------------- | -------- | -------------------- |
| main_TransD.m      | Script   | 可变维反演主程序 | forward_func.m    | Function | 一维正演             |
| TransD.m           | Function | 可变维反演函数   | likelihood_func.m | Function | 似然函数             |
| bostick_func.m     | Function | BOSTICK反演      | perturb_func.m    | Function | 扰动模型参数         |
| generate_cov.m     | Function | 生成协方差矩阵   | mesh_func.m       | Function | 将模型参数舍入到网格 |
| plot_chain.m       | Script   | 绘制各项图件     | test_model.m      | Script   | 合成测试观测数据     |
| plot_mesh.m        | Script   | 绘制各项图件     | calc_indicators.m | Script   | 计算各项反演结论指标 |
| plot_result_pics.m | Script   | 绘制各项图件     | Data_mats         | Folder   | 部分测试用数据       |
