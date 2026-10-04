# T1 二分图最大权匹配 本地测评

本目录是自建的本地测评程序，**不参与正式测评**，仅用于在提交前核对正确性与
性能。测评方式与公开的 `testing/T0_ab` 保持一致：只收集
`problems/T1_bipartite_matching/src/` 下的全部 `.cpp`，用本目录 `include/`
中的官方接口副本编译，不使用题目自带的 `CMakeLists.txt` 与 `main.cpp`。

用例构成：

| 档位 | 用例 | 核对方式 |
| --- | --- | --- |
| 基础 | 10 个边界用例（空图、无边、单条负边、零权边、题面样例、全负边等） | 位压缩 DP 暴力 + 参考实现 |
| 对拍 | 400 组小规模随机图（含负权、零权、混合权） | 位压缩 DP 暴力 |
| 拓展 1 | 5 组 `L,R <= 200` 稠密图 | 独立 SPFA 费用流参考实现 |
| 拓展 2 | 5 组 `L,R <= 500` 稠密图 | 独立 SPFA 费用流参考实现 |
| 拓展 3 | 5 组 `L,R <= 2000`、`E <= 50000` 稀疏图 | 独立 SPFA 费用流参考实现 |

运行方式：

```bash
cmake -S testing/T1_bipartite_matching -B /tmp/pip27-t1-grader
cmake --build /tmp/pip27-t1-grader -j
/tmp/pip27-t1-grader/grader          # 可选参数：随机种子
```
