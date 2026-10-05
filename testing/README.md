# 本地测评

本目录下的 `T1_bipartite_matching/`、`T2_dynamic_kth/` 是自建的本地测评程序，
**不参与正式测评**，仅用于提交前核对正确性与性能。它们与官方公开的
`T0_ab/` 测评示例保持同样的组织方式：

- 只收集 `problems/<题号>/src/` 下的全部 `.cpp`，与正式测评一致；
- 用各目录 `include/` 中的官方接口头文件副本编译，考生无法影响测评接口；
- 不使用题目自带的 `CMakeLists.txt` 与 `main.cpp`。

各题的用例构成、核对方式与运行命令见对应目录下的 `README.md`。

运行全部本地测评：

```bash
cmake -S testing/T1_bipartite_matching -B /tmp/pip27-t1-grader
cmake --build /tmp/pip27-t1-grader -j && /tmp/pip27-t1-grader/grader

cmake -S testing/T2_dynamic_kth -B /tmp/pip27-t2-grader
cmake --build /tmp/pip27-t2-grader -j && /tmp/pip27-t2-grader/grader
```

说明：`testing/` 目录会被仓库 `.gitignore` 中的 `Testing/` 规则在大小写不敏感
的文件系统上误伤，新增文件需 `git add -f` 强制加入版本控制。
