# MSA-Kit Profile-Based Star Refinement（实验版 V4）

这是基于此前 MSA-Kit-V2 的独立试验版本，面向 DNA/RNA；保留原有 C++17、命令行、FASTA/FASTQ 输入输出与 EasyMSA `benchmark/build.sh`、`benchmark/run.sh` 接口。**这是研究原型，并非已有论文方法的完整复现，也不保证超过 MAFFT 等工具。**

## 新增算法

在 `src/align.cpp` 的 V2 初始自适应中心星比对和插入槽局部优化后，调用 `profile_iterative_refine(result)`。

1. 将已有MSA按3种确定性分组方式划为两个profile（偶数/奇数行、前半/后半行、每3行取1行）。
2. 删除各组内的全gap列，同时保留组内所有序列之间的既有列对应。
3. 使用全局动态规划重新对齐两个profile。DP的列-列得分等于两组所有序列对的匹配得分之和，匹配+2、错配-1、gap与碱基配对-2、N配对0，T/U视为同种碱基。
4. 仅当跨组sum-of-pairs目标**严格提高**时接受新MSA。组内列关系不变，因此全局相同SP目标也严格提高。
5. 最多两轮；为限制资源，仅在序列数3–64、MSA长度不超过1200列时尝试；单次DP最多700000个状态，超限会保留原版结果而不是运行无限搜索。

与V2只在中心插入槽移动gap不同，该版本允许重新排列两组profile的完整列对应关系；不使用量子退火和外部比对软件。

**局限：** 这是有限的确定性profile划分精修，而不是完整MAFFT渐进比对或任意子树搜索；内部SP不保证外部Q、TC提升。大型任务可能因为预算跳过精修。V2自身的全矩阵NW仍可能耗时/耗内存。不同工具的评分体系不一定一致。

## 编译和运行

在仓库根目录：

```bash
bash benchmark/build.sh
python3 benchmark/test_adapter.py
bash benchmark/run.sh examples/tiny.fasta profile_tiny.fasta 1 dna
cat profile_tiny.fasta
```

通用契约测试：

```bash
cmake -S . -B build-learning -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=ON
cmake --build build-learning --parallel 2
(cd build-learning && ctest -L contract --output-on-failure)
```

`reference` 分类中的某些断言要求旧星比对的特定gap布局，算法变化时应以 `contract` 和新方法的质量回归验证为准，不能把旧布局一致性当成新方法正确性的必要条件。

## EasyMSA 上传

上传仓库**完整文件**：尤其 `src/`、`include/`、`extern/`（包含 `spdlog/sinks/sink.h`）、`benchmark/`、`CMakeLists.txt` 和许可证。不要上传本地 `build-*` 产物。运行 `bash benchmark/build.sh` 验证独立构建，然后提交 GitHub 仓库链接至 EasyMSA 平台。V0、V1、V2、V4应在相同数据集与配置下比较Q/TC/SP、内存和耗时。

## 已完成的本地验证

2026-11-08：在容器内 GCC/CMake 重新构建、Contract 2/2、adapter 6/6、官方tiny FASTA均通过。80组固定种子合成数据与V2对照：新版56组输出不同，所有80组满足等长、去gap可恢复、没有全gap列且上述内部SP目标不下降。没有进行EasyMSA真实平台验证，亦不代表实际DNA数据的质量提升。
