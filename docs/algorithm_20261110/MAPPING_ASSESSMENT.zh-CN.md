# minimap2 是否是最合适的定位器？

评估日期：2026-11-11。方法与源码审查，不是不同 mapper 的性能 benchmark；本次没有更换生产调用器或重跑全队列。

## 结论

**minimap2 是当前可用的基线，没有证据证明它在该任务上最优，也没有证据支持现在直接换掉它。** 应将“候选物理位置是否找全”和“给定两段序列怎样精确比对”分开。更精确的局部比对不能自动解决重复拷贝定位，最优编辑得分也不等于真实生物学对应。

当前应保留 minimap2，先评估候选保留策略及 graph occurrence 对来源判定的实际贡献；若改算法，优先考虑在图上下文约束下比较多个候选来源。后者是未来方向，不是当前实现已完成的功能。

## 当前依赖与信息损失

`source.rs` 调用 `minimap2 -cx asm20`，默认 k=19，并覆盖若干 preset 设置，包括 `-s 40`、`--end-bonus 100`、`--secondary=yes -N 20`。`context.rs` 实际只解析 `tp:A:P` 且 MAPQ≥20 的记录。secondary 参与 minimap2 内部主次判定，但其记录不作为 PathAllele 来源候选保留。

这是两个不同的限制：mapper 可能根本没有找到某位置；也可能输出了竞争位置，但当前 parser 丢弃了它。仅提高 `-N` 不会自动让这些候选进入 PathAllele。现有 graph-only 补充只能在唯一节点、方向、序列及间距条件下恢复部分未覆盖区域。

候选定位器接口目前耦合 minimap2 的 PAF、`tp:A:P`、MAPQ 和 `=/X` CIGAR 语义。替换程序名不等于兼容。不能把其他工具输出任意赋为 MAPQ=60、标记 primary，或把仅含 M 的 CIGAR 当作精确匹配。

## 可比较的路线

| 路线 | 能回答的问题 | 需要控制的限制 |
|---|---|---|
| minimap2 asm5/asm10/asm20 与当前参数对照 | 预设、种子和 gap 代价是否影响候选覆盖/错误配对 | 改变的参数须逐项记录；不能只按 VCF 行数或可调用率选优 |
| MUMmer4 / nucmer | 使用不同种子/聚类策略能否补回候选，尤其重复区域 | 默认种子偏唯一参考匹配；`--maxmatch` 增加候选，也可能增加歧义；delta 与当前 PAF 接口不同 |
| wfmash | 不同粗定位和 WFA 精比对组合的覆盖与来源一致性 | 仍有候选生成阶段；需要检查 PAF tags、精确 CIGAR 和评分含义；若 PGGB 本身用它构图，方法误差并不独立 |
| 已定位窗口内的 WFA2 / BiWFA | 边界已明确时，局部序列比对与表示是否稳定 | 本身不是全基因组位置搜索器，不是 minimap2 定位阶段的直接替代品 |
| 图 occurrence 先产生多个受限候选窗口，再局部比对 | 降低对单个 primary 选择的依赖 | 需要新的来源竞争、重复访问和不可辨识状态规则；可能组合爆炸；尚未实现/验证 |

minimap2 官方文档按分歧程度给出不同 assembly presets；`asm20` 面向数个百分点差异，不等于对每个 EBV 局部区间最优。[minimap2 手册](https://github.com/lh3/minimap2/blob/master/minimap2.1)

nucmer 的默认唯一参考匹配及可选 `--maxmatch` 见 [MUMmer4 文档](https://github.com/mummer4/mummer)。wfmash 将近似同源定位与 WFA 碱基比对分开，见 [wfmash 文档](https://github.com/waveygang/wfmash)。WFA2 的比对模式见 [WFA2 文档](https://github.com/Marco-Sola-Lab/WFA2-lib)。上表关于如何用于 PathAllele 的建议是接口和方法层面的推断，不是这些工具对 EBV 的已测性能结论。

## oriP 对本问题的启示

对 E5_sLCL 的固定结果做序列回放：PathAllele 的一条记录、vg/wave 的两条记录和 Linear 的 21 条记录都能重建相同 104-bp 序列。原始 vg snarl ALT 已含这条序列；保存的 Linear PAF 将该区域置于一个正向长比对中。这里的主要差别是局部反向结构是否被识别及如何表示，**不是 minimap2 丢掉了该 ALT**。不能据此断言更换 mapper 就能改善该位点；完整序列保留与通用 RC 分类同样重要。详见 [oriP 案例](ORIP_CASE_STUDY.zh-CN.md)。

## 建议的验证设计

1. 固定参考、原始组装、PGGB 图、样本集合、输入版本与 downstream caller 规则。
2. 先测现有 mapper 参数；再测候选保留/筛选；最后比较替代 mapper，避免把多个改动的收益混在一起。
3. 用独立已知序列/来源真值区分单拷贝、重复访问、近同源多拷贝、局部倒位附加替换、断裂 contig、低复杂度及确实不可辨识状态。已知 oriP 只作回归控制。
4. 分开报告正确来源召回与错配、完整 allele 正确性、结构分类、callability、错误参考调用、missingness 校准、时间和峰值内存。不可将更多 ALT 或更少 missing 直接视为更好。
5. 对 graph occurrence、候选筛选、共有锚点和 RC 分类分别做消融；用完整序列/统一事件单位比较，避免将 DEL+INS 和一个 replacement 的行数差异当作检出率差异。

候选更丰富仍不能解决所有信息上不可辨识的重复拷贝，正确结果有时就是保留不确定。
