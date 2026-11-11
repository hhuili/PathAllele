# 方法定位与表述边界

更新：2026-11-11。当前实现基线 `057d45a9`。

推荐英文：

> PathAllele combines oriented graph-path occurrences with linear-reference mapping to extract complete alleles and conservatively project variants from haploid assembly pangenomes.

针对当前应用：

> Developed primarily for small-genome PGGB graphs, PathAllele distinguishes repeated visits to a graph node from the node's identity and uses mapping and path context to resolve their reference correspondence.

这里的 occurrence 是某节点在某条路径上的一次有方向的访问，由 path、step 和序列偏移区别。PGGB 图可以有环，路径可能因重复序列等结构重访节点，但不是每个图/路径都必然如此。循环拓扑不自动证明生物学重复拷贝。序列重复也可能分散在不同节点，不能只检查重复节点 ID。

## 与其他方法的区别

| 方法 | 组织变异的起始单位 | 边界 |
|---|---|---|
| 直接解析 pairwise alignment | 序列比对中的差异/比对块 | 具体 caller 也可以保留完整 allele、回放和 missingness |
| PathAllele | 由定位与 occurrence 证据限定的完整来源区间 | 仍依赖单线性参考及候选定位，不提供全图嵌套坐标系统 |
| vg deconstruct | snarl/site 及实际路径穿越 | 同样使用参考和样本路径；支持参考循环的上下文分配及嵌套导出 |
| pantree | reference tree 与非树 variant edges | 定义可覆盖线性参考之外的序列差异，与平面参考投影目标不同 |

完整 allele、序列重建、保守 missingness、使用样本路径及支持参考循环，都不是可以单独宣称为 PathAllele 独占的新颖性。当前可辩护的工具特点是 **occurrence-aware source resolution 与完整 allele / event projection 的组合**。是否形成新的方法贡献以及相对性能优势，需要 prior-art review、独立真值和消融。

## 应用范围

当前开发重点是 EBV 等小型、单倍体组装泛基因组；其他病毒是预期用途而非已经完成验证的所有物种。没有硬性“小基因组长度阈值”，也没有大型真核图的可扩展性证据。支持有环图不代表支持环形参考坐标；当前后者不支持。

## 主图与文档分工

主图使用 8-bp 倒位加邻近替换的短合成模拟，画出 mapping/occurrence 两路证据、可信侧翼、完整 allele 和结构分类/事件投影。图形借鉴 oriP 图的节点与曲线，但不复用真实位点或宣称此例证明重复节点消歧。occurrence 定义与重访节点示例、阈值、graph-only 补充及条件 WFA2 放方法或图注。README 保持用途、安装、使用、输出和限制的简洁结构。比对器分析见 [MAPPING_ASSESSMENT.zh-CN.md](MAPPING_ASSESSMENT.zh-CN.md)。

来源（2026-11-11 核对）：[vg deconstruct 官方文档](https://github.com/vgteam/vg/wiki/VCF-export-with-vg-deconstruct)、[pantree 仓库](https://github.com/oclb/pantree)、[pantree 论文](https://doi.org/10.1016/j.xgen.2026.101327)、[PGGB 仓库](https://github.com/pangenome/pggb)。
