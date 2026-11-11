# PathAllele：路径 occurrence、完整等位序列与保守投影

更新：2026-11-11。源码基线 `057d45a9c51aba959ecd88eece279a047a02681d`（v0.2.0-dev.20261110）；本地与部署端的 graph/source/context/anchors 四个模块校验和一致。本文描述当前实现，未来方案见[比对器评估](MAPPING_ASSESSMENT.zh-CN.md)。

## 定位与适用范围

PathAllele 联合图路径 occurrence 与线性参考定位，提取来源明确的完整等位序列，并保守地投影为变异及基因型。当前开发重点是 PGGB 构建的小型单倍体组装泛基因组，主要评估对象为 EBV；其他病毒等小基因组是预期应用场景，不代表已逐一验证。大型真核泛基因组尚无规模和性能结论。小基因组不是硬性长度上限。

支持的输入是具有显式节点序列、P paths 和 0M links 的 GFA1，以及与完整参考路径一致的单条 A/C/G/T FASTA。一个单倍体样本可有多个 contig。当前不支持 GFA W walks、segment overlap、二倍体/多倍体基因型及环形参考坐标。**图中存在循环或路径重访节点，与采用环形参考坐标是两个概念。**

## 1. 节点与 occurrence 分开处理

PGGB 图可以包含循环，重复序列等结构可使一条路径多次经过同一节点。不能因此声称所有 PGGB 图或每条路径都一定重复经过节点，也不能把每个图循环自动解释为真实重复拷贝。

一个 occurrence 是：**某个有方向的节点在某条具体路径上的一次出现**。可用 `(path, step, oriented_node, sequence_offset)` 描述，其中 step 和 offset 区别同一路径上的重复访问。参考 occurrence 另外对应一个具体参考序列位置。

例：图中只有一个 X 节点，参考路径为 `L-X-U-X-R`，其 X 出现在 step 1 和 step 3。样本路径 `L-X-R` 的 X，不能仅凭节点 ID 与参考 X 配对，也不能用“样本第一次 X = 参考第一次 X”。上下文不足时，对应关系仍不确定。即使 L/R 已确定、完整 X-U-X/X allele 可以保存，也不必推断哪一个祖先 X 拷贝被删除。

程序沿已有 P 路径读取节点并按方向拼接序列，保留路径步数、节点方向与碱基偏移。它不搜索所有图中可能走法，不预先将所有节点分类为 anchor，也不先做 snarl 或 reference-tree 分解。L links 在读入时用于验证路径连接；下游主要使用已观察的路径与 occurrence 索引。实现：`src/graph.rs`。

## 2. 参考 FASTA 与参考路径

参考路径提供图节点、方向、重复访问及累积坐标；FASTA 提供参考序列和 contig 名称，供定位、参考索引及 VCF 检查/规范化使用。读入时强制校验：

```text
spell(reference_path) == reference_FASTA_sequence
```

因此外部 FASTA 在序列信息上不增加新的参考。理论上可从指定 REF path 自动生成，但当前 CLI 仍要求 `-r`，尚未实现自动生成输入参考的接口。显式 FASTA 也用于检查图参考与预期参考版本是否一致。参考路径的默认选择与 FASTA contig 名称相关，必要时显式使用 `--reference-path`；不要将当前实现描述成通过全图序列搜索自动选择参考路径。

## 3. 线性定位与图 occurrence 联合确定来源

`source::build` 将每个样本的每条 contig 路径拼为 `paths.fa`，再将完整 contig 序列比对到 `reference.fa`。不是把 graph 本身或仅候选锚点比对到 FASTA。

当前调用包括 `minimap2 -cx asm20 -k 19 -s 40 --end-bonus 100 --no-hash-name --eqx --secondary=yes -N 20`（部分值来自配置）。`context::read_alignments` 只保留达到 MAPQ 阈值（默认 20）的 `tp:A:P` 记录；**secondary 虽然输出，却不作为来源候选进入当前解析器**。CIGAR 中的精确匹配段建立 query→reference 位置索引，不直接把 I/D/X 当作最终变异。

`anchors::occurrence_seeds` 在精确匹配上下文内检查样本节点是否对应参考上同方向、同节点且坐标一致的 occurrence。`context_runs` 合并连续匹配、执行长度及序列唯一性筛选，再求保守共线链。接受的来源还要处理重叠 query ownership、竞争物理拷贝和结构关系。

未被已接受来源覆盖的区域，可通过唯一节点、方向、唯一序列及间距限制形成 graph-only 补充链。该分支不是对所有漏定位的全图穷举恢复，也不是完全独立于 mapping 的全图来源判定。当前候选空间受 minimap2、primary/MAPQ 筛选及 graph-only 条件共同限制。

## 4. 锚点是候选链共同必经的匹配，不是共享节点的别名

候选锚点为 `(r, q, n)`。共线链要求参考与 query 位置均有序且不重叠，得分为匹配长度之和。`anchors::obligatory` 用前向/反向动态规划和 BigUint 链计数，保留所有最高分候选链共同经过的锚点。

这个保证只针对**已构建并通过资格筛选的候选集合与当前评分规则**，不是所有生物学同源关系、所有基因组位置或所有最优碱基比对。重复节点的单次出现不会自动成为锚点；映射、相邻上下文和序列唯一性条件仍适用。局部内部图锚点及唯一序列锚点另用于细分完整 allele。

## 5. 先保留完整 allele，再处理 event 表示

`source::Collector::accept` 从相邻可信锚点之间提取完整 REF/ALT，并保存样本、原始路径、query/ref 区间、路径步区间和方向。倒位/复杂结构证据使用节点方向及序列关系，不由长度差单独决定。

`alleles::prepare/decompose` 做锚点分解、正向与反向互补方向处理、参考规范化及完整序列回放。当前代码确实使用 `bcftools norm` 左规范化；不要沿用其他流程“不独立左对齐”的约定。要求：

```text
apply(reference_window, edits) == observed_ALT
```

边界扩展、方向差异或编辑组合不稳定时，程序可保留更大的 compound allele，并对无法确定的局部投影保留不确定 mask。回放只验证输出表示与观察序列一致；错误组装或错误来源也可能满足回放。

`partial()` 检测消耗参考碱基的编辑区间是否缺乏完整的 known.unique 覆盖，必要时调用内置 vcflib/WFA2 做局部投影。它是条件分支，不是每个完整 allele 必经的发现步骤。WFA2/BiWFA 是既有比对实现，不是 PathAllele 新发明的比对算法。

事件分类还有独立的通用序列证据层：`classify::Classifier::annotate` 检查局部反向互补片段。因此 `NATIVE=SEQUENCE; INVBASIS=LOCAL_RC` 表示来源阶段没有给出原生图倒位类别，后续局部 RC 证据支持倒位成分；不是“最初已检测到 path-native INV”。精确同区间反向互补可为 INV，伴随边界/附加变化可为 INV_COMPLEX。当前生产代码没有针对 oriP 名称、坐标或样本的条件分支。已知案例及其与 vg/Linear 的区别见 [oriP 核验](ORIP_CASE_STUDY.zh-CN.md)；此案例是回归控制，不是独立准确率证据。

## 6. 基因型及缺失状态

| 输出 | 语义 |
|---|---|
| `Graph.loci.vcf.gz` / `Graph.loci.tsv.gz` | 完整局部序列身份；GT/FT 是 complete-sequence 状态；PG/PF 是 event-projection 状态 |
| `Graph.vcf.gz` | exact event 的 GT/FT；某记录 GT=0 不必意味着整个局部窗口都是参考序列 |
| `Graph.sv.vcf.gz` / members | emitted SV alleles 的兼容性聚类及样本状态，不是祖先突变事件 |

无法确定不等于参考。`NO_DATA`、`MULTIPLE_OCCURRENCES`、`AMBIGUOUS_ALIGNMENT`、`OVERLAPPING_EDIT`、`PARTIAL_SEQUENCE`、`UNRESOLVED_STRUCTURE` 等状态区分不同证据缺口。完整序列身份与事件投影是不同层次，某事件投影不确定不必抹去已确定的完整序列身份。

## 7. 方法贡献的边界

完整 allele、序列回放和保守 missingness 都可由其他 assembly-based caller 实现。PathAllele 当前可辨认的特点，是把这些处理与图路径 occurrence、来源定位及保守锚点链相结合，尤其针对 PGGB 路径重访节点时的物理位置歧义。其他工具也可处理参考循环，不能将“支持有环图”单独当作独占创新。

源码证明行为，不证明创新优先权或准确率优势。需要独立真值与消融，分别测量来源对应准确性、完整 allele 正确性、callability、missingness 和资源消耗；已知 EBV 位点仅作为回归控制。工具比较见[方法定位](METHOD_POSITIONING.zh-CN.md)。
