# oriP：序列保留与倒位识别是不同层次

核验：2026-11-11；样本 E5_sLCL；参考窗口 B95_8_Raji:8404–8507（1-based，104 bp）。固定输出、逐样本 GT、原始 vg snarl ALT、PathAllele source allele 和原生 PGGB 路径均已核对。结果及文件哈希保存在 `oriP_case_audit.json`，可用 `python3 audit_orip_example.py --project-root /path/to/project` 重做本地序列回放。

## 实际观察

| 流程 | E5_sLCL 的携带记录 | 重建 104-bp ALT | 倒位表达 |
|---|---|---|---|
| PathAllele | 1 条，v2383，GT=2 | 与样本序列完全相同 | INV_COMPLEX；NATIVE=SEQUENCE；INVBASIS=LOCAL_RC |
| 当前 vg/wave | 2 条：82-bp DEL + 82-bp INS | 同一条序列 | 分解为两个长度事件 |
| 当前 Linear | 21 条 SNP / 小 indel | 同一条序列 | 没有在该窗口输出 INV |

回放只使用实际携带的 ALT，并核对 REF、坐标、非重叠及重建结果。没有把 missing 改为 0；该结果不构成全基因组 callability 结论。

vg 的原始 `>56518>65227` snarl 位于 6823，REF 长 2231 bp，E5_sLCL 的 ALT（GT=107）长 2383 bp，**已包含完整的 104-bp 样本序列**。因此不能说 vg 没有发现这段序列，更不能说图拓扑阻止它导出该 allele。其后分解后的表示为：8411 的 82-bp 缺失与 8495 的 82-bp 插入，两者共同回放正确。

Linear 的 21 条记录也回放正确。保存的比较 PAF 在该位点是正向 primary 比对，MAPQ 60，覆盖参考 [0,12302)，内部用多个 I/D 及 mismatch 表达差异。该 PAF 是保留的比较定位文件，不冒充当时 SVIM-asm 输入 BAM。

## 为什么没有作为 INV 输出

**vg/wave：** 导出阶段保留了序列。该项目的下一步使用 vcflib 1.0.15 `vcfwave -I 1000`。其反向判断比较整个 REF/ALT 的正向/反向 k-mer sketch，而不是先在一个较大的 allele 内搜索所有短局部倒位；未选反向时以 WFA 分解序列。这里的 84-bp 倒位嵌在超过 2 kb 的复合位点中，所以“局部倒位存在”不等于“整个 allele 更像反向参考”。这是源码与观察结果支持的机制解释，尚未通过修改该设置的消融实验证明唯一因果。

尤其不能简单说“84<1000，所以被 -I 1000 过滤”：原始外层 REF/ALT 都大于 1000 bp，该门槛作用于输入 allele，而不是其中倒位核心的长度。源码见 [vcfwave](https://github.com/vcflib/vcflib/blob/v1.0.15/src/vcfwave.cpp) 和 [vcf-wfa](https://github.com/vcflib/vcflib/blob/v1.0.15/src/vcf-wfa.cpp)。

**Linear：** 当前流程使用 minimap2/paftools 输出短变异，SVIM-asm 输出 SV。paftools 的序列差异表示可以保留 ALT，但不会自动将这一组短编辑重解释为本地复杂倒位。SVIM-asm 的倒位候选依赖相邻 alignment segments 的方向变化；如果 mapper 用一条正向长比对容纳局部差异，而没有足够的反向 split-alignment signature，就没有同等 INV 证据。这与保留 PAF 和最终调用一致；原始 SVIM BAM 未在本次检查中取得，因此不将其写成已经逐条证明的当时 BAM 因果。源码见 [SVIM-asm collect](https://github.com/eldariont/svim-asm/blob/master/src/svim_asm/SVIM_COLLECT.py) 和 [segment analysis](https://github.com/eldariont/svim-asm/blob/master/src/svim_asm/SVIM_inter.py)。

## PathAllele 是在哪一步识别的？

E5_sLCL 的来源 PC009732 是 mapped、正向、SEQUENCE 类型，内部未先被认定为 path-native INV。后续通用分类器 `src/classify.rs::Classifier::annotate` 搜索局部反向互补片段：使用参考唯一 canonical seed、检查 ALT 中种子唯一性、扩展反向匹配，并要求有效反向核心达到阈值。只有整段编辑恰好为同区间反向互补才是 INV；局部反向成分伴随边界/附加变化则为 INV_COMPLEX，证据标为 LOCAL_RC。

这套规则对所有 allele 执行。生产源码中没有 oriP、E5_sLCL、相关坐标或特定序列分支。已有随机序列测试覆盖精确/边界位移倒位、附加替换、远端反向插入、重复及回文反例；本次 `cargo test --locked --test classification` 的 6 个测试通过。该事实支持“通用规则而非特例条件”，不证明独立生物学准确率，也不证明该功能的开发动机与已知案例无关。

## 原生 PGGB 与展示图不能混为一谈

原生图中该窗口参考路径有 51 个相交节点访问，E5_sLCL 有 19 个；它们在各自路径内都只访问一次，方向均为 `+`。共同节点 64580+（CGG）分别位于参考 step 4000/offset 8403 和样本 step 4138/offset 8562（均 0-based）。这个窗口不是“重复节点消歧产生准确率提升”的证据。

用户提供的 oriP 图是由已验证变异重建的局部可视化图，其反向路径展示结构关系，不能当作原生 PGGB 的发现证据。实际 104-bp 序列可描述为：共同左侧 8 bp、TT→AC、84-bp 核心反向互补、共同右侧 10 bp。PathAllele 规范化记录边界与这个可视化分段不完全相同，不代表唯一历史断点。

算法主图最终只采用这个结构关系的**短合成模拟**：AGTC–CCGA → TCGG–GACT，加 TT→AC；不是原始序列的首尾缩写，也不标 84 bp、真实坐标或样本名。短模拟低于实际阈值，只解释信息流，不宣称它能在默认参数下触发 INV 分类。

可用于论文的结论：三种流程在该样本中保留了相同局部序列，但表示粒度和结构识别不同；PathAllele 的通用局部 RC 分类在完整 allele 上识别出复杂倒位。不能由一个已知实例推导整体灵敏度优势。
