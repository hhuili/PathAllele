# Panel a：紧凑算法示意

更新：2026-11-11。

沿用用户提供 oriP 图中的矩形碱基节点与曲线连接，只借鉴图形语言和“倒位伴随邻近替换”的结构关系。使用短合成序列：AGTC–CCGA → TCGG–GACT，另有 TT→AC。不得复用真实 84-bp 序列、样本名、坐标，或用省略号假装完成简化。

180 × 58.5 mm 单面板，三组对象从左向右排列：graph paths、mapping + occurrences、complete alleles。青色表示参考，紫色表示样本，黑色标出抽象锚点。输入为两条正向路径，后续完整 allele 的局部反向互补关系用交叉虚线和方向箭头表达，避免把原生图误画成已发现倒位。

图内只保留对象名称和必要步骤。occurrence 的严格定义、重复访问的歧义、缺失状态、证据阈值和分支条件放入图注/方法。不使用 Spelled contigs、REF FASTA = spell(REF path)、Same node ≠ same occurrence、Position · orientation · context 等用户已否定的标签。

此例不展示重复节点拷贝消歧的效果，也不宣称 vg/Linear 丢失了同一 ALT。短序列只解释结构关系，不作为默认参数下可运行的测试。
