# 逐页中文讲稿

## 1. iNKT cytotoxicity and Gene Ontology
本轮新增的是细胞毒性与GO计算。输入来自9月5日最新分析对象；9月15日仅是展示更新。

## 2. Analysis design and reference mapping
采用论文中的方法和图形组织方式。没有运行SCORPION、随机森林或声称重现PAGER专有网络；GO采用官方MGI小鼠注释。

## 3. Cytotoxic gene coverage and scoring
分数是转录表达而不是杀伤实验读数。不同基因集的绝对分数不可视为相同量纲的功能强弱。

## 4. Bone marrow: cytotoxic expression and condition effects
按组织和原分群展示细胞毒性；红色为T2高，蓝色为Ctrl高。星号是细胞层面探索性检验，不能代表动物重复验证。

## 5. Spleen: cytotoxic expression and condition effects
按组织和原分群展示细胞毒性；红色为T2高，蓝色为Ctrl高。星号是细胞层面探索性检验，不能代表动物重复验证。

## 6. Thymus: cytotoxic expression and condition effects
按组织和原分群展示细胞毒性；红色为T2高，蓝色为Ctrl高。星号是细胞层面探索性检验，不能代表动物重复验证。

## 7. Bone marrow: GO biological processes: cluster identity
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 8. Spleen: GO biological processes: cluster identity
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 9. Thymus: GO biological processes: cluster identity
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 10. Bone marrow: GO biological processes: condition response
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 11. Spleen: GO biological processes: condition response
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 12. Thymus: GO biological processes: condition response
GO气泡图只展示显著且至少两个命中基因的代表条目。T2/ Ctrl列分别表示该方向上调基因的富集，不把富集直接解释成通路激活。完整BP、MF、CC结果见数据表。

## 13. Bone marrow: local cytotoxicity along pseudotime
重新计算了组织内diffusion map，避免复用全局轨迹。用前体marker选择根群，三个根细胞做敏感性。平滑曲线不是置信区间，条件曲线差别不等于分化速率差别。

## 14. Spleen: local cytotoxicity along pseudotime
重新计算了组织内diffusion map，避免复用全局轨迹。用前体marker选择根群，三个根细胞做敏感性。平滑曲线不是置信区间，条件曲线差别不等于分化速率差别。

## 15. Thymus: local cytotoxicity along pseudotime
重新计算了组织内diffusion map，避免复用全局轨迹。用前体marker选择根群，三个根细胞做敏感性。平滑曲线不是置信区间，条件曲线差别不等于分化速率差别。

## 16. Core signature changes do not imply a full killing pathway
骨髓C4和脾脏C3的核心3基因模块及原4基因效应模块有一致局部变化；加入Ctla2a后的论文E1代表模块不显著。杀伤相关GO也未显著，不能推断整个细胞毒性通路增强或减弱。

## 17. Cellular components of cytotoxic module genes
这里只画真实GO细胞组分注释矩阵，不把共同定位误画为PPI相互作用。

## 18. Robustness: QC matching and alternative clustering
绿色范围表示五次QC匹配结果范围，不是置信区间。完整稳定分群分析保存在表里，不把S群直接等同原来的C群。

## 19. Evidence summary and next interpretation
总结按表达差异绝对值排序。优先看基因组成、QC匹配、背景校正评分是否同向，再讨论GO机制；不要把转录结果写成已验证杀伤增强或减弱。
