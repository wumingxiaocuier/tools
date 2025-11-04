# 加载必要的包
library(UpSetR)

# 您的数据 - 按您要求的顺序排列
expression_data <- c(
  "pFind-new&MSFragger&MaxQuant&Comet" = 115099,  # 所有集合均为True
  "pFind-new&MSFragger&MaxQuant" = 5370,  # pFind-new、MSFragger、MaxQuant为True
  "pFind-new&MSFragger&Comet" = 10888,  # pFind-new、MSFragger、Comet为True
  "pFind-new&MaxQuant&Comet" = 2773,  # pFind-new、MaxQuant、Comet为True
  "MSFragger&MaxQuant&Comet" = 922,  # MSFragger、MaxQuant、Comet为True
  "pFind-new&MSFragger" = 6578,  # pFind-new、MSFragger为True
  "pFind-new&MaxQuant" = 4653,  # pFind-new、MaxQuant为True
  "pFind-new&Comet" = 1357,  # pFind-new、Comet为True
  "MSFragger&MaxQuant" = 1508,  # MSFragger、MaxQuant为True
  "MSFragger&Comet" = 709,  # MSFragger、Comet为True
  "MaxQuant&Comet" = 531,  # MaxQuant、Comet为True
  "pFind-new" = 28783,  # 仅pFind-new为True
  "MSFragger" = 4105,  # 仅MSFragger为True
  "MaxQuant" = 2421,  # 仅MaxQuant为True
  "Comet" = 1596  # 仅Comet为True
)

# 创建颜色向量
# 第1个（4个集合）：d颜色
# 第2-5个（3个集合）：c颜色
# 第6-11个（2个集合）：b颜色
# 第12-15个（1个集合）：a颜色
bar_colors <- c(
  rep("purple", 1),      # 第1个 - 4个集合
  rep("green", 4),       # 第2-5个 - 3个集合
  rep("blue", 6),        # 第6-11个 - 2个集合
  rep("red", 4)          # 第12-15个 - 1个集合
)

# 绘制图形
upset(fromExpression(expression_data),
      nsets = 4,
      nintersects = 15,
      sets = c("pFind-new", "MSFragger", "MaxQuant", "Comet"),  # 左侧集合顺序
      order.by = "degree",
      mainbar.y.label = "Intersection Size",
      sets.x.label = "Set Size",
      text.scale = c(1.3, 1.3, 1, 1, 1.5, 1.5),
      # main.bar.col = bar_colors,  # 设置柱子颜色
      # matrix.color = bar_colors,
      queries = list(
        # 高亮1：同时包含pFind-new和MSFragger的交集（红色）
        list(query = intersects, params = list("pFind-new", "MSFragger", "MaxQuant", "Comet"), 
            color = "#22AC4A", active = TRUE),
        # 高亮2：仅包含Comet的交集（绿色）
        list(query = intersects, params = list("pFind-new", "MSFragger", "MaxQuant"), 
            color = "#87C441", active = TRUE),
        list(query = intersects, params = list("pFind-new", "MSFragger", "Comet"), 
            color = "#87C441", active = TRUE),
        list(query = intersects, params = list("pFind-new", "MaxQuant", "Comet"), 
            color = "#87C441", active = TRUE),
        list(query = intersects, params = list("MSFragger", "MaxQuant", "Comet"), 
            color = "#87C441", active = TRUE),
        # 高亮3：同时包含MaxQuant、Comet的交集（蓝色）
        list(query = intersects, params = list("pFind-new", "MSFragger"), 
            color = "#FBC611", active = TRUE),
        list(query = intersects, params = list("pFind-new", "MaxQuant"), 
            color = "#FBC611", active = TRUE),
        list(query = intersects, params = list("pFind-new", "Comet"), 
            color = "#FBC611", active = TRUE),
        list(query = intersects, params = list("MSFragger", "MaxQuant"), 
            color = "#FBC611", active = TRUE),
        list(query = intersects, params = list("MSFragger", "Comet"), 
            color = "#FBC611", active = TRUE),
        list(query = intersects, params = list("MaxQuant", "Comet"), 
            color = "#FBC611", active = TRUE),
        # 可继续添加：想高亮哪个交集，就复制一个list改参数
        list(query = intersects, params = list("pFind-new"), 
            color = "#F26723", active = TRUE),
        list(query = intersects, params = list("MSFragger"), 
            color = "#F26723", active = TRUE),
        list(query = intersects, params = list("MaxQuant"), 
            color = "#F26723", active = TRUE),
        list(query = intersects, params = list("Comet"), 
            color = "#F26723", active = TRUE)
      ),
      point.size = 5,             # 点的大小
      line.size = 1.2               # 线的粗细
)