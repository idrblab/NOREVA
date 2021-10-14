#' @title Plotcircularplot Metabolomic Study with dataset without QCSs and ISs.
#' @description plot circular plot
#' @param data String, the input file name
#' @param outputfile Format string, together with cutoff value and type to generate formatted file name
#' @param cutoff Integer, which means to filter the results, the default is 100
#' @param outputtype Double-precision floating-point number, representing the characteristic value represented by the maximum length of the rectangle, the default is 40
#' @param maxValue Double-precision floating-point number, representing the characteristic value represented by the maximum length of the rectangle, the default is 40
#' @param colorSet Hexadecimal color string group, representing the four-layer color setting of the graph from the inside to the outside, the default is red (#EA4335), blue (#4285F4), yellow (#FBBC05), purple (#800080)
#' @param bgColor Hexadecimal color string, representing the background color of the graphic drawing, the default is white (#FFFFFF).
#' @param fontColor Hexadecimal color string, representing the font color, the default is black (#000000).
#' @param totalAngle Double-precision floating-point number, representing the total angle of rotation of the drawing, in degrees, the default value is 340
#' @import rJava
#' @importFrom utils read.csv
#' @usage norvisualization(data, outputfile,cutoff="100",outputtype,maxValue,
#' colorSet, totalAngle,bgColor, fontColor)
#' @export norvisualization
#' @examples
#' library(NOREVA)
#' \donttest{multi_qc_data <- PrepareInuputFiles(dataformat = 1,
#' rawdata = "Multiclass_with_QCS.csv")}
#' \donttest{normulticlassqcall(fileName = multi_qc_data,
#' assum_a="Y", assum_b="Y", assum_c="N")}
#' \donttest{norvisualization(data = "OUTPUT-NOREVA-Overall.Ranking.Data.csv",
#' cutoff = "100")}

#rownames(data) <- data[,1]
#data <- read.csv("OUTPUT-NOREVA-Overall.Ranking.Data.csv", header = TRUE)

norvisualization <- function(data, outputfile="NOREVA-Ranking-Top.%d.workflows.%s",cutoff="100",
                             outputtype="pdf",  maxValue="40",
                             colorSet = c("#FFA6A6", "#B5B5FF", "#C2DFB1", "#FFE699"), totalAngle = "340",
                             bgColor = "#FFFFFF", fontColor="#000000"
) {
  path_1 <- data
  #path_1 <- rawdataset
  pre_file2_1 <- readLines(path_1, n = 2)
  loc <- which.max(c(length(unlist(strsplit(pre_file2_1, ","))), length(unlist(strsplit(pre_file2_1, ";"))), length(unlist(strsplit(pre_file2_1, "\t")))))
  sep_seq <- c(",", ";", "\t")
  data <- read.csv(path_1,header=TRUE,sep=sep_seq[loc])
  data <- na.omit(data)
  ## 数据映射转化
  data<-data[,-c(2:6)]
  #data <- data[c("X","Criteria.Ca-Value", "Criteria.Cb-Value" ,"Criteria.Cc-Value", "Criteria.Cd-Value")]
  data["Criteria.Ca.Value"][data["Criteria.Ca.Value"]>=0.7]<-5
  data["Criteria.Ca.Value"][data["Criteria.Ca.Value"]<0.7&data["Criteria.Ca.Value"]>=0.3]<-10
  data["Criteria.Ca.Value"][data["Criteria.Ca.Value"]<0.3]<-40

  data["Criteria.Cb.Value"][data["Criteria.Cb.Value"]>=0.8]<-40
  data["Criteria.Cb.Value"][data["Criteria.Cb.Value"]<0.8&data["Criteria.Cb.Value"]>=0.5]<-10
  data["Criteria.Cb.Value"][data["Criteria.Cb.Value"]<0.5]<-5

  data["Criteria.Cc.Value"][data["Criteria.Cc.Value"]>=0.3]<-40
  data["Criteria.Cc.Value"][data["Criteria.Cc.Value"]<0.3&data["Criteria.Cc.Value"]>=0.15]<-10
  data["Criteria.Cc.Value"][data["Criteria.Cc.Value"]<0.15]<-5

  data["Criteria.Cd.Value"][data["Criteria.Cd.Value"]>=0.9]<-40
  data["Criteria.Cd.Value"][data["Criteria.Cd.Value"]<0.9&data["Criteria.Cd.Value"]>=0.7]<-10
  data["Criteria.Cd.Value"][data["Criteria.Cd.Value"]<0.7]<-5

  data <- as.matrix(data)


  .jinit(classpath = "inst/java/BioJar-3.8_jdk_1.8.jar", force.init = FALSE) ##启动JVM，没必要在函数里启动浪费效率

  BarDiagram <- J("biojar.function.graphics.CircularHistogram") ## 定义java的BarDiagram类，自定义类，绘图调用类
  Boolean <- J("java.lang.Boolean") ## 定义java的Boolean类，系统标准类，用于创建布尔对象参数
  String <- J("java.lang.String") ## 定义java的String类，系统标准类，用于创建字符串对象参数
  Integer <- J("java.lang.Integer") ## 定义java的Integer类，系统标准类，用于创建整型对象参数
  Double <- J("java.lang.Double") ## 定义java的Double类，系统标准类，用于创建双精度浮点数对象参数
  ## 将R矩阵转化为所需Java的ArrayList对象
  ArrayList <- J("java.util.ArrayList")
  criteriaList <- new(ArrayList)
  for (index in 1:length(data[,1])) {
    ObjectArray <- .jarray(data[index,])
    criteriaList$add(ObjectArray)
  }
  ## 实例化java绘图核心类BarDiagram
  bd <- new(BarDiagram)
  ## 装载输入数据
  bd$loadData(criteriaList,.jarray(c("Legend","Criteria.Ca.Value","Criteria.Cb.Value","Criteria.Cc.Value","Criteria.Cd.Value")),new(Integer, cutoff))
  ## 设置柱形图上限值
  bd$setMaxValue(new(Double, maxValue))
  ## 设置旋转总角度
  bd$setTotalAngle(new(Double, totalAngle))
  ## 设置图形填充颜色组
  bd$setBarColorSet(.jarray(colorSet))
  ## 设置背景色
  bd$setBackgroundColor(new(String, bgColor))
  ## 设置字体颜色
  bd$setFontColor(new(String, fontColor))
  ## 绘制图像
  bd$drawFigure(new(String, outputfile), new(Boolean, "false"), new (String, outputtype))

  message("The 'figure' has been successfully saved in the current path.")
  message("Please use 'getwd()' to find the current path!")

}

