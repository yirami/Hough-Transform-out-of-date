# Hough-Transform
The source code for detecting and tracking the lane lines.
车道线检测与追踪的源代码，部分资源亦用于论文。

---
## 说明
+ MATLAB程序分两部分，`直线检测` `车道线检测与跟踪`
+ `./src/`文件夹下为用于编译生成`mex`文件的`c`源代码  
+ `需要`将程序`根目录`加入MATLAB路径  
+ `./LineDetection/`文件夹下为改进的直线检测算法DCHT的MATLAB脚本
+ `./LaneD&T/`文件夹下为基于DCHT的车道线检测与追踪的MATLAB脚本
+ `./`根目录下是编译好了的`mex`文件
---
### 第一部分：`直线检测`
改进传统Hough变换->DCHT
#### 展示
1. 投票量（右侧为改进型）  
![道路图片投票量](https://i.loli.net/2017/09/12/59b7b4bc0aa28.jpg)
+ 效果
  + 道路  
  ![道路检测效果](https://i.loli.net/2017/09/12/59b7b4c09283b.jpg)
  + 房子  
  ![房子检测效果](https://i.loli.net/2017/09/12/59b7b4bc0c210.jpg)
  + 室内  
  ![室内检测效果](https://i.loli.net/2017/09/12/59b7b4bc11297.jpg)
+ 耗时  
![检测耗时](https://i.loli.net/2017/09/12/59b7b4bbee8d2.jpg)
---
### 第二部分：`车道线检测与跟踪`
基于DCHT算法，实时处理在高速公路上手机采集的视频。  
#### 展示
+ GIF片段  
![片段2](https://i.loli.net/2017/09/12/59b7a39881ff9.gif)  
![片段1](https://i.loli.net/2017/09/12/59b7a397a80ea.gif)
+ 源视频  
[百度云](https://pan.baidu.com/s/1gfzmme3)，*提取密码*：`jayc`
