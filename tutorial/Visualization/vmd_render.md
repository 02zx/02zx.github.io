# 如何在VMD中保存当前视角
```tcl
set view [molinfo top get {center_matrix rotate_matrix scale_matrix global_matrix}]
molinfo top set {center_matrix rotate_matrix scale_matrix global_matrix} $view
```

# 如何在VMD中添加字幕[^1]
```tcl
color Display Background white
display projection Orthographic
display depthcue off
display rendermode GLSL
display shadows on
display ambientocclusion on 
light 2 on
light 3 on
axes location Off


trace variable vmd_frame(1) w sdf
mol new
proc sdf {args} {global vmd_frame
graphics 0 delete all
graphics 0 color 16
graphics 0 text {0 0 0} [format "%4dps" [expr $vmd_frame(1)*1]] size 6 thickness 4}
```
color 16指黑色, 打开vmd后先输入以上代码然后再导入轨迹, 轨迹所在的ID必须是1. 

# ffmpeg制作视频：

ffmpeg -r 24 -start_number 0 -i pic.%05d.jpg -vframes 401 -crf 22 -c:v libx264 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" -pix_fmt yuv420p  video.mp4

若图片的长宽不是2的倍数，应使用 -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" 进行调整

-start_number 0 和 -vframes 401 控制了图片的范围，从0到400号用于视频合成

# ffmpeg中增加字幕
ffmpeg -r 24 -i noH.%05d.ppm -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2, drawtext=text='Time %{eif\:floor((n)/10)\:d}.%{eif\:mod(n\,10)\:d} ps': x=w-tw-10: y=h-th-10: fontsize=24: fontcolor=black" -c:v libx264 -pix_fmt yuv420p output.mp4

# ffmpeg将两个视频横向拼接
ffmpeg -i C2.mp4 -i output.mp4 -filter_complex "[0:v][1:v]hstack=inputs=2[v]" -map "[v]" -c:v libx264 combined.mp4


# 几何变换
```tcl
set sel [atomselect top all]
$sel move [transaxis z 90] #绕z转90度
set vec {0 -1 0}
$sel move [transvecinv $vec] #使体系中的一个向量平行于x轴
#注意VMD中c对应的是z轴方向的量
```

[^1]:http://sobereva.com/13
