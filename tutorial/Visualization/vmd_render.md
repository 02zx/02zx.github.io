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


[^1]:http://sobereva.com/13
