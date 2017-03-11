%2016年8月3日10:33:18
%为ORB-SLAM提供数据集，tum数据集格式
%VFR=vision.VideoFileReader('H:\Interactivity\New Aerial Film\20160801相机相对位姿实验\DJI_0047.mp4');
%% Get frame from video
clc;
clear;
%% 读取视频
 video_file='G:\ComputerVision2016\3D Program\DataSet\MyDataset\Cellphone\VID_20170115_141356.mp4';
 video=VideoReader(video_file);
 frame_number=floor(video.Duration * video.FrameRate);
%% 提取
GETFRAME=true;
fid=fopen('rgb.txt','wt');
if(~GETFRAME)
fprintf(fid,'%s\n','# ZhangPeike 2016年8月3日11:23:33');
fprintf(fid,'%s\n','# DJI 0047');
fprintf(fid,'%s\n','timestamp name');
end
 for i=1:frame_number
     image_name=num2str(i);
     image_name=strcat(image_name,'.jpg');
     %读出图片
     I=read(video,i);
     if GETFRAME
     %写图片，分辨率0.5, 焦距和主点坐标1.5                           
     imwrite(imresize(I,0.5),image_name,'jpeg');
     I=[];
     else
        disp(video.CurrentTime);
        image_name=strcat('rgb/',image_name);
        fprintf(fid,'%.7f %s\n',video.CurrentTime,image_name);
     end
 end
fclose(fid);
