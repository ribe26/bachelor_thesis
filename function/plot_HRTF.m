function [f] = plot_HRTF(HRTF,fs)
%UNTITLED2 この関数の概要をここに記述
%  HRTFをプロットする
HRTF_abs=20*log10(abs(HRTF));
plot_y=HRTF_abs(1:length(HRTF_abs)/2);
plot_x=(1:length(HRTF_abs)/2);
f=figure;
plot(plot_x,plot_y);
ylim([-70 20]);
end