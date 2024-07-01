function [HRTF,HRTF_abs] = earlyHRTF(HRIR,N,fs)
%UNTITLED この関数の概要をここに記述
%   earlyHRTFを求める
[early]=BlackHarris(HRIR,N);
HRTF=fft(early,fs);
HRTF_abs=abs(HRTF);
end