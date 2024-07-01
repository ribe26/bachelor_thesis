function [output] = LSD(HRTF_1,HRTF_2,NFFT,Fs)
%UNTITLED2 この関数の概要をここに記述
%二つのHTRFの振幅スペクトル上での誤差(LSD)を計算する

resolution=Fs/NFFT;

k1=0;
k2=0;
count=0;

while count<3000
    count=count+resolution;
    k1=k1+1;
end

k2=k1;
while count<16000
    count=count+resolution;
    k2=k2+1;
end

k2=k2;

disp("k1:");
disp(k1);

disp("k2");
disp(k2)


calc=0;

for i=k1:k2
    calc=calc+(20*log10(abs(HRTF_1(i)/HRTF_2(i))))^2;
end
output=sqrt(calc/(k2-k1+2));
end