% スペクトラルキューに基づく手法で推定したHRTFと実測HRTFをプロットする
clear;
addpath('./function');
addpath('./');

% サンプリング周波数とFFTパラメータの設定
Fs = 48000;
NFFT = 512;
Nout = NFFT / 2 + 1;
f = ((1:Nout-1)' ./ NFFT) .* Fs;
f2 = (1:Fs/2);

% パラメータの設定
person = "subject1";
convoice_list = [0, 0; 10, 0; 20, 0; 30, 0; 40, 0; 50, 0; 60, 0; 70, 0; 80, 0; 90, 0; 80, 180; 70, 180; 60, 180; 50, 180; 40, 180; 30, 180; 20, 180; 10, 180; 0, 180];
person_list = ["subject1"];

% HRTFデータの初期化
HRTF_L = [];
HRTF_R = [];
TF_L = [];
TF_R = [];

% HRTFデータの読み込み
for i = 1:length(convoice_list(:,1))
    elevation = convoice_list(i, 1);
    azimuth = convoice_list(i, 2);
    
    IRfname = sprintf("./hrtf/%s/elev%d/R%de%03da_new.dat", person, elevation, elevation, azimuth);
    ILfname = sprintf("./hrtf/%s/elev%d/L%de%03da_new.dat", person, elevation, elevation, azimuth);
    
    Rdat = fopen(IRfname, 'r', 'b');
    IRdata = fread(Rdat, 'float');
    IR_spec = fft(IRdata, NFFT);
    fclose(Rdat);
    
    Ldat = fopen(ILfname, 'r', 'b');
    ILdata = fread(Ldat, 'float');
    IL_spec = fft(ILdata, NFFT);
    fclose(Ldat);

    HRTF_L = [HRTF_L, IL_spec];
    HRTF_R = [HRTF_R, IR_spec];
end

% 低角度と高角度のHRTFから任意の角度のHRTFを推定
lowAngle = 100;
highAngle = 110;
elevation = 105;

low_index=fix(lowANgle/10)+1;
high_index=fix(highAngle/10)+1;

new_HRTF = calc_angle_HRTF(HRTF_L, elevation, lowAngle, highAngle);

% 実測HRTFと推定HRTFのプロット
f3 = figure;
hold on;
semilogx(f, 20*log10(abs(HRTF_L(1:length(IR_spec)/2, 11))), "-.r", "LineWidth", 0.8, 'DisplayName', 'HRTF_{90}');
semilogx(f, 20*log10(abs(new_HRTF(1:length(IR_spec)/2))), "-b", "LineWidth", 1, 'DisplayName', 'HRTF_{synthesised}');
semilogx(f, 20*log10(abs(HRTF_L(1:length(IR_spec)/2, 12))), "-.m", "LineWidth", 0.8, 'DisplayName', 'HRTF_{100}');
legend;
axis([200 Fs/2 -60 20]);
xlabel('Frequency [Hz]');
ylabel('Relative level [dB]');

% プロットの保存
path = sprintf("./graph/%s/linear/", person);
mkdir(path);
savename = sprintf("%s%dHRTF.png", path, elevation);
saveas(f3, savename);
