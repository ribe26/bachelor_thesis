%スペクトラルキューに基づく手法で推定したHRTFを実測値のをプロットするコード
clear;
addpath('./function');
addpath('./');

% FFTのパラメータの設定
Fs = 48000;
NFFT = 512;
Nout = NFFT / 2 + 1;
f = ((1:Nout-1)' ./ NFFT) .* Fs;
f2 = (1:Fs/2);

% HRTFを読み込む角度のリストの設定
person = "subject1";
convoice_list = [0, 0; 10, 0; 20, 0; 30, 0; 40, 0; 50, 0; 60, 0; 70, 0; 80, 0; 90, 0; 80, 180; 70, 180; 60, 180; 50, 180; 40, 180; 30, 180; 20, 180; 10, 180; 0, 180];
convoice_angle = [0:10:180];

HRTF_L = [];
HRTF_R = [];
TF_L = [];
TF_R = [];

% 提案手法でパラメータと仰角の関係を近似した関数のパラメータの読み込み
P1_L_path = sprintf("./parameters/%s/P1_L.txt", person);
N1_L_path = sprintf("./parameters/%s/N1_L.txt", person);
N2_L_path = sprintf("./parameters/%s/N2_L.txt", person);
N3_L_path = sprintf("./parameters/%s/N3_L.txt", person);
P1_R_path = sprintf("./parameters/%s/P1_R.txt", person);
N1_R_path = sprintf("./parameters/%s/N1_R.txt", person);
N2_R_path = sprintf("./parameters/%s/N2_R.txt", person);
N3_R_path = sprintf("./parameters/%s/N3_R.txt", person);

P1_L = readmatrix(P1_L_path);
N1_L = readmatrix(N1_L_path);
N2_L = readmatrix(N2_L_path);
N3_L = readmatrix(N3_L_path);
P1_R = readmatrix(P1_R_path);
N1_R = readmatrix(N1_R_path);
N2_R = readmatrix(N2_R_path);
N3_R = readmatrix(N3_R_path);

% 各仰角で近似したHRTFを実測値HRTFを比較してプロットする
for elevation = 1:length(convoice_list(:,1))
    % 右耳HRIRを読み込みしてHRTFに変換する
    IR_R_fname = sprintf("./hrtf/%s/elev%d/R%de%03da_new.dat", person, convoice_list(elevation,1), convoice_list(elevation,1), convoice_list(elevation,2));
    Rdat = fopen(IR_R_fname, 'r', 'b');
    IR_R_data = fread(Rdat, 'float');
    fclose(Rdat);
    IR_R_spec = fft(IR_R_data, Fs);

    % 左耳HRIRを読み込みしてHRTFに変換する
    IR_L_fname = sprintf("./hrtf/%s/elev%d/L%de%03da_new.dat", person, convoice_list(elevation,1), convoice_list(elevation,1), convoice_list(elevation,2));
    Ldat = fopen(IR_L_fname, 'r', 'b');
    IR_L_data = fread(Ldat, 'float');
    fclose(Ldat);
    IR_L_spec = fft(IR_L_data, Fs);

    % パラメータからHRTFを近似した双二次フィルタの係数を計算する
    P_L = calc_parameters(convoice_angle(elevation), P1_L, N1_L, N2_L, N3_L);
    P_R = calc_parameters(convoice_angle(elevation), P1_R, N1_R, N2_R, N3_R);

    [~, sos_L] = PK(P_L(1,1:3), Fs);
    [~, sos_R] = PK(P_R(1,1:3), Fs);

    for i = 2:length(P_L(:,1))
        [~, sos_PK_L] = PK(P_L(i,1:3), Fs);
        sos_L = [sos_L; sos_PK_L];
    end

    for i = 2:length(P_R(:,1))
        [~, sos_PK_R] = PK(P_R(i,1:3), Fs);
        sos_R = [sos_R; sos_PK_R];
    end

    %各ピークやノッチを近似した双二次フィルタを重ね合わせて近似HRTFを生成
    TF_L = ones(Fs/2, 1);
    TF_R = ones(Fs/2, 1);

    for i = 1:length(P_L(:,1))
        TF_L = TF_L .* freqz(sos_L(i,1:3), sos_L(i,4:6), Fs/2);
    end

    for i = 1:length(P_R(:,1))
        TF_R = TF_R .* freqz(sos_R(i,1:3), sos_R(i,4:6), Fs/2);
    end

    % 右耳のHRTFをプロットする
    figure;
    hold on;
    semilogx(f2, 20*log10(abs(IR_R_spec(1:length(IR_R_spec)/2))), 'DisplayName', 'HRTF');
    semilogx(f2, 20*log10(abs(TF_R(1:Fs/2))), 'DisplayName', 'parametric HRTF');
    title(sprintf("HRTF Right Ear Elevation %d°", convoice_angle(elevation)));
    legend;
    axis([90 Fs/2 -60 20]);
    xlabel('Frequency [Hz]');
    ylabel('Relative level [dB]');
    path = sprintf("./graph/%s/parametric/R/", person);
    mkdir(path);
    savename = sprintf("%s%d.png", path, elevation);
    saveas(gcf, savename);
    close(gcf);

    % 左耳のHRTFをプロットする
    figure;
    hold on;
    semilogx(f2, 20*log10(abs(IR_L_spec(1:length(IR_L_spec)/2))), 'DisplayName', 'HRTF');
    semilogx(f2, 20*log10(abs(TF_L(1:Fs/2))), 'DisplayName', 'parametric HRTF');
    title(sprintf("HRTF Left Ear Elevation %d°", convoice_angle(elevation)));
    legend;
    axis([90 Fs/2 -60 20]);
    xlabel('Frequency [Hz]');
    ylabel('Relative level [dB]');
    path = sprintf("./graph/%s/parametric/L/", person);
    mkdir(path);
    savename = sprintf("%s%d.png", path, elevation);
    saveas(gcf, savename);
    close(gcf);
end
