%HRTFのスペクトラルキューを仰角についての関数として表現する手法

% 初期化とパスの追加
clear;
addpath('./function');
addpath('./');

% サンプリング周波数と周波数ベクトルの設定
Fs = 48000;
Nout = Fs / 2;
f = (1:Fs/2);
f2 = (1:Fs);

% パラメータリストの設定
convoice_list = [0, 0; 10, 0; 20, 0; 30, 0; 40, 0; 50, 0; 60, 0; 70, 0; 80, 0; 90, 0; 80, 180; 70, 180; 60, 180; 50, 180; 40, 180; 30, 180; 20, 180; 10, 180; 0, 180];
person_list = ["subject1"];
time_list = [5, 10, 20, 30, 40, 50, 60];

% 各人物の処理
for pr = 1:length(person_list)
    person = person_list(pr);

    % パラメータファイルのパス設定
    P1_L_path = sprintf("./parameters/%s/P1_L.txt", person);
    N1_L_path = sprintf("./parameters/%s/N1_L.txt", person);
    N2_L_path = sprintf("./parameters/%s/N2_L.txt", person);
    N3_L_path = sprintf("./parameters/%s/N3_L.txt", person);

    P1_R_path = sprintf("./parameters/%s/P1_R.txt", person);
    N1_R_path = sprintf("./parameters/%s/N1_R.txt", person);
    N2_R_path = sprintf("./parameters/%s/N2_R.txt", person);
    N3_R_path = sprintf("./parameters/%s/N3_R.txt", person);

    % パラメータファイルの読み込み
    P1_L = readmatrix(P1_L_path);
    N1_L = readmatrix(N1_L_path);
    N2_L = readmatrix(N2_L_path);
    N3_L = readmatrix(N3_L_path);

    P1_R = readmatrix(P1_R_path);
    N1_R = readmatrix(N1_R_path);
    N2_R = readmatrix(N2_R_path);
    N3_R = readmatrix(N3_R_path);

    % 各時間リストの処理
    for tm = 1:length(time_list)
        time = time_list(tm);
        span = time / 180; % 角度を切り替えるスパン
        min_angle = 1;
        convoice = [0, 0];

        % 各仰角の処理
        for elevation = 0:min_angle:180
            % 白色雑音の生成
            white_noise = randn([round(Fs * span), 1]);
            white_noise = white_noise / max(abs(white_noise));
            vspec = fft(white_noise);
            pFs = round(length(white_noise) / 2);

            % パラメータ計算
            P_L = calc_parameters(elevation, P1_L, N1_L, N2_L, N3_L);
            P_R = calc_parameters(elevation, P1_R, N1_R, N2_R, N3_R);

            % PKフィルタの設計
            [~, sos_L] = PK(P_L(1, 1:3), Fs);
            [~, sos_R] = PK(P_R(1, 1:3), Fs);

            for i = 2:length(P_L(:, 1))
                [~, sos_PK_L] = PK(P_L(i, 1:3), Fs);
                sos_L = [sos_L; sos_PK_L];
            end

            for i = 2:length(P_R(:, 1))
                [~, sos_PK_R] = PK(P_R(i, 1:3), Fs);
                sos_R = [sos_R; sos_PK_R];
            end

            % 伝達関数の計算
            TF_L = ones(pFs, 1);
            TF_R = ones(pFs, 1);

            for i = 1:length(P_L(:, 1))
                TF_L = TF_L .* freqz(sos_L(i, 1:3), sos_L(i, 4:6), pFs);
            end

            for i = 1:length(P_R(:, 1))
                TF_R = TF_R .* freqz(sos_R(i, 1:3), sos_R(i, 4:6), pFs);
            end

            % HRTFの生成
            gen_HRTF_L = [TF_L; flip(TF_L)];
            gen_HRTF_R = [TF_R; flip(TF_R)];

            if mod(length(vspec), 2) == 1
                gen_HRTF_L = gen_HRTF_L(1:end-1);
                gen_HRTF_R = gen_HRTF_R(1:end-1);
            end

            % 逆FFTによるコンボリューション
            convoiceL = real(ifft(vspec .* gen_HRTF_L));
            convoiceR = real(ifft(vspec .* gen_HRTF_R));

            % 正規化
            convoiceL = convoiceL / max(abs(convoiceL));
            convoiceR = convoiceR / max(abs(convoiceR));

            % 合成音声の生成
            convoice = [convoice; convoiceL, convoiceR];
        end

        % 音声ファイルの保存
        save_path = sprintf("sounds/%s/%d/", person, time);
        mkdir(save_path);
        save_file = sprintf("sounds/%s/%d/0to180_parametric.wav", person, time);
        audiowrite(save_file, convoice, Fs);
    end
end

% 音声の再生（必要に応じてコメント解除）
%sound(convoice, Fs);
