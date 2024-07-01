%二つのHRTFを線形2点補完することで、任意仰角のHRTFを推定する手法
% 初期化とパスの追加
clear;
addpath('./function');
addpath('./');

% サンプリング周波数とFFTパラメータの設定
Fs = 48000;
NFFT = 512;
Nout = NFFT / 2 + 1;
f = ((1:Nout-1)' ./ NFFT) .* Fs;
f2 = (1:Fs/2);

% パラメータリストの設定
convoice_list = [0, 0; 10, 0; 20, 0; 30, 0; 40, 0; 50, 0; 60, 0; 70, 0; 80, 0; 90, 0; 80, 180; 70, 180; 60, 180; 50, 180; 40, 180; 30, 180; 20, 180; 10, 180; 0, 180];
person_list = ["subject1"];
%time_list = [5, 10, 20, 30, 40, 50, 60];
time_list=[30];
% 各人物の処理
for pr = 1:length(person_list)
    person = person_list(pr);

    % 各時間リストの処理
    for tm = 1:length(time_list)
        time = time_list(tm);
        span = time / 180; % 角度を切り替えるスパン
        min_angle = 1;

        HRTF_L = [];
        HRTF_R = [];
        TF_L = [];
        TF_R = [];
        
        pfs = round(Fs * span);

        % HRTFデータの読み込み
        for i = 1:length(convoice_list(:, 1))
            elevation = convoice_list(i, 1);
            azimuth = convoice_list(i, 2);
            
            IRfname = sprintf("./hrtf/%s/elev%d/R%de%03da_new.dat", person, elevation, elevation, azimuth);
            ILfname = sprintf("./hrtf/%s/elev%d/L%de%03da_new.dat", person, elevation, elevation, azimuth);
            
            Rdat = fopen(IRfname, 'r', 'b');
            IRdata = fread(Rdat, 'float');
            IR_spec = fft(IRdata, pfs);
            fclose(Rdat);
            
            Ldat = fopen(ILfname, 'r', 'b');
            ILdata = fread(Ldat, 'float');
            IL_spec = fft(ILdata, pfs);
            fclose(Ldat);

            HRTF_L = [HRTF_L, IL_spec];
            HRTF_R = [HRTF_R, IR_spec];
        end

        % コンボリューション音声の生成
        convoice = [0, 0];

        for elevation = 0:min_angle:179
            % 白色雑音の生成と正規化
            white_noise = randn([round(Fs * span), 1]);
            white_noise = white_noise / max(abs(white_noise));
            vspec = fft(white_noise);
            pfs = length(vspec);

            % 線形補間により仰角=elevationの角度のHRTFを推定
            lowAngle=fix(elevation/10)*10;
            highAngle=lowAngle+10;
            ely_L = calc_angle_HRTF(HRTF_L, elevation,lowAngle,highAngle);
            ely_R = calc_angle_HRTF(HRTF_R, elevation,lowAngle,highAngle);

            % 白色雑音にHRTFを畳み込み
            convoiceL = real(ifft(vspec .* ely_L));
            convoiceR = real(ifft(vspec .* ely_R));

            % 正規化
            convoiceL = convoiceL / max(abs(convoiceL));
            convoiceR = convoiceR / max(abs(convoiceR));

            % 合成音声の更新
            convoice = [convoice; convoiceL, convoiceR];
        end

        % 音声ファイルの保存
        save_path = sprintf("sounds/%s/%d/", person, time);
        mkdir(save_path);
        save_file = sprintf("sounds/%s/%d/0to180_linear.wav", person, time);
        audiowrite(save_file, convoice, Fs);
    end
end
