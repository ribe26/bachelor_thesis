%スペクトラルキューのパラメータ(ゲイン、バンド幅、中心周波数)を仰角についての関数として近似した際の
%係数を保存しておく手法 保存先は/parameter内
clear;
addpath('./function');

% 定数の設定
Fs = 48000;
NFFT = 512;
Nout = NFFT / 2 + 1;
f = ((1:Nout-1)' ./ NFFT) * Fs;
f2 = (1:Fs);
person_list = ["subject1"];
dirc_list = ["R", "L"];
convoice_list = [0,0; 10,0; 20,0; 30,0; 40,0; 50,0; 60,0; 70,0; 80,0; 90,0; 80,180; 70,180; 60,180; 50,180; 40,180; 30,180; 20,180; 10,180; 0,180];
angle_total = length(convoice_list(:,1));

% メインループ
for pr = 1:length(person_list)
    person = person_list(pr);
    for dr = 1:length(dirc_list)
        dirc = dirc_list(dr);
        
        P1 = [];
        N1 = [];
        N2 = [];
        N3 = [];

        % 各エレベーション角度で処理
        for elevation = 1:length(convoice_list(:,1))
            [P1, N1, N2, N3] = process_elevation(person, dirc, convoice_list, elevation, Fs, NFFT, P1, N1, N2, N3);
        end

        % N2, N3が見つからなかった区間のパラメータを線形に補間する
        originalN2=N2;
        originalN3=N3;
        N2 = interpolate_missing_params(convoice_list, N2, angle_total);
        N3 = interpolate_missing_params(convoice_list, N3, angle_total);
        
        N = 10; % 次数
        save_parameters(P1,N,person,dirc,'P1');
        save_parameters(N1,N,person,dirc,'N1');
        save_parameters(N2,N,person,dirc,'N2');
        save_parameters(N3,N,person,dirc,'N3');
        
        
        x = linspace(0, 180, 181);
        plot_param(P1, N, x, person, dirc, 'P1', 'W');
        plot_param(P1, N, x, person, dirc, 'P1', 'B');
        plot_param(P1, N, x, person, dirc, 'P1', 'G');
        
        plot_param(N1, N, x, person, dirc, 'N1', 'W');
        plot_param(N1, N, x, person, dirc, 'N1', 'B');
        plot_param(N1, N, x, person, dirc, 'N1', 'G');
        
        plot_param(N2, N, x, person, dirc, 'N2', 'W', originalN2);
        plot_param(N2, N, x, person, dirc, 'N2', 'B', originalN2);
        plot_param(N2, N, x, person, dirc, 'N2', 'G', originalN2);
        
        plot_param(N3, N, x, person, dirc, 'N3', 'W', originalN3);
        plot_param(N3, N, x, person, dirc, 'N3', 'B', originalN3);
        plot_param(N3, N, x, person, dirc, 'N3', 'G', originalN3);
    end
end

% 各エレベーション角度の処理を行う関数
function [P1, N1, N2, N3] = process_elevation(person, dirc, convoice_list, elevation, Fs, NFFT, P1, N1, N2, N3)
    IRfname=sprintf("./hrtf/%s/elev%d/%s%de%03da_new.dat",person,convoice_list(elevation,1),dirc,convoice_list(elevation,1),convoice_list(elevation,2));
    Rdat = fopen(IRfname, 'r', 'b');
    IRdata = fread(Rdat, 'float');
    IR_spec = fft(IRdata, NFFT);

    [ely, ely_abs] = earlyHRTF(IRdata, 44, NFFT);
    BF = 20*log10(ely_abs(1:NFFT/2));

    % ピークとノッチを近似した双二次フィルタのパラメータの読み込み
    datafile = sprintf("./nums/%s/nums_%d_%d_%s.txt", person, convoice_list(elevation,1), convoice_list(elevation,2), dirc);
    P = readmatrix(datafile);
    P(:,1) = P(:,1) * Fs / (2 * pi);

    % 一番中心周波数が小さいピークとノッチをP1, N1とする
    P1 = [P1; (elevation-1)*10, convoice_list(elevation,1), convoice_list(elevation,2), P(1,:)];
    N1 = [N1; (elevation-1)*10, convoice_list(elevation,1), convoice_list(elevation,2), P(2,:)];

    % N2, N3が存在するかの探索
    [N2, N3] = find_additional_notches(P, N2, N3, elevation, convoice_list);
end

% N2, N3が存在するかの探索を行う関数
function [N2, N3] = find_additional_notches(P, N2, N3, elevation, convoice_list)
    foundN2 = 0;
    foundN3 = 0;
    
    if length(P(:,1)) > 2
        for i = 3:length(P(:,1))
            if P(i,1) < 17000 && foundN2 == 0
                N2 = [N2; (elevation-1)*10, convoice_list(elevation,1), convoice_list(elevation,2), P(i,:)];
                foundN2 = 1;
            elseif P(i,1) >= 17000 && foundN3 == 0
                N3 = [N3; (elevation-1)*10, convoice_list(elevation,1), convoice_list(elevation,2), P(i,:)];
                foundN3 = 1;
            end
        end
    end

    if foundN2 == 0
        N2 = [N2; -1, convoice_list(elevation,1), convoice_list(elevation,2), 0, 0, 0];
    end
    if foundN3 == 0
        N3 = [N3; -1, convoice_list(elevation,1), convoice_list(elevation,2), 0, 0, 0];
    end
end

% N2, N3が見つからなかった区間のパラメータを線形に補間する関数
function param = interpolate_missing_params(convoice_list, param, angle_total)
    originalParam = param;
    startParam = 0;
    endParam = 0;

    for i = 1:length(convoice_list(:,1))
        if param(i,1) == -1
            if startParam == 0
                startParam = i;
            end
        end
        if ((param(i,1) ~= -1 && startParam ~= 0 && endParam == 0) || (param(i,1) == -1 && startParam ~= 0 && i == length(convoice_list(:,1))))
            endParam = i;
            if startParam == 1
                param(1,:) = [0, 0, 0, param(endParam, 4), 0, 0];
            else
                startParam = startParam - 1;
            end
            if endParam == angle_total
                param(angle_total,:) = [180, 0, 180, param(startParam, 4), 0, 0];
            end
            for j = startParam+1:endParam-1
                dist_S = j - startParam;
                dist_E = endParam - j;
                len = endParam - startParam;
                omomiS = dist_E / len;
                omomiE = dist_S / len;
                new_W = omomiS * param(startParam, 4) + omomiE * param(endParam, 4);
                new_B = omomiS * param(startParam, 5) + omomiE * param(endParam, 5);
                new_G = omomiS * param(startParam, 6) + omomiE * param(endParam, 6);
                param(j,:) = [(j-1)*10, originalParam(j,2), originalParam(j,3), new_W, new_B, new_G];
            end
            startParam = 0;
            endParam = 0;
        end

        if startParam == angle_total
            param(angle_total,:) = [180, 0, 180, param(angle_total-1, 4), 0, 0];
        end
    end
end

function save_parameters(data,N,person,dirc,prefix)
    f_G = polyfit(data(:,1), data(:,6), N);
    f_B = polyfit(data(:,1), data(:,5), N);
    f_W = polyfit(data(:,1), data(:,4), N);
    output_martix=[f_W.',f_B.',f_G.'];
    path=sprintf("parameters/%s",person);
    mkdir(path)
    save_path=sprintf("parameters/%s/%s_%s.txt",person,prefix,dirc);
    writematrix(output_martix,save_path);
end




function plot_param(data, N, x, person, dirc, prefix, suffix, originalData)
    if(suffix=="G")
        suf_ind=6;
    elseif(suffix=="B")
        suf_ind=5;
    else
        suf_ind=4;
    end
    f_plot = polyfit(data(:,1), data(:,suf_ind), N);
    y_plot = f_plot(N+1);
    for i = 1:N
        y_plot = y_plot + f_plot(i) * x.^(N+1-i);
    end
    
    f = figure;
    hold on
    scatter(data(:,1), data(:,suf_ind), 50, "DisplayName", sprintf('%s %s', prefix, suffix), "MarkerEdgeColor", "blue", "MarkerFaceColor", "blue");
    plot(x, y_plot, "LineWidth", 2, "DisplayName", "approximate curve");
    if exist('originalData', 'var')
        scatter(originalData(:,1), originalData(:,suf_ind), 40, "DisplayName", sprintf('original %s %s', prefix, suffix), "MarkerEdgeColor", "red", "MarkerFaceColor", "red");
    end
    legend;
    xlabel('elevation angle[degree]');
    ylabel(sprintf('%s [Hz]', suffix));
    xlim([1 (19-1)*10]);
    
    path = sprintf("./graph/%s/parameters/%s/", person, dirc);
    mkdir(path)
    savename = sprintf("%s%s_%s_%s.png", path, prefix, suffix, dirc);
    saveas(f, savename);
    close(f);
end
