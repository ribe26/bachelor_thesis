%HRTFを双二次フィルタの重ね合わせで近似するコード
%最適化部分はオープンソースのベイズ最適化コードを拝借した
clear;
addpath('./function');
addpath('./BADS')


person_list=["sato"];
elev=[0,10,20,30,40,50,60,70,80,90];
azim=[0,180];
dirc_list=["R","L"];

Fs=48000;
NFFT=512;
Nout=NFFT/2+1;
f = ((0:Nout-1)'./NFFT).*Fs;


for peo=1:length(person_list)
    person=person_list(peo);
    for direction=1:2
        dirc=dirc_list(direction);
        for elevation=1:length(elev)
            for az=1:length(azim)
                %IRfname=new_hrtfpath(person,elev(elevation),azim(az),dirc);
                IRfname=sprintf("./hrtf/%s/elev%d/%s%de%03da_new.dat",person,elev(elevation),dirc,elev(elevation),azim(az));
                Rdat = fopen(IRfname,'r','b');
                IRdata = fread(Rdat,'float');
                IR_spec=fft(IRdata,NFFT);

                [ely,ely_abs]=earlyHRTF(IRdata,66,NFFT);
                BF=20*log10(ely_abs(1:Nout));

                [peaks,peak_index]=findpeaks(20*log10(ely_abs(1:length(ely_abs)/2)),"MinPeakProminence",3);
                [notch,notch_index]=findpeaks(-20*log10(ely_abs(1:length(ely_abs)/2)),"MinPeakProminence",2,"MinpeakHeight",0);
                notch=notch*-1;



                peak_index=peak_index(f(peak_index)>2000);
                peak_index=peak_index(f(peak_index)<25000);
                labelp=ones(length(peak_index),1);

                notch_index=notch_index(f(notch_index)>4000);
                notch_index=notch_index(f(notch_index)<25000);
                labeln=zeros(length(notch_index),1)-ones(length(notch_index),1);


                pnum=1;
                nnum=3;
                if length(peak_index)<pnum
                    pnum=length(peak_index);
                end
                if length(notch_index)<nnum
                    nnum=length(notch_index);
                end


                %% ---------------------ノッチとピークを結合-------------------------
                %ピークとノッチにラベルをつけたからピーク，ノッチが交互に並ぶようにする
                index = zeros(pnum + nnum, 2);
                i = 1; j = 1; k = 1;
                while 1
                    if j <= pnum
                        index(i,:) = [peak_index(j), labelp(j)];
                        i = i + 1;
                        j = j + 1;
                    end
                    if k <= nnum
                        index(i,:) = [notch_index(k), labeln(k)];
                        i = i + 1;
                        k = k + 1;
                    end
                    if i > pnum + nnum
                        break;
                    end
                end
                IDTable = array2table(index);
                IDTable = sortrows(IDTable);
                index = table2array(IDTable);

                %% フィルタ最適化
                % P     :パラメータ[角周波数w, バンド幅B, ゲインG]の保存用
                Q_l = 2*pi*100/Fs;  %バンド幅Bの下限値(ベイズ最適化)
                Q_h = 2*pi*10000/Fs; %バンド幅Bの上限値(ベイズ最適化)

                NPK=pnum+nnum;
                Nout=NFFT/2+1;
                %f =(0:Fs/2-1);
                f = ((0:Nout-1)'./NFFT).*Fs;
                %----------------------PKのパラメータ最適化-----------------------
                for N = 1 : NPK
                    if NPK <= length(index)
                        %---------------------初期値を設定-----------------------
                        Freq = f(index(N,1));
                        Pin = [2*pi*Freq/Fs, 2*pi*1000/Fs, BF(index(N,1))];
                        w_l = 2*pi*Freq/Fs;
                        w_h = 2*pi*Freq/Fs;
                        %----------------------最適化部分-----------------------
                        if N == 1
                            E1 = @(P)JE1_W(P, BF, Fs);
                            G_l = BF(index(N,1));
                            G_h = BF(index(N,1));
                            lb = [w_l, Q_l, G_l];   %badsとランダム探索に渡す下限値の配列
                            ub = [w_h, Q_h, G_h];   %badsとランダム探索に渡す上限値の配列
                            PK_best = bads(E1, Pin, lb, ub);
                            [~, sos] = PK(PK_best, Fs);
                            P = [index(N,2), PK_best];
                        else
                            E2 = @(P)JE2_W(P, BF, sos, Fs);
                            G_l = -80;
                            G_h = 30;
                            lb = [w_l, Q_l, G_l];   %badsとランダム探索に渡す下限値の配列
                            ub = [w_h, Q_h, G_h];   %badsとランダム探索に渡す上限値の配列
                            PK_best = bads(E2, Pin, lb, ub);

                            %あまりにもノッチが低すぎる場合は調節する。
                            %{
                            for i=1:length(PK_best(:,1))
                                if(PK_best(i,3)<-30)
                                    PK_best(i,3)=20*log10(abs(IR_spec(index(N,1))));
                                end
                            end
                            %}



                            [~, sos_PK] = PK(PK_best, Fs);
                            sos = [sos; sos_PK];
                            P = [P ;index(N,2), PK_best];
                        end
                    else
                        flag = 1;
                    end
                    clear sos_PK PK lb ub Hd;
                end


                TABLE = [P, sos];
                TABLE = array2table(TABLE);
                TABLE = sortrows(TABLE,2);
                P = table2array(TABLE(1:NPK,1:4));
                sos = table2array(TABLE(1:NPK,5:10));

                %----------------------------全体の最適化------------------------
                G_l = -80;
                G_h = 30;
                0;
                for j = 1 : NPK - 3
                    for n = j : j + 3
                        E_fin = @(P)JE4_W(P, BF, sos, n, Fs);
                        P_fin = P(n,2:4);
                        w_l = P_fin(1);
                        w_h = P_fin(1);
                        lb = [w_l, Q_l, G_l];     %bads関数に渡す下限値の配列
                        ub = [w_h, Q_h, G_h];     %bads関数に渡す上限値の配列
                        P(n,2:4) = bads(E_fin, P_fin, lb, ub);

                        %あまりにもノッチが低すぎる場合は調節する。
                        %{
                            for i=1:length(P(:,1))
                                if(P(i,4)<-30)
                                    P(i,4)=20*log10(abs(IR_spec(index(n,1))));
                                end
                            end
                        %}

                        [~, sos_PK] = PK(P(n,2:4), Fs);
                        sos(n,:) = sos_PK;
                        clear sos_PK lb ub;
                    end
                end
                %Pに最終的な中心周波数、バンド幅、ゲインが入っている
                %P=[?,w1,B1,G1]
                %  [?,w2,B2,G2]...

                %------------------------周波数でソートを行う----------------------
                TABLE = array2table(P);
                TABLE = sortrows(TABLE);
                POUT = table2array(TABLE(1:NPK,1:4));
                %--------------------最終的な係数からフィルタ作成--------------------

                %Hd = freqz(sos,Nout);
                %Hd = 20*log10(abs(Hd));
                P_out_matrix=zeros(length(P(:,1)),3);
                for loop=1:length(P(:,1))
                    P_out_matrix(loop,:)=P(loop,2:4);
                end
                %P_out_matrix=[P(1,2:4);P(2,2:4);P(3,2:4);P(4,2:4);P(5,2:4);P(6,2:4)];
                save_path=sprintf("nums/%s/",person);
                mkdir(save_path);

                save_path=sprintf("nums/%s/nums_%d_%d_%s.txt",person,elev(elevation),azim(az),dirc);
                writematrix(P_out_matrix,save_path);
                save_path=sprintf("nums/%s/sos_%d_%d_%s.txt",person,elev(elevation),azim(az),dirc);
                writematrix(sos,save_path);
            end
        end
    end
end


%最適化結果の一部をプロットする
%{
f1=figure
hold on
for i = 1 : NPK
                    semilogx(f(1:end), 20*log10(freqz(sos(i,1:3),sos(i,4:6),Nout)));
                    axis([90 Fs/2 -60 20]);
                    xlabel('Frequency [Hz]');
                    ylabel('Relative level [dB]');
end
                
TF=ones(Nout,1);
for i = 1:length(sos(:,1))
    TF=TF.*freqz(sos(i,1:3),sos(i,4:6),Nout);
end



f2=figure
hold on
plot(f,20*log10(abs(TF)),'DisplayName','parametric early HRTF');
plot(f,BF,'DisplayName','early HRTF');
legend
axis([90 Fs/2 -60 20]);
xlabel('Frequency [Hz]');
ylabel('Relative level [dB]');
                

f3=figure
hold on
plot(f,20*log10(abs(IR_spec(1:length(IR_spec)/2))),'DisplayName','HRTF');
plot(f,BF(2:Fs/2+1),'DisplayName','early HRTF');
legend
axis([90 Fs/2 -60 20]);
xlabel('Frequency [Hz]');
ylabel('Relative level [dB]');
%}

                





