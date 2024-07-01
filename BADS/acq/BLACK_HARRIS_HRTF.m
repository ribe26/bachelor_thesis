function HRTF = BLACK_HARRIS_HRTF(HRIR, Ns)
    %%Ns:ブラックマンハリス窓のポイント数, HRIR:頭部インパルス応答
    
    [HRIR_MAX, Index_MAX] = max(abs(HRIR));
    
    HRIR_1 = HRIR(Index_MAX-Ns+1:Index_MAX+Ns);
    w = blackmanharris(2*Ns);  %ブラックマンハリス窓の作成
    HRIR_2 = HRIR_1 .* w;      %ブラックマンハリス窓で振幅スペクトルを切り出す     

    HRIR_New = zeros(512,1);   %512サンプルの空配列を用意

    HRIR_New(257-Ns+1:257+Ns) = HRIR_2;

    HRTF = abs(fft(HRIR_New));
end