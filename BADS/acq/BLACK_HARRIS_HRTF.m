function HRTF = BLACK_HARRIS_HRTF(HRIR, Ns)
    %%Ns:�u���b�N�}���n���X���̃|�C���g��, HRIR:�����C���p���X����
    
    [HRIR_MAX, Index_MAX] = max(abs(HRIR));
    
    HRIR_1 = HRIR(Index_MAX-Ns+1:Index_MAX+Ns);
    w = blackmanharris(2*Ns);  %�u���b�N�}���n���X���̍쐬
    HRIR_2 = HRIR_1 .* w;      %�u���b�N�}���n���X���ŐU���X�y�N�g����؂�o��     

    HRIR_New = zeros(512,1);   %512�T���v���̋�z���p��

    HRIR_New(257-Ns+1:257+Ns) = HRIR_2;

    HRTF = abs(fft(HRIR_New));
end