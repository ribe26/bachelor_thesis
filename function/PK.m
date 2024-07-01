function [Hd, sos] = PK(P, Fs)
%ゲイン，バンド幅,中心周波数から双二次フィルタの係数を求める
NFFT = 512;
if mod(NFFT,2)==0
    Nout = (NFFT/2)+1;
else
    Nout = (NFFT+1)/2;
end
w = P(1);
B = P(2);
G = P(3);
alpha = sin(w)*sinh((log(2)/2)*B*(w/(sin(w)))); g = 10^(G/40);
bL0 = 1 + alpha * g;
bL1 = -2 * cos(w);
bL2 = 1 - alpha * g;
aL0 = 1 + alpha / g;
aL1 = -2 * cos(w);
aL2 = 1 - alpha / g;
b = [bL0 bL1 bL2]; a = [aL0 aL1 aL2];
Hd = freqz(b,a,Nout);
Hd = 20*log10(abs(Hd));
sos = [b, a];
end