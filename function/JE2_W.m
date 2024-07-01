function J = JE2_W(P, HF, sos, Fs)
NFFT = 512;
if mod(NFFT,2)==0
      Nout = (NFFT/2)+1;
else
      Nout = (NFFT+1)/2;
end
f = ((0:Nout-1)'./NFFT).*Fs;
[~, sos_new] = PK(P, Fs);
sos = [sos; sos_new];
Hd = freqz(sos,Nout);
Hd = 20*log10(abs(Hd));

E = sum(abs(Hd(2:end) - HF(2:end)) ./ (2*pi*f(2:end)./Fs));
%E = sum(abs(Hd(2:end) - HF(2:end)));
J = mean(E);
end