function J = JE1_W(P, HF, Fs)
NFFT = 512;
if mod(NFFT,2)==0
      Nout = (NFFT/2)+1;
else
      Nout = (NFFT+1)/2;
end
f = ((0:Nout-1)'./NFFT).*Fs;
[Hd,~] = PK(P, Fs);

disp(length(Hd(2:end)));
disp(length(HF(2:end)));
disp(length(f(2:end)));

E = sum(abs(Hd(2:end) - HF(2:end)) ./ (2*pi*f(2:end)./Fs));
%E = sum(abs(Hd(2:end) - HF(2:end)));
J = mean(E);
end