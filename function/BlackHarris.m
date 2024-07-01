function [output] = BlackHarris(HRIR,N)
%earlyHRTFのためにHRIRを処理する関数
%1:振幅の最大値のインデックスを取り出す
%2:そこを中心として２Nサンプルのブラックマンハリス窓で切り抜く
%3:512サンプルのゼロ行列の中心に窓ので切り抜いた配列の中心を合わせて加算する
 
%振幅の最大値を見つける
[~,Index_max]=max(abs(HRIR));
%最大値のインデックスを中心に2Nサンプル切り抜く
target=HRIR(Index_max-N+1:Index_max+N);
%ブラックマンハリス窓を作る
window=blackmanharris(2*N);
%窓とかけ合わせる
produced=target.*window;
%出力用の配列を用意する
output=zeros(512,1);
output(257-N+1:257+N)=produced;
end