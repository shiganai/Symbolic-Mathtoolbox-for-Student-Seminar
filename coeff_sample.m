
% 文字と関数を用意する
syms x y
f = x^2 + y + x*y + 2*x + 3*x

% coeffs によって係数を取得できる
[Coeffs, Terms] = coeffs(f, [x,y])

% Coeffs .* Terms の総和が元の表現に等しくなる
isequal( f, sum( Coeffs .* Terms ) )

% 入力によって取り出す変数を選べる
[Coeffs, Terms] = coeffs(f, [x])

% 関数ベクトルに対しては適応できない
f_vec = [
    x^2 + y + x*y + 2*x + 3*x;
    x^3 + x^2*y + y^2 + 2*x
    ]
try
    [Coeffs, Terms] = coeffs(f_vec, [x,y])
catch ME
    ME.message
end

% 関数ベクトルに対して定義できないのは不便なので作った
[Coeffs, Terms] = coeffs_Vector( f_vec, [x,y])
isequal( f_vec, sum( Coeffs .* Terms, 2 ) )


















