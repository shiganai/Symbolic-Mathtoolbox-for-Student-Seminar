
% t をパラメータとする角度 theta と独立な距離 r を定義する
% 仮定は宣言できないことは declare_sample.m を参照
syms theta(t)
syms r real

% 位置ベクトルを横ベクトルで定義する
position = r * [cos(theta), sin(theta)]

% 速度ベクトルは位置ベクトルの微分で計算できる
% diff( position ) でも diff( position, t )でも構わない
velocity = diff( position )
isequal_diff_vs_difft = isequal( diff( position ), diff( position, t ) )

% 速度ベクトルの2乗を計算する。(r\omega)^2を期待する
% しかしそうはならない
simplify( velocity * velocity' )

% conj(theta) が邪魔
% パラメータを持たない文字として別の文字を定義し、置き換える
% 微分も正しく計算されているのが分かる
syms theta_R omega_R real
subs_R = @(input) subs( input, [theta, diff(theta,t)], [theta_R, omega_R])
velocity_R = subs_R( velocity )
simplify( subs_R( velocity * velocity' ) )

% 加速度ベクトルを定義して、そのあとに置き換える
% diff( X, t, t ) によって2階微分が計算できる
acceleration = diff( position, t, t )
syms alpha_R real
subs_R = @(input) subs( input, [theta, diff(theta,t), diff(theta,t,t)], [theta_R, omega_R, alpha_R])
acceleration_R = simplify( subs_R( acceleration ) )

% 動径方向と回転方向に分ける
% ベクトルを分けるのは1組の直行ベクトルとの内積が良い
r_vec = [cos(theta), sin(theta)];
rot_vec = [-sin(theta), cos(theta)];

acceleration_r = acceleration * r_vec';
acceleration_rot = acceleration * rot_vec';

acceleration_r_R = simplify( subs_R( acceleration_r ) )
acceleration_rot_R = simplify( subs_R( acceleration_rot ) )

% 端点における質量と端点での力を動径方向と回転方向に定義すれば
% 円運動の方程式が導ける
syms m F_r F_rot real
Netwon_eq = [
    m * acceleration_r_R == F_r
    m * acceleration_rot_R == F_rot
    ]

% ちなみに、Netwon_eq(t) の各項目にアクセスするためには Netwon_eq(2,1) は正しくない
try
    Netwon_eq(2,1)
catch ME
    ME.message
end

% Netwon_eq(t) の各項目にアクセスするためにはまず formula を取る
Netwon_eq_f = formula( Netwon_eq )
Netwon_eq_f(2,1)

%% 合成関数の微分について

% パラメータに対して1回しか経由しない場合は問題ない
syms x(t)
f = x^2;
diff(f,x)
diff(f,t)

% パラメータに対して2回またぐと0になる
syms y(x)
fy = y^2;
diff(fy,y)
diff(fy,x)
diff(fy,t)





























