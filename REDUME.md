- [1. 作業メモ](#1-作業メモ)
- [2. モデルでやりたいことリスト](#2-モデルでやりたいことリスト)
- [3. 行動方程式、上から計算](#3-行動方程式上から計算)
- [4. 定義式](#4-定義式)
  - [4.1. 名目値と実質値と価格の関係の恒等式](#41-名目値と実質値と価格の関係の恒等式)
- [5. 会計恒等式](#5-会計恒等式)
  - [5.1. 取引フロー表](#51-取引フロー表)
  - [5.2. バランスシート表](#52-バランスシート表)
  - [5.3. 完全統合表](#53-完全統合表)
  - [5.4. フローの整合性](#54-フローの整合性)
  - [5.5. ストックとフローの接続の整合性](#55-ストックとフローの接続の整合性)
  - [5.6. ストックの整合性](#56-ストックの整合性)
- [6. パラメータの値](#6-パラメータの値)

# 1. 作業メモ

- 新規株式発行と自社株買いを追加する。売上規模(名目値)と発行部数が比例する状態が望ましいかもしれない。企業の資金調達方法が追加される
  - 行動方程式を考える
  - シミュレーションに加える
- 政府部門内で資本ストックを作る。公的固定資本形成を作る
  - 表を更新する
  - 恒等式を更新する
  - 投資の行動方程式を考える
  - シミュレーションに加える
- 資本価格と消費財価格を分離する
  - 資本財価格の行動方程式を考える
  - シミュレーションに加える

# 2. モデルでやりたいことリスト

- 政府支出の増加率の変化の影響をシミュレーションしてまとめる
  - 技術や制度などの生産面の制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 消費需要の制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 企業の資金繰りの制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 政府の支出増加傾向の分類
    - 基準シナリオ（ほかの制約が影響しないことを確認）
    - 政府支出一定シナリオ
    - 政府支出急増シナリオ
- カレツキアンモデルの賃金主導型と利潤主導型が再現されているかどうか確認
  - 再現されないなら、それはなぜかを考察
    - 企業の内部資金がネックになって投資量が決まるときは利潤を増やすほど経済成長し、消費需要がネックになって投資量が決まるときは賃金を増やすほど経済成長する、というのが自分の直感ではある。内部資金と需要の一次結合で投資水準を決めるという仮定をカレツキアンモデルは採用しているようだが（要確認）カレツキアンモデルが示す性質はその仮定に決定的に依存しているのではないか？
      - $\epsilon_3$ を大きくすると賃金分配率が上がる。このパラメータで調整する
    - そもそも技術や制度に由来する一人当たりの生産能力自体がネックになっている場合は、利潤分配率や労働分配率と技術進歩の速度を関連付けない限り、利潤主導型とも賃金主導型とも異なる成長が起こるはず。
  - 賃金主導型成長だけ再現されるなら、その理由を考察
    - 賃金主導型成長が起こる条件を帰納的に探す
  - 利潤主導型成長だけ再現されるなら、その理由を考察
    - 利潤主導型成長が起こる条件を帰納的に探す
  - 現実社会が、利潤主導型と賃金主導型のどちらできれいに説明できるか、文献をあたる


# 3. 行動方程式、上から計算

- $\zeta_2 = \zeta_{2-1} (1 + \zeta_3 \cdot abs(randn()))$
- $p = p_{-1} \exp\{\mu_3 \min(1, \max(u^e_{-1}, \frac{C_{w-1}^D + C_{i-1}^D + G_{-1}^D}{p_{-1} ζ_{2-1}}) - u^T)\}$
- $W_f = (1 - \epsilon_1)W_{f-1} + \epsilon_1 \cdot \max(0, \epsilon_3(W_{f-1} + \Pi_{-1} - I_{-1}))$
- $W_b = (1 - \epsilon_1)W_{b-1} + \epsilon_1 \cdot \max(0, \epsilon_2(\Pi_{b-1} + r L_{-1}) + \epsilon_4 (H_{b-1} - M_{-1}))$
- $w_g = w_{g-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $g^D = g^D_{-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $C_w^D = ((1 - \alpha_5) C_{w-1}^D + \alpha_5 \cdot \max\{0, \alpha_1 (W^e-T_{iw}^e-T_{ew}^e-r L_{w-1}) + \alpha_2 (M_{w-1} + H_{w-1} - L_{w-1})\}) \frac{C_{w-1}}{C_{w-1}-\Delta C_{w-1}}$
- $C_i^D = ((1 - \alpha_5) C_{i-1}^D + \alpha_5 \cdot \max(0, \alpha_3 (\Pi_i^e-T_{ii}^e-T_{ei}^e) + (\alpha_4 + \alpha_6 \frac{E_{i-1}}{E_{-1}}) M_{i-1})) \frac{C_{w-1}}{C_{w-1}-\Delta C_{w-1}}$
- $C_w = \frac{C_w^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $C_i = \frac{C_i^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $G = \frac{G^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $T_{ef} = \gamma_1 K_{f-1}$
- $T_{ei} = \gamma_2 (M_{i-1} + E_{i-1})$
- $T_{ew} = \gamma_2 (M_{w-1} + H_{w-1})$
- $i = btw(0, (u_{-1} - u^T)k_{f-1} + \beta_1 k_{f-1}, \beta_2 \frac{M_{f-1} - L_{f-1}}{p})$
- $\Delta e = $ $c$ に比例する体系になるように
- $T_{ii} = \tau_1 \Pi_{i-1}$
- $T_{iw} = \tau_1 W_{-1}$
- $T_{ff} = \tau_2 (C+I+G-W_f-T_{ef}-r L_{f-1})$
- $\Pi_i = max(0, \frac{E_{i-1}}{E_{-1}}(\theta_1 (\Pi - I) + \theta_2 (M_{f-1} - L_{f-1})))$
- $\Pi_b = max(0, \frac{E_{b-1}}{E_{-1}}(\theta_1 (\Pi - I) + \theta_2 (M_{f-1} - L_{f-1})))$
- $T_{fb} = max(0, \tau_2 (\Pi_b + r L_{-1} - W_b + \Delta p_e e_{b-1}))$
- $H_w = \iota_1 C_w$
- $L_w = \iota_2 W$
- $E_i^T= \iota_3 NW_i^e$
- $E_b^T = (1 - \iota_7) E_{b-1} + \iota_7 (\frac{Π_i^e + Π_b^e}{E^e} \frac{L^e}{r})$
- $p_e = \frac{E_i^T + E_b^T}{e}$
- $E_i = min(\iota_3 NW_i, E)$
- $M_i = NW_i - E_i$
- $L_f^D = (1 - \iota_4)L_{f-1}^D + \iota_4 max(0, I + W_f + T_{ef} + T_{ff} + r L_{f-1} - M_{f-1} - NL_f)$
- $L_f = btw(0, L_f^D, max(0, \iota_5 NL_f))$
- $\Delta M_f = NL_f + \Delta L_f$ 新規株式発行や自社株買いを行わないという仮定

# 4. 定義式

- $x^e = ((1 - \lambda_e)*x^e[t-1] + \lambda_e*x[t-1])*x[t-1]/(x[t-1]-Δx[t-1])$
- $u = \frac{c + g}{\zeta_1 k_{f-1}}$
- $NL_w = -C_w+W-T_{iw}-T_{ew}-r L_{w-1}$
- $NL_i = -C_i-T_{ii}-T_{ei}+\Pi_i$
- $NL_f = -I + \Pi_f$
- $NL_b = -W_b - T_{fb} + \Pi_b + r L_{-1}$
- $NL_g = -G-W_g+T_i+T_e+T_f$
- $btw(A, B, C) = max(A, min(B, C))$ ただし( $A \leq C$ )

## 4.1. 名目値と実質値と価格の関係の恒等式

- $G = p g$
- $C_w = p c_w$
- $C_i = p c_i$
- $K_f = p k_f$
- $I = p i$
- $E_i = p_e e_i$
- $E = p_e e$
- $E_b = p_e e_b$

# 5. 会計恒等式

## 5.1. 取引フロー表

|                      |    労働者     |      資本家       |   企業(経常)    |  企業(資本)   |       銀行        | 統合政府(経常) | 統合政府(資本) | 合計 |
| :------------------- | :-----------: | :---------------: | :-------------: | :-----------: | :---------------: | :------------: | :------------: | ---- |
| 消費                 |    $-C_w$     |      $-C_i$       |      $+C$       |               |                   |                |                | $0$  |
| 投資                 |               |                   |      $+I$       |     $-I$      |                   |                |                | $0$  |
| 公的固定資本形成     |               |                   |                 |               |                   |                |                | $0$  |
| 政府支出（賃金除く） |               |                   |      $+G$       |               |                   |      $-G$      |                | $0$  |
| 賃金                 |     $+W$      |                   |     $-W_f$      |               |      $-W_b$       |     $-W_g$     |                | $0$  |
| 所得税               |   $-T_{iw}$   |     $-T_{ii}$     |                 |               |                   |     $+T_i$     |                | $0$  |
| 資産税               |   $-T_{ew}$   |     $-T_{ei}$     |    $-T_{ef}$    |               |                   |     $+T_e$     |                | $0$  |
| 法人税               |               |                   |    $-T_{ff}$    |               |     $-T_{fb}$     |     $+T_f$     |                | $0$  |
| 企業利潤             |               |     $+\Pi_i$      |     $-\Pi$      |   $+\Pi_f$    |     $+\Pi_b$      |                |                | $0$  |
| 借入金金利           | $-r L_{w-1}$  |                   |  $-r L_{f-1}$   |               |    $+r L_{-1}$    |                |                | $0$  |
| [メモ:Net Lending]   |    $NL_w$     |      $NL_i$       |                 |    $NL_f$     |      $NL_b$       |                |     $NL_g$     | $0$  |
| 預金                 | $-\Delta M_w$ |   $-\Delta M_i$   |                 | $-\Delta M_f$ |    $+\Delta M$    |                |                | $0$  |
| 借入                 | $+\Delta L_w$ |                   |                 | $+\Delta L_f$ |    $-\Delta L$    |                |                | $0$  |
| 株式                 |               | $-p_e \Delta e_i$ | $+p_e \Delta e$ |               | $-p_e \Delta e_b$ |                |                | $0$  |
| 現金                 | $-\Delta H_w$ |                   |                 |               |   $-\Delta H_b$   |                |  $+\Delta H$   | $0$  |
| 合計                 |      $0$      |        $0$        |       $0$       |      $0$      |        $0$        |      $0$       |      $0$       |      |

## 5.2. バランスシート表

|        | 労働者  | 資本家  |  企業   |  銀行   | 統合政府 | 合計 |
| :----- | :-----: | :-----: | :-----: | :-----: | :------: | :--: |
| 資本   |         |         |  $+K_f$   |         |  $+K_g$  | $+K$ |
| 預金   | $+M_w$  | $+M_i$  | $+M_f$  |  $-M$   |          | $0$  |
| 株式   |         | $+E_i$  |  $-E$   | $+E_b$  |          | $0$  |
| 借入   | $-L_w$  |         | $-L_f$  |  $+L$   |          | $0$  |
| 現金   | $+H_w$  |         |         | $+H_b$  |   $-H$   | $0$  |
| 純資産 | $-NW_w$ | $-NW_i$ | $-NW_f$ | $-NW_b$ | $-NW_g$  | $-K$ |
| 合計   |   $0$   |   $0$   |   $0$   |   $0$   |   $0$    | $0$  |

## 5.3. 完全統合表

|                        |    労働者     |           資本家            |            企業            |            銀行             |         統合政府          |             合計             |
| :--------------------- | :-----------: | :-------------------------: | :------------------------: | :-------------------------: | :-----------------------: | :--------------------------: |
| 期首純資産             |  $NW_{w-1}$   |         $NW_{i-1}$          |         $NW_{f-1}$         |         $NW_{b-1}$          |        $NW_{g-1}$         |           $K_{-1}$           |
| 資本のキャピタルゲイン |               |                             | $+\Delta p \cdot k_{f-1}$  |                             | $+\Delta p \cdot k_{g-1}$ |     $+\Delta p \cdot k$      |
| 資本の増減             |               |                             |    $+p \cdot \Delta k_f$     |                             |   $+p \cdot \Delta k_g$   | $+p \cdot \Delta k$ |
| 預金の増減             | $+\Delta M_w$ |        $+\Delta M_i$        |       $+\Delta M_f$        |         $-\Delta M$         |                           |             $0$              |
| 株式のキャピタルゲイン |               | $+\Delta p_e \cdot e_{i-1}$ | $-\Delta p_e \cdot e_{-1}$ | $+\Delta p_e \cdot e_{b-1}$ |                           |             $0$              |
| 株式の増減             |               |   $+p_e \cdot \Delta e_i$   |   $+p_e \cdot \Delta e$    |   $+p_e \cdot \Delta e_b$   |                           |             $0$              |
| 借入の増減             | $-\Delta L_w$ |                             |       $-\Delta L_f$        |         $+\Delta L$         |                           |             $0$              |
| 現金の増減             | $+\Delta H_w$ |                             |                            |        $+\Delta H_b$        |        $-\Delta H$        |             $0$              |
| 期末純資産             |    $NW_w$     |           $NW_i$            |           $NW_f$           |           $NW_b$            |          $NW_g$           |             $K$              |

## 5.4. フローの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $-C_w+W-T_{iw}-T_{ew}-r L_{w-1} = \Delta M_w - \Delta L_w + \Delta H_w$
- [x] $-C_i-T_{ii}-T_{ei}+P_i = \Delta M_i + p_e \Delta e_i$
- [x] $\Pi = C+I+G-W_f-T_{ef}-T_{ff}-r L_{f-1}$
- [x] $-I + \Pi_f = \Delta M_f - \Delta L_f - \Delta E$
- [x] $-W_b - T_{fb} + \Pi_b + r L_{-1} = -\Delta M + \Delta L + p_e \Delta e_b + \Delta H_b$
- [ ] $-G-W_g+T_i+T_e+T_f = -\Delta H$
- [x] $C = C_w + C_i$
- [x] $W = W_f + W_g + W_b$
- [x] $T_i = T_{iw} + T_{ii}$
- [x] $T_e = T_{ew} + T_{ei} + T_{ef}$
- [x] $T_f = T_{ff} + T_{fb}$
- [x] $\Pi_f = \Pi - \Pi_i - \Pi_b$
- [ ] $NL_w + NL_i + NL_f + NL_b + NL_g = 0$
- [x] $\Delta M = \Delta M_w + \Delta M_i + \Delta M_f$
- [x] $\Delta L = \Delta L_w + \Delta L_f$
- [ ] $\Delta e_i + \Delta e_b - \Delta e = 0$
- [x] $\Delta H = \Delta H_w + \Delta H_b$

## 5.5. ストックとフローの接続の整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $k_f = (1 - \beta_1)k_{f-1} + i$
- [x] $\Delta k = k - k_{-1}$
- [x] $K_f = K_{f-1} + \Delta K_f = p_{-1} k_{f-1} + \Delta p \cdot k_{f-1} + p \cdot \Delta k_f$
- [x] $\Delta M_w = M_w - M_{w-1}$
- [x] $\Delta M_i = M_i - M_{i-1}$
- [x] $M_f = M_{f-1} + \Delta M_f$
- [x] $M = M_{-1} + \Delta M$
- [x] $\Delta L_w = L_w - L_{w-1}$
- [x] $\Delta L_f = L_f - L_{f-1}$
- [ ] $L = L_{-1} + \Delta L$
- [x] $E_i = E_{i-1} + \Delta E_i = p_{e-1} e_{i-1} + \Delta p_e \cdot e_{i-1} + p_e \cdot \Delta e_i$
- [x] $E_b = E_{b-1} + \Delta E_b = p_{e-1} e_{b-1} + \Delta p_e \cdot e_{b-1} + p_e \cdot \Delta e_b$
- [x] $E = E_{-1} + \Delta E = p_{e-1} e_{-1} + \Delta p_e \cdot e_{-1} + p_e \cdot \Delta e$
- [x] $\Delta H_w = H_w - H_{w-1}$
- [x] $H_b = H_{b-1} + \Delta H_b$
- [x] $H = H_{-1} + \Delta H$
- [x] $NW_w = NW_{w-1} + \Delta M_w - \Delta L_w + \Delta H_w (=NW_{w-1} + NL_w)$
- [x] $NW_i = NW_{i-1} + \Delta M_i + \Delta p_e \cdot e_{i-1}  + p_e \cdot \Delta e_i (= NW_{i-1} + NL_i + \Delta p_e \cdot e_{i-1})$
- [x] $NW_f = NW_{f-1} + \Delta K_f + \Delta M_f - \Delta E - \Delta L_f (= NW_{f-1} + NL_f - \Delta E + \Delta K_f)$
- [x] $NW_b = NW_{b-1} - \Delta M + \Delta E_b + \Delta L + \Delta H_b (= NW_{b-1} + NL_b + \Delta p_e \cdot e_{b-1})$
- [x] $NW_g = NW_{g-1} - \Delta H (= NW_{g-1} + NL_g)$

## 5.6. ストックの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $M = M_w + M_i + M_f$
- [ ] $E = E_i + E_b$
- [x] $e = e_i + e_b$
- [x] $L = L_w + L_f$
- [x] $H = H_w + H_f$
- [ ] $NW_w + NW_i + NW_f + NW_b + NW_g = 0$
- [ ] $NW_w = M_w - L_w + H_w$
- [ ] $NW_i = M_i + E_i$
- [ ] $NW_f = K_f + M_f - E - L_f$
- [ ] $NW_b = -M + E_b + L + H_b$
- [ ] $NW_g = -H$

# 6. パラメータの値

- $\alpha_1 = 0.9$
- $\alpha_2 = 0.1$
- $\alpha_3 = 0.1$
- $\alpha_4 = 0.02$
- $\alpha_5 = 0.5$
- $\alpha_6 = 0.2$
- $\beta_1 = 0.05$
- $\beta_2 = 0.5$
- $\gamma_1 = 0.015$
- $\gamma_2 = 0.02$
- $\delta_1 = 0.02$
- $\delta_2 = 0.3$
- $\epsilon_1 = 0.5$
- $\epsilon_2 = 0.7$
- $\epsilon_3 = 0.8$
- $\epsilon_4 = 0.05$
- $\zeta_1 = 1.0$
- $\zeta_2 = 300.0$ (シミュレーションのための初期値)
- $\zeta_3 = 0.023$
- $\theta_1 = 0.4$
- $\theta_2 = 0.1$
- $\iota_1 = 0.1$
- $\iota_2 = 1.0$
- $\iota_3 = 0.9$
- $\iota_4 = 0.5$
- $\iota_5 = 10.0$
- $\iota_6 = 0.5$
- $\iota_7 = 0.5$
- $\lambda_e = 0.5$
- $\mu_1 = 0.3$
- $\mu_2 = 0.1$
- $\mu_3 = 0.5$
- $\tau_1 = 0.2$
- $\tau_2 = 0.2$
- $u^T = 0.8$
- $G_0 = 1$
- $r = 0.01$

