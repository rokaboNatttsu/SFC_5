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

- 資本財生産企業と消費財生産企業を分ける
- 企業から統合政府部門への支払いは、税金以外も追加するべきか？
- 新規株式発行と自社株買いを追加する。売上規模(名目値)と発行部数が比例する状態が望ましいかもしれない。企業の資金調達方法が追加される
  - 行動方程式を考える
  - シミュレーションに加える
- 政府部門内で資本ストックを作る。公的固定資本形成を作る
  - 投資の行動方程式を考える
  - シミュレーションに加える
- 資本価格と消費財価格を分離する
  - 資本財価格の行動方程式を考える
  - シミュレーションに加える
- 家計の人数を導入
- パラメータの整理

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
- $\zeta_4 = \zeta_{4-1} (1 + \zeta_3 \cdot abs(randn()))$
- $p = p_{-1} \exp\{\mu_3 \min(1, \max(u_{c-1}^e, \frac{C_{w-1}^D + C_{i-1}^D + G_{-1}^D}{p_{-1} \zeta_{2-1}}) - u^T)\}$
- $p_k = p_{k-1} \exp\{\mu_3 \min(1, \max(u_{k-1}^e, \frac{I_{f-1}^D + I_{g-1}^D}{p_{k-1} \zeta_{4-1}}) - u^T)\}$
- $W_c = (1 - \epsilon_1)W_{c-1} + \epsilon_1 \cdot \max(0, \epsilon_3(W_{c-1} + \Pi_{c-1}))$
- $W_k = (1 - \epsilon_1)W_{k-1} + \epsilon_1 \cdot \max(0, \epsilon_3(W_{k-1} + \Pi_{k-1}))$
- $W_b = (1 - \epsilon_1)W_{b-1} + \epsilon_1 \cdot \max(0, \epsilon_2(\Pi_{b-1} + r L_{-1}) + \epsilon_4 (H_{b-1} - M_{-1}))$
- $w_g = w_{g-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $g^D = g^D_{-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $C_w^D = ((1 - \alpha_5) C_{w-1}^D + \alpha_5 \cdot \max\{0, \alpha_1 (W^e-T_{iw}^e-T_{ew}^e-r L_{w-1}) + \alpha_2 (M_{w-1} + H_{w-1} - L_{w-1})\}) \frac{C_{w-1}}{C_{w-1}-\Delta C_{w-1}}$
- $C_i^D = ((1 - \alpha_5) C_{i-1}^D + \alpha_5 \cdot \max(0, \alpha_3 (\Pi_{ci}^e+\Pi_{ki}^e-T_{ii}^e-T_{ei}^e) + (\alpha_4 + \alpha_6 (\frac{E_{ci-1}}{E_{c-1}}+\frac{E_{ki-1}}{E_{k-1}})) M_{i-1})) \frac{C_{w-1}}{C_{w-1}-\Delta C_{w-1}}$
- $C_w = \frac{C_w^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $C_i = \frac{C_i^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $G = \frac{G^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{f-1}, \zeta_2)})}$
- $T_{ec} = \gamma_1 K_{c-1}$
- $T_{ek} = \gamma_1 K_{k-1}$
- $T_{ei} = \gamma_2 (M_{i-1} + E_{i-1})$
- $T_{ew} = \gamma_2 (M_{w-1} + H_{w-1})$
- $i_c^D = btw(0, (u_{c-1} - u^T)k_{c-1} + \beta_1 k_{c-1}, \beta_2 \frac{M_{c-1} - L_{c-1}}{p})$
- $i_k^D = btw(0, (u_{k-1} - u^T)k_{k-1} + \beta_1 k_{k-1}, \beta_2 \frac{M_{k-1} - L_{k-1}}{p_k})$
- $i_g^D = i_{g-1}(1 + \delta_1 - \delta_2 \frac{p_k-p_{k-1}}{p_{k-1}})$
- $i_c = \frac{i_c^D}{\max(1, \frac{i_f^D + i_g^D}{p \cdot \min(\zeta_5 k_{f-1}, \zeta_4)})}$
- $i_k = \frac{i_k^D}{\max(1, \frac{i_f^D + i_g^D}{p \cdot \min(\zeta_5 k_{f-1}, \zeta_4)})}$
- $i_g = \frac{i_g^D}{\max(1, \frac{i_f^D + i_g^D}{p \cdot \min(\zeta_5 k_{f-1}, \zeta_4)})}$
- $\Delta e = \max(0, \frac{\iota_8 (I_f-\beta_1 K_{f-1})}{p_{e-1}}) - \max(0, \frac{\iota_9 (M_{f-1} - L_{f-1})}{p_{e-1}} - \iota_{10} c)$
- $T_{ii} = \tau_1 \Pi_{i-1}$
- $T_{iw} = \tau_1 W_{-1}$
- $T_{fc} = \tau_2 (C+G-W_c-T_{ec}-r L_{c-1})$
- $T_{fk} = \tau_2 (I_f+I_g-W_k-T_{ek}-r L_{k-1})$
- $\Pi_{ci} = max(0, \frac{E_{ci-1}}{E_{c-1}}(\theta_1 (\Pi_c - I_c) + \theta_2 (M_{c-1} - L_{c-1})))$
- $\Pi_{cb} = max(0, \frac{E_{cb-1}}{E_{c-1}}(\theta_1 (\Pi_c - I_c) + \theta_2 (M_{c-1} - L_{c-1})))$
- $\Pi_{ki} = max(0, \frac{E_{ki-1}}{E_{k-1}}(\theta_1 (\Pi_k - I_k) + \theta_2 (M_{k-1} - L_{k-1})))$
- $\Pi_{kb} = max(0, \frac{E_{kb-1}}{E_{k-1}}(\theta_1 (\Pi_k - I_k) + \theta_2 (M_{k-1} - L_{k-1})))$
- $T_{fb} = max(0, \tau_2 (\Pi_{cb} + \Pi_{kb} + r L_{-1} - W_b + \Delta p_{ec} e_{cb-1} + \Delta p_{ek} e_{kb-1}))$
- $H_w = \iota_1 C_w$
- $L_w = \iota_2 W$
- $E_{ci}^T= \iota_3 NW_i^e \frac{C+G}{C+G+I_c+I_g}$
- $E_{ki}^T= \iota_3 NW_i^e \frac{I_c+I_g}{C+G+I_c+I_g}$
- $E_{cb}^T = ((1 - \iota_7) E_{cb-1} + \iota_7 \frac{r_{E-1}}{r_{-1} + r_{E-1}} (L + E_{cb-1} + E_{kb-1}))\frac{E_{c-1}}{E_{c-1}+E_{k-1}}$
- $E_{kb}^T = ((1 - \iota_7) E_{kb-1} + \iota_7 \frac{r_{E-1}}{r_{-1} + r_{E-1}} (L + E_{cb-1} + E_{kb-1}))\frac{E_{k-1}}{E_{c-1}+E_{k-1}}$
- $p_{ec} = \frac{E_{ci}^T + E_{cb}^T}{e_c}$
- $p_{ek} = \frac{E_{ki}^T + E_{kb}^T}{e_k}$
- $E_{ci} = min(\iota_3 NW_i \frac{C+G}{C+G+I_c+I_g}, E_c)$
- $E_{ki} = min(\iota_3 NW_i \frac{I_c+I_g}{C+G+I_c+I_g}, E_k)$
- $M_i = NW_i - E_{ci} - E_{ki}$
- $L_c^D = (1 - \iota_4)L_{c-1}^D + \iota_4 max(0, \iota_{11} (I_c + W_c + T_{ec} + T_{fc} + r L_{c-1}) - M_{c-1})$
- $L_c = btw(0, L_c^D, max(0, \iota_5 NL_c))$
- $L_k^D = (1 - \iota_4)L_{k-1}^D + \iota_4 max(0, \iota_{11} (W_k+T_{ek}+T_{fk}+r L_{k-1}) - M_{k-1})$
- $L_k = btw(0, L_k^D, max(0, \iota_5 NL_k))$
- $\Delta M_c = NL_c + \Delta L_c + p_{ec} \cdot \Delta e_c$
- $\Delta M_k = NL_k + \Delta L_k + p_{ek} \cdot \Delta e_k$

# 4. 定義式

- $x^e = ((1 - \lambda_e) x^e_{-1} + \lambda_e x_{-1}) \frac{x_{-1}}{x_{-1}-\Delta x_{-1}}$
- $u_c = \frac{c + g}{\zeta_1 k_{f-1}}$
- $u_k = \frac{i_f}{\zeta_1 k_{f-1}}$
- $NL_w = -C_w+W-T_{iw}-T_{ew}-r L_{w-1}$
- $NL_i = -C_i-T_{ii}-T_{ei}+\Pi_i$
- $NL_f = -I_f + \Pi_f$
- $NL_b = -W_b - T_{fb} + \Pi_b + r L_{-1}$
- $NL_g = -I_g + GS$
- $btw(A, B, C) = max(A, min(B, C))$ ただし( $A \leq C$ )
- $r_E = \frac{\Pi_{ci} + \Pi_{cb} + \Pi_{ki} + \Pi_{kb}}{E_c + E_k}$

## 4.1. 名目値と実質値と価格の関係の恒等式

- $G = p g$
- $C_w = p c_w$
- $C_i = p c_i$
- $I_f = p_k i_f$
- $I_c = p_k i_c$
- $I_k = p_k i_k$
- $I_g = p_k i_k$
- $K_c = p_k k_c$
- $K_k = p_k k_k$
- $K_g = p_k k_g$
- $E_i = p_e e_i$
- $E = p_e e$
- $E_b = p_e e_b$

# 5. 会計恒等式

## 5.1. 取引フロー表

|                      |    労働者     |         資本家          | 消費財生産企業(経常) | 消費財生産企業(資本) | 資本財生産企業(経常) | 資本財生産企業(資本) |          銀行           | 統合政府(経常) | 統合政府(資本) | 合計 |
| :------------------- | :-----------: | :---------------------: | :------------------: | :------------------: | :------------------: | :------------------: | :---------------------: | :------------: | :------------: | ---- |
| 消費                 |    $-C_w$     |         $-C_i$          |         $+C$         |                      |                      |                      |                         |                |                | $0$  |
| 投資                 |               |                         |                      |        $-I_c$        |        $+I_f$        |        $-I_k$        |                         |                |                | $0$  |
| 公的固定資本形成     |               |                         |                      |                      |        $+I_g$        |                      |                         |                |     $-I_g$     | $0$  |
| 政府支出（賃金除く） |               |                         |         $+G$         |                      |                      |                      |                         |      $-G$      |                | $0$  |
| 賃金                 |     $+W$      |                         |        $-W_c$        |                      |        $-W_k$        |                      |         $-W_b$          |     $-W_g$     |                | $0$  |
| 所得税               |   $-T_{iw}$   |        $-T_{ii}$        |                      |                      |                      |                      |                         |     $+T_i$     |                | $0$  |
| 資産税               |   $-T_{ew}$   |        $-T_{ei}$        |      $-T_{ec}$       |                      |      $-T_{ek}$       |                      |                         |     $+T_e$     |                | $0$  |
| 法人税               |               |                         |      $-T_{fc}$       |                      |       $-T_{fk}$       |                      |        $-T_{fb}$        |     $+T_f$     |                | $0$  |
| 消費財生産企業利潤   |               |       $+\Pi_{ci}$       |       $-\Pi_c$       |     $+\Pi_{cc}$      |                      |                      |       $+\Pi_{cb}$       |                |                | $0$  |
| 資本財生産企業利潤   |               |       $+\Pi_{ki}$       |                      |                      |       $-\Pi_k$        |      $+\Pi_{kk}$      |       $+\Pi_{kb}$       |                |                | $0$  |
| 借入金金利           | $-r L_{w-1}$  |                         |     $-r L_{c-1}$     |                      |     $-r L_{k-1}$     |                      |       $+r L_{-1}$       |                |                | $0$  |
| 政府の貯蓄           |               |                         |                      |                      |                      |                      |                         |     $-GS$      |     $+GS$      | $0$  |
| [メモ:Net Lending]   |    $NL_w$     |         $NL_i$          |                      |        $NL_c$        |                      |        $NL_k$        |         $NL_b$          |                |     $NL_g$     | $0$  |
| 預金                 | $-\Delta M_w$ |      $-\Delta M_i$      |                      |    $-\Delta M_c$     |                      |    $-\Delta M_k$     |       $+\Delta M$       |                |                | $0$  |
| 借入                 | $+\Delta L_w$ |                         |                      |    $+\Delta L_c$     |                      |    $+\Delta L_k$     |       $-\Delta L$       |                |                | $0$  |
| 消費財生産企業株式   |               | $-p_{ec} \Delta e_{ci}$ |                      | $+p_{ec} \Delta e_c$ |                      |                      | $-p_{ec} \Delta e_{cb}$ |                |                | $0$  |
| 資本財生産企業株式   |               | $-p_{ek} \Delta e_{ki}$ |                      |                      |                      | $+p_{ek} \Delta e_k$ | $-p_{ek} \Delta e_{kb}$ |                |                | $0$  |
| 現金                 | $-\Delta H_w$ |                         |                      |                      |                      |                      |      $-\Delta H_b$      |                |  $+\Delta H$   | $0$  |
| 合計                 |      $0$      |           $0$           |         $0$          |         $0$          |         $0$          |         $0$          |           $0$           |      $0$       |      $0$       |      |

## 5.2. バランスシート表

|                    | 労働者  |  資本家   | 消費財生産企業 | 資本財生産企業 |   銀行    | 統合政府 | 合計 |
| :----------------- | :-----: | :-------: | :------------: | :------------: | :-------: | :------: | :--: |
| 資本               |         |           |     $+K_c$     |     $+K_k$     |           |  $+K_g$  | $+K$ |
| 預金               | $+M_w$  |  $+M_i$   |     $+M_c$     |     $+M_k$     |   $-M$    |          | $0$  |
| 消費財生産企業株式 |         | $+E_{ci}$ |     $-E_c$     |                | $+E_{cb}$ |          | $0$  |
| 資本財生産企業株式 |         | $+E_{ki}$ |                |     $-E_k$     | $+E_{kb}$ |          | $0$  |
| 借入               | $-L_w$  |           |     $-L_c$     |     $-L_k$     |   $+L$    |          | $0$  |
| 現金               | $+H_w$  |           |                |                |  $+H_b$   |   $-H$   | $0$  |
| 純資産             | $-NW_w$ |  $-NW_i$  |    $-NW_c$     |    $-NW_k$     |  $-NW_b$  | $-NW_g$  | $-K$ |
| 合計               |   $0$   |    $0$    |      $0$       |      $0$       |    $0$    |   $0$    | $0$  |

## 5.3. 完全統合表

|                                      |    労働者     |           資本家            |       消費財生産企業        |       資本財生産企業        |    銀行     |          統合政府           |            合計            |
| :----------------------------------- | :-----------: | :-------------------------: | :-------------------------: | :-------------------------: | :---------: | :-------------------------: | :------------------------: |
| 期首純資産                           |  $NW_{w-1}$   |         $NW_{i-1}$          |         $NW_{c-1}$          |         $NW_{k-1}$          | $NW_{b-1}$  |         $NW_{g-1}$          |          $K_{-1}$          |
| 資本のキャピタルゲイン               |               |                             | $+\Delta p_k \cdot k_{c-1}$ | $+\Delta p_k \cdot k_{k-1}$ |             | $+\Delta p_k \cdot k_{g-1}$ | $+\Delta p_k \cdot k_{-1}$ |
| 資本の増減                           |               |                             |   $+p_k \cdot \Delta k_c$   |   $+p_k \cdot \Delta k_k$   |             |   $+p_k \cdot \Delta k_g$   |   $+p_k \cdot \Delta k$    |
| 預金の増減                           | $+\Delta M_w$ |        $+\Delta M_i$        |        $+\Delta M_c$        |        $+\Delta M_k$        | $-\Delta M$ |                             |            $0$             |
| 消費財生産企業株式のキャピタルゲイン |               | $+\Delta p_{ec} \cdot e_{ci-1}$ | $-\Delta p_{ec} \cdot e_{c-1}$  | | $+\Delta p_{ec} \cdot e_{cb-1}$ |             |             $0$             |
| 資本財生産企業株式のキャピタルゲイン |               | $+\Delta p_{ek} \cdot e_{ki-1}$ |  | $-\Delta p_{ek} \cdot e_{k-1}$ | $+\Delta p_{ek} \cdot e_{kb-1}$ |             |             $0$             |
| 消費財生産企業株式の増減                           |               |   $+p_{ec} \cdot \Delta e_{ci}$   |    $+p_{ec} \cdot \Delta e_c$    | |  $+p_{ec} \cdot \Delta e_{cb}$   |             |             $0$             |
| 資本財生産企業株式の増減                           |               |   $+p_{ek} \cdot \Delta e_{ki}$   |     |  $+p_{ek} \cdot \Delta e_k$  |  $+p_{ek} \cdot \Delta e_{kb}$   |             |             $0$             |
| 借入の増減                           | $-\Delta L_w$ |                             |        $-\Delta L_c$        | $-\Delta L_k$        |         $+\Delta L$         |             |             $0$             |
| 現金の増減                           | $+\Delta H_w$ |                     |        |                             |        $+\Delta H_b$        | $-\Delta H$ |             $0$             |                            |
| 期末純資産                           |    $NW_w$     |           $NW_i$            |           $NW_c$            | $NW_k$            |           $NW_b$            |   $NW_g$    |             $K$             |                            |

## 5.4. フローの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $-C_w+W-T_{iw}-T_{ew}-r L_{w-1} = \Delta M_w - \Delta L_w + \Delta H_w$
- [x] $-C_i-T_{ii}-T_{ei}+P_i = \Delta M_i + p_e \Delta e_i$
- [x] $\Pi = C+I_f+I_{gf}+G-W_f-T_{ef}-T_{ff}-r L_{f-1}$
- [x] $-I_f + \Pi_f = \Delta M_f - \Delta L_f - \Delta E$
- [x] $-W_b - T_{fb} + \Pi_b + r L_{-1} = -\Delta M + \Delta L + p_e \Delta e_b + \Delta H_b$
- [x] $-G+I_{gg}-W_g+T_i+T_e+T_f-GS = 0$
- [ ] $-I_g+GS= \Delta H$
- [x] $C = C_w + C_i$
- [x] $I_g = I_{gf} + I_{gg}$
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

- [x] $k_f = (1 - \beta_1)k_{f-1} + i_f$
- [x] $k_g = (1 - \beta_1)k_{g-1} + i_g$
- [x] $\Delta k_f = k_f - k_{f-1}$
- [x] $\Delta k_g = k_g - k_{g-1}$
- [x] $K_f = K_{f-1} + \Delta K_f = p_{-1} k_{f-1} + \Delta p \cdot k_{f-1} + p \cdot \Delta k_f$
- [x] $K_g = K_{g-1} + \Delta K_g = p_{-1} k_{g-1} + \Delta p \cdot k_{g-1} + p \cdot \Delta k_g$
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
- [x] $NW_g = NW_{g-1} + \Delta K_g - \Delta H (= NW_{g-1} + NL_g)$

## 5.6. ストックの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $K = K_f + K_g$
- [x] $M = M_w + M_i + M_f$
- [ ] $E = E_i + E_b$
- [x] $e = e_i + e_b$
- [x] $L = L_w + L_f$
- [x] $H = H_w + H_f$
- [ ] $NW_w + NW_i + NW_f + NW_b + NW_g = K$
- [ ] $NW_w = M_w - L_w + H_w$
- [ ] $NW_i = M_i + E_i$
- [ ] $NW_f = K_f + M_f - E - L_f$
- [ ] $NW_b = -M + E_b + L + H_b$
- [ ] $NW_g = K_g-H$

# 6. パラメータの値

- $\alpha_1 = 0.9$
- $\alpha_2 = 0.1$
- $\alpha_3 = 0.1$
- $\alpha_4 = 0.02$
- $\alpha_5 = 0.5$
- $\alpha_6 = 0.1$
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
- $\zeta_4 = 1.0$
- $\theta_1 = 0.4$
- $\theta_2 = 0.1$
- $\iota_1 = 0.1$
- $\iota_2 = 1.0$
- $\iota_3 = 0.9$
- $\iota_4 = 0.5$
- $\iota_5 = 5.0$
- $\iota_6 = 0.5$
- $\iota_7 = 0.5$
- $\iota_8 = 0.5$
- $\iota_9 = 0.02$
- $\iota_{10} = 0.02$
- $\iota_{11} = 2.0$
- $\lambda_e = 0.5$
- $\mu_1 = 0.3$
- $\mu_2 = 0.1$
- $\mu_3 = 0.5$
- $\tau_1 = 0.2$
- $\tau_2 = 0.2$
- $u^T = 0.8$
- $G_0 = 1$
- $r = 0.01$

