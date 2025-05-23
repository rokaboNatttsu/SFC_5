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

要修正まとめ

- 資本財生産企業と消費財生産企業の株式の新規発行及び自社株買いのアルゴリズムを見直す。長期的に、実質の生産量の推移と概ね比例するようなアルゴリズムを作りたい
- 振動が収束するように

# 2. モデルでやりたいことリスト

- 政府支出の増加率の変化の影響をシミュレーションしてまとめる
  - 技術や制度などの生産面の制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 消費需要の制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 企業の資金繰りの制約が経済規模を決めている場合、政府の支出増加傾向がどのように影響するか
  - 政府の支出増加傾向の分類
    - 基準シナリオ（ほかの制約が影響しないことを確認）
    - 政府支出一定シナリオ
    - 政府支出急増シナリオ
- 投資需要のネックが企業の内部資金になる用のパラメータ調整すると、カレツキアンモデルとは別のアプローチで、「利潤主導型成長」と解釈されてきた現象を再現できる



# 3. 行動方程式、上から計算

- $\zeta_2 = \zeta_{2-1} (1 + \zeta_4)$
- $\zeta_3 = \zeta_{3-1} (1 + \zeta_4)$
- $W_c = (1 - \epsilon_7) W_{c-1} + \epsilon_7 ((1 - \epsilon_6)\Omega_c + \epsilon_6 W_{c-1} \exp(\epsilon_2 ((C_{w-1}^D + C_{i-1}^D + G_{-1}^D)/(p_{-1} \min(\zeta_1 (k_{c-1}-Δk_{c-1}), \zeta_{2-1})) - u^T)))$
- $W_k = (1 - \epsilon_7) W_{k-1} + \epsilon_7 ((1 - \epsilon_6)\Omega_k + \epsilon_6 W_{k-1} \exp(\epsilon_2 ((i_{c-1}^D + i_{k-1}^D + i_{g-1}^D)/(min(\zeta_5 (k_{k-1}-Δk_{k-1}), \zeta_{3-1})) - u^T)))$
- $W_b = (1 - \epsilon_1)W_{b-1} + \epsilon_1 \cdot \max(0, \epsilon_3(\Pi_{cb-1} + \Pi_{kb-1} + r L_{-1}) + \epsilon_4 (H_{b-1} - M_{-1}))$
- $T_{ec} = \gamma_1 K_{c-1}$
- $T_{ek} = \gamma_1 K_{k-1}$
- $T_{ei} = \gamma_2 (M_{i-1} + E_{ci-1} + E_{ki-1})$
- $T_{ew} = \gamma_2 (M_{w-1} + H_{w-1})$
- $p = (1 + mu_c) \frac{p_k i_c^e +  W_c+T_{ec}+T_{fc}^e+r L_{c-1}}{c^e+g^e}$
  - $p = p_{-1} \exp\{\mu_3 \min(1, \max(u_{c-1}^e, \frac{C_{w-1}^D + C_{i-1}^D + G_{-1}^D}{p_{-1} \zeta_{2-1}}) - u^T)\}$
- $p_k = (1 + mu_k) \frac{W_k+T_{ek}+T_{fk}^e+r L_{k-1}}{i_c^e+i_g^e}$
  - $p_k = p_{k-1} \exp\{\mu_3 \min(1, \max(u_{k-1}^e, \frac{I_{f-1}^D + I_{g-1}^D}{p_{k-1} \zeta_{3-1}}) - u^T)\}$
- $w_g = w_{g-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $g^D = g^D_{-1}(1 + \delta_1 - \delta_2 \frac{p-p_{-1}}{p_{-1}})$
- $C_w^D = ((1 - \alpha_6) C_{w-1}^D + \alpha_6 \cdot \max\{0, \alpha_1 (W^e-T_{iw}^e-T_{ew}^e-r L_{w-1}) + \alpha_2 (M_{w-1} + H_{w-1} - L_{w-1})\}) \frac{C_{w-1}}{C_{w-1}-\Delta C_{w-1}}$
- $C_i^D = ((1 - \alpha_6) C_{i-1}^D + \alpha_6 \cdot \max(0, \alpha_3 (\Pi_{ci}^e+\Pi_{ki}^e-T_{ii}^e-T_{ei}^e) + (\alpha_4 + \alpha_5 \frac{E_{ci-1} + E_{ki-1}}{E_{c-1} + E_{k-1}}) M_{i-1})) \frac{C_{i-1}}{C_{i-1}-\Delta C_{i-1}}$
- $C_w = \frac{C_w^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{c-1}, \zeta_2)})}$
- $C_i = \frac{C_i^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{c-1}, \zeta_2)})}$
- $G = \frac{G^D}{\max(1, \frac{C_w^D + C_i^D + G^D}{p \cdot \min(\zeta_1 k_{c-1}, \zeta_2)})}$
- $i_c^D = btw(0, (u_{c-1} - u^T)k_{c-1} + \beta_1 k_{c-1}, max(0, \beta_2 \frac{M_{c-1} - L_{c-1}}{p_k}))$
- $i_k^D = btw(0, (u_{k-1} - u^T)k_{k-1} + \beta_1 k_{k-1}, max(0, \beta_2 \frac{M_{k-1} - L_{k-1}}{p_k}))$
- $i_g^D = i_{g-1}^D(1 + \delta_1 - \delta_2 \frac{p_k-p_{k-1}}{p_{k-1}})$
- $i_c = \frac{i_c^D}{\max(1, \frac{i_c^D + i_k^D + i_g^D}{\min(\zeta_5 k_{k-1}, \zeta_3)})}$
- $i_k = \frac{i_k^D}{\max(1, \frac{i_c^D + i_k^D + i_g^D}{\min(\zeta_5 k_{k-1}, \zeta_3)})}$
- $i_g = \frac{i_g^D}{\max(1, \frac{i_c^D + i_k^D + i_g^D}{\min(\zeta_5 k_{k-1}, \zeta_3)})}$
- $\Delta e_c = \max(0, \frac{\kappa_1 (p_k i_c^D-\beta_1 K_{c-1})}{p_{ec-1}}) - \max(0, \frac{(\kappa_2 + \max(0, \kappa_3 (u^T - u_{c-1}))) (M_{c-1} - L_{c-1})}{p_{ec-1}})$
- $\Delta e_k = \max(0, \frac{\kappa_1 (p_k i_k^D-\beta_1 K_{k-1})}{p_{ek-1}}) - \max(0, \frac{(\kappa_2 + \max(0, \kappa_3*(u^T - u_{k-1}))) (M_{k-1} - L_{k-1})}{p_{ek-1}})$
- $T_{ii} = \tau_1 (\Pi_{ci-1} + \Pi_{ki-1})$
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
- $E_{ci}^T= \kappa_4 NW_i^e \frac{C+G}{C+G+I_c+I_g}$
- $E_{ki}^T= \kappa_4 NW_i^e \frac{I_c+I_g}{C+G+I_c+I_g}$
- $E_{cb}^T = ((1 - \kappa_5) E_{cb-1} + \kappa_5 \frac{r_{E-1}}{r + r_{E-1}} (L_{-1} + E_{cb-1} + E_{kb-1}))\frac{E_{c-1}}{E_{c-1}+E_{k-1}}$
- $E_{kb}^T = ((1 - \kappa_5) E_{kb-1} + \kappa_5 \frac{r_{E-1}}{r + r_{E-1}} (L_{-1} + E_{cb-1} + E_{kb-1}))\frac{E_{k-1}}{E_{c-1}+E_{k-1}}$
- $p_{ec} = \frac{E_{ci}^T + E_{cb}^T}{e_c}$
- $p_{ek} = \frac{E_{ki}^T + E_{kb}^T}{e_k}$
- $E_{ci} = min(\kappa_4 NW_i^e \frac{C+G}{C+G+I_c+I_g}, E_c)$
- $E_{ki} = min(\kappa_4 NW_i^e \frac{I_c+I_g}{C+G+I_c+I_g}, E_k)$
- $M_i = NW_i - E_{ci} - E_{ki}$
- $L_c^D = (1 - \iota_4)L_{c-1}^D + \iota_4 max(0, \iota_3 (I_c + W_c + T_{ec} + T_{fc} + r L_{c-1}) - M_{c-1})$
- $L_c = btw(0, L_c^D, max(0, \iota_5 NL_c))$
- $L_k^D = (1 - \iota_4)L_{k-1}^D + \iota_4 max(0, \iota_3 (W_k+T_{ek}+T_{fk}+r L_{k-1}) - M_{k-1})$
- $L_k = btw(0, L_k^D, max(0, \iota_5 NL_k))$

# 4. 定義式

- $x^e = ((1 - \lambda_e) x^e_{-1} + \lambda_e x_{-1})$
- $u_c = \frac{c + g}{\zeta_1 k_{c-1}}$
- $u_k = \frac{i_f + i_g}{\zeta_5 k_{k-1}}$
- $NL_w = -C_w+W-T_{iw}-T_{ew}-r L_{w-1}$
- $NL_i = -C_i-T_{ii}-T_{ei}+\Pi_{ci} + \Pi_{ki}$
- $NL_c = -I_c + \Pi_c$
- $NL_k = -I_k + \Pi_k$
- $NL_b = -W_b - T_{fb} + \Pi_{cb} + \Pi_{kb} + r L_{-1}$
- $NL_g = -I_g + GS$
- $btw(A, B, C) = max(A, min(B, C))$ ただし( $A \leq C$ )
- $r_E = \frac{\Pi_{ci} + \Pi_{cb} + \Pi_{ki} + \Pi_{kb}}{E_c + E_k}$
- $Ωc = ϵ5*p_{-1} (c_{-1} + g_{-1})$
- $Ωk = ϵ5*p_{-1}*(i_{c-1} + i_{g-1})$

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
- $E_{ci} = p_{ce} e_{ci}$
- $E_c = p_{ce} e_c$
- $E_{cb} = p_{ce} e_{cb}$
- $E_{ki} = p_{ke} e_{ki}$
- $E_k = p_{ke} e_k$
- $E_{kb} = p_{ke} e_{kb}$

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
| 消費財生産企業株式の増減                           |               |   $+p_{ec} \cdot \Delta e_{ci}$   |    $-p_{ec} \cdot \Delta e_c$    | |  $+p_{ec} \cdot \Delta e_{cb}$   |             |             $0$             |
| 資本財生産企業株式の増減                           |               |   $+p_{ek} \cdot \Delta e_{ki}$   |     |  $-p_{ek} \cdot \Delta e_k$  |  $+p_{ek} \cdot \Delta e_{kb}$   |             |             $0$             |
| 借入の増減                           | $-\Delta L_w$ |                             |        $-\Delta L_c$        | $-\Delta L_k$        |         $+\Delta L$         |             |             $0$             |
| 現金の増減                           | $+\Delta H_w$ |                     |        |                             |        $+\Delta H_b$        | $-\Delta H$ |             $0$             |                            |
| 期末純資産                           |    $NW_w$     |           $NW_i$            |           $NW_c$            | $NW_k$            |           $NW_b$            |   $NW_g$    |             $K$             |                            |

## 5.4. フローの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $-C_w+W-T_{iw}-T_{ew}-r L_{w-1} = \Delta M_w - \Delta L_w + \Delta H_w$
- [ ] $-C_i-T_{ii}-T_{ei}+\Pi_{ci}+\Pi_{ki} = \Delta M_i + p_{ec} \Delta e_{ci} + p_{ek} \Delta e_{ki}$
- [x] $\Pi_c = C+G-W_c-T_{ec}-T_{fc}-r L_{c-1}$
- [x] $\Delta M_c = -I_c + \Pi_{cc} + \Delta L_c + p_{ec} \Delta e_c$
- [x] $\Pi_k = I_f + I_g-W_k-T_{ek}-T_{fk}-r L_{k-1}$
- [x] $\Delta M_k = -I_k + \Pi_{kk} + \Delta L_k + p_{ek} \Delta e_k$
- [x] $-W_b - T_{fb} + \Pi_{cb} + \Pi_{kb} + r L_{-1} = -\Delta M + \Delta L + p_{ec} \Delta e_{cb} + p_{ek} \Delta e_{kb} + \Delta H_b$
- [x] $-G-W_g+T_i+T_e+T_f-GS = 0$
- [ ] $-I_g+GS+\Delta H=0$
- [x] $C = C_w + C_i$
- [x] $I_f = I_c + I_k$
- [x] $W = W_f + W_g + W_b$
- [x] $T_i = T_{iw} + T_{ii}$
- [x] $T_e = T_{ew} + T_{ei} + T_{ec} + T_{ek}$
- [x] $T_f = T_{fc} + T_{fk} + T_{fb}$
- [x] $\Pi_{cf} = \Pi_c - \Pi_{ci} - \Pi_{cb}$
- [x] $\Pi_{kf} = \Pi_k - \Pi_{ki} - \Pi_{kb}$
- [ ] $NL_w + NL_i + NL_c + NL_k + NL_b + NL_g = 0$
- [x] $\Delta M = \Delta M_w + \Delta M_i + \Delta M_c + \Delta M_k$
- [x] $\Delta L = \Delta L_w + \Delta L_c + \Delta L_k$
- [ ] $\Delta e_{ci} + \Delta e_{cb} - \Delta e_c = 0$
- [ ] $\Delta e_{ki} + \Delta e_{kb} - \Delta e_k = 0$
- [x] $\Delta H = \Delta H_w + \Delta H_b$
- [x] $\Delta k = \Delta k_c + \Delta k_k + \Delta k_g$

## 5.5. ストックとフローの接続の整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $k_c = (1 - \beta_1)k_{c-1} + i_c$
- [x] $k_k = (1 - \beta_1)k_{k-1} + i_k$
- [x] $k_g = (1 - \beta_1)k_{g-1} + i_g$
- [x] $\Delta k_c = k_c - k_{c-1}$
- [x] $\Delta k_k = k_k - k_{k-1}$
- [x] $\Delta k_g = k_g - k_{g-1}$
- [x] $\Delta K_c = K_c - K_{c-1}$
- [x] $\Delta K_k = K_k - K_{k-1}$
- [x] $\Delta K_g = K_g - K_{g-1}$
- [ ] $K = K_{-1} + \Delta K$
- [ ] $k = k_{-1} + \Delta k$
- [ ] $\Delta K_c = \Delta p_k \cdot k_{c-1} + p_k \cdot \Delta k_c$
- [ ] $\Delta K_k = \Delta p_k \cdot k_{k-1} + p_k \cdot \Delta k_k$
- [ ] $\Delta K_g = \Delta p_k \cdot k_{g-1} + p_k \cdot \Delta k_g$
- [ ] $\Delta K = \Delta p_k \cdot k_{-1} + p_k \cdot \Delta k$
- [x] $\Delta M_w = M_w - M_{w-1}$
- [x] $\Delta M_i = M_i - M_{i-1}$
- [x] $M_c = M_{c-1} + \Delta M_c$
- [x] $M_k = M_{k-1} + \Delta M_k$
- [ ] $M = M_{-1} + \Delta M$
- [x] $\Delta L_w = L_w - L_{w-1}$
- [x] $\Delta L_c = L_c - L_{c-1}$
- [x] $\Delta L_k = L_k - L_{k-1}$
- [ ] $L = L_{-1} + \Delta L$
- [x] $\Delta E_{ci} = E_{ci} - E_{ci-1}$
- [x] $\Delta E_{ki} = E_{ki} - E_{ki-1}$
- [x] $\Delta E_{cb} = E_{cb} - E_{cb-1}$
- [x] $\Delta E_{kb} = E_{kb} - E_{kb-1}$
- [x] $\Delta E_c = E_c - E_{c-1}$
- [x] $\Delta E_k = E_k - E_{k-1}$
- [ ] $\Delta E_{ci} = \Delta p_{ec} \cdot e_{ci-1} + p_{ec} \cdot \Delta e_{ci}$
- [ ] $\Delta E_{ki} = \Delta p_{ek} \cdot e_{ki-1} + p_{ek} \cdot \Delta e_{ki}$
- [ ] $\Delta E_{cb} = \Delta p_{ec} \cdot e_{cb-1} + p_{ec} \cdot \Delta e_{cb}$
- [ ] $\Delta E_{kb} = \Delta p_{ek} \cdot e_{kb-1} + p_{ek} \cdot \Delta e_{kb}$
- [ ] $\Delta E_c = \Delta p_{ec} \cdot e_{c-1} + p_{ec} \cdot \Delta e_c$
- [ ] $\Delta E_k = \Delta p_{ek} \cdot e_{k-1} + p_{ek} \cdot \Delta e_k$
- [x] $e_c = e_{c-1} + \Delta e_c$
- [x] $\Delta e_{ci} = e_{ci} - e_{ci-1}$
- [x] $\Delta e_{cb} = e_{cb} - e_{cb-1}$
- [x] $e_k = e_{k-1} + \Delta e_k$
- [x] $\Delta e_{ki} = e_{ki} - e_{ki-1}$
- [x] $\Delta e_{kb} = e_{kb} - e_{ck-1}$
- [x] $\Delta H_w = H_w - H_{w-1}$
- [ ] $H_b = H_{b-1} + \Delta H_b$
- [ ] $H = H_{-1} + \Delta H$
- [ ] $NW_w = NW_{w-1} + \Delta M_w - \Delta L_w + \Delta H_w$
  - [x] $NW_w = NW_{w-1} + NL_w$
- [ ] $NW_i = NW_{i-1} + \Delta M_i + \Delta p_{ec} \cdot e_{ci-1} + \Delta p_{ek} \cdot e_{ki-1} + p_{ec} \cdot \Delta e_{ci} + p_{ek} \cdot \Delta e_{ki}$
  - [x] $NW_i = NW_{i-1} + NL_i + \Delta p_{ec} \cdot e_{ci-1} + \Delta p_{ek} \cdot e_{ki-1}$
- [ ] $NW_c = NW_{c-1} + p_k \cdot \Delta k_c + \Delta p_k \cdot k_{c-1} + \Delta M_c - \Delta p_{ec} \cdot e_{c-1} - p_{ec} \cdot \Delta e_c - \Delta L_c$
  - [x] $NW_c = NW_{c-1} + NL_c - \Delta p_{ec} \cdot e_{c-1} + \Delta K_c$
- [ ] $NW_k = NW_{k-1} + p_k \cdot \Delta k_k + \Delta p_k \cdot k_{k-1} + \Delta M_k - \Delta p_{ek} \cdot e_{k-1} - p_{ek} \cdot \Delta e_k - \Delta L_k$
  - [x] $NW_k = NW_{k-1} + NL_k - \Delta p_{ek} \cdot e_{k-1} + \Delta K_k$
- [ ] $NW_b = NW_{b-1} - \Delta M + p_{ec} \cdot \Delta e_{cb} + \Delta p_{ec} \cdot e_{cb-1} + p_{ek} \cdot \Delta e_{kb} + \Delta p_{ek} \cdot e_{kb-1} + \Delta L + \Delta H_b$
  - [x] $NW_b = NW_{b-1} + NL_b + \Delta p_{ec} \cdot e_{cb-1} + \Delta p_{ek} \cdot e_{kb-1}$
- [ ] $NW_g = NW_{g-1} + p_k \cdot \Delta k_g + \Delta p_k \cdot k_{g-1} - \Delta H$
  - [x] $NW_g = NW_{g-1} + NL_g + \Delta K_g$

## 5.6. ストックの整合性

モデルで使う恒等式にチェックを入れ、隠れた恒等式にはチェックを入れない

- [x] $K = K_c + K_k + K_g$
- [x] $k = k_c + k_k + k_g$
- [x] $M = M_w + M_i + M_c + M_k$
- [x] $E_{cb} = E_c - E_{ci}$
- [x] $E_{kb} = E_k - E_{ki}$
- [ ] $e_c = e_{ci} + e_{cb}$
- [ ] $e_k = e_{ki} + e_{kb}$
- [x] $L = L_w + L_c + L_k$
- [x] $H = H_w + H_b$
- [ ] $NW_w + NW_i + NW_c + NW_k + NW_b + NW_g = K$
- [ ] $NW_w = M_w - L_w + H_w$
- [x] $NW_i = M_i + E_{ci} + E_{ki}$
- [ ] $NW_c = K_c + M_c - E_c - L_c$
- [ ] $NW_k = K_k + M_k - E_k - L_k$
- [x] $NW_b = -M + E_{cb} + E_{kb} + L + H_b$
- [ ] $NW_g = K_g-H$

# 6. パラメータの値

- $\alpha_1 = 0.9$
- $\alpha_2 = 0.1$
- $\alpha_3 = 0.1$
- $\alpha_4 = 0.02$
- $\alpha_5 = 0.1$
- $\alpha_6 = 0.5$
- $\beta_1 = 0.05$
- $\beta_2 = 1.0$
- $\gamma_1 = 0.015$
- $\gamma_2 = 0.02$
- $\delta_1 = 0.02$
- $\delta_2 = 0.3$
- $\epsilon_1 = 0.5$
- $\epsilon_2 = 0.5$
- $\epsilon_3 = 0.7$
- $\epsilon_4 = 0.05$
- $\epsilon_5 = 0.695$
- $\epsilon_6 = 0.7$
- $\epsilon_7 = 0.3$
- $\zeta_1 = 1.0$
- $\zeta_2 = 150.0$ (バーンイン期間の初期値)
- $\zeta_3 = 150.0$ (バーンイン期間の初期値)
- $\zeta_4 = 0.02$
- $\zeta_5 = 1.0$
- $\theta_1 = 0.4$
- $\theta_2 = 0.1$
- $\iota_1 = 0.1$
- $\iota_2 = 1.0$
- $\iota_3 = 2.0$
- $\iota_4 = 0.5$
- $\iota_5 = 5.0$
- $\kappa_1 = 0.5$
- $\kappa_2 = 0.01$
- $\kappa_3 = 0.1$
- $\kappa_4 = 0.9$
- $\kappa_5 = 0.5$
- $\lambda_e = 0.5$
- $\mu_1 = 0.3$
- $\mu_2 = 0.1$
- $\mu_3 = 0.5$
- $\tau_1 = 0.2$
- $\tau_2 = 0.2$
- $u^T = 0.8$
- $G_0 = 1$
- $r = 0.01$
- $mu_c = 0.3$
- $mu_k = 0.5$