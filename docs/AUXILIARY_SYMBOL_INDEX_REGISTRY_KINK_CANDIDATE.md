# Auxiliary Symbol Index — spatial registry kink candidate

Status: **CANDIDATE auxiliary index.** These symbols belong to `CANDIDATE_SPATIAL_REGISTRY_KINK_FEEDBACK.md` and are not promoted to the active normal-only theory.

### $j$
- English: repeat index along the registry/slip direction
- 한국어: registry/slip 방향 반복단위 index
- Mathematical meaning: labels neighboring repeats in the moving reference row
- Status: CANDIDATE DISCRETE INDEX

### $s_j$
- English: local registry displacement of repeat $j$
- 한국어: $j$번째 반복단위의 국소 registry 변위
- Mathematical definition: the moving-row position is $x_j=jb+s_j$
- Unit: length
- Status: CANDIDATE STRUCTURAL STATE

### $z_j$
$$
s_j=s_0+z_j b+\tilde s_j,\qquad z_j\in\mathbb Z.
$$
- English: local unwrapped registry-well index
- 한국어: 국소 unwrapped registry well index
- Physical meaning: integer count of completed lattice-period registry shifts at repeat $j$
- Unit: dimensionless
- Status: CANDIDATE TOPOLOGICAL LABEL

### $\tilde s_j$
- English: intrawell registry offset
- 한국어: well 내부 registry 편차
- Mathematical meaning: residual part of $s_j-s_0$ after subtracting $z_jb$
- Unit: length
- Status: CANDIDATE STATE

### $q_j^{\mathrm{k}}$
$$
q_j^{\mathrm{k}}=z_{j+1}-z_j.
$$
- English: discrete kink/topological charge
- 한국어: 이산 kink/위상 전하
- Physical meaning: detects a boundary between adjacent regions occupying different registry wells
- Unit: dimensionless
- Status: CANDIDATE DIAGNOSTIC

### $E_{\parallel}[\{s_j\}]$
$$
E_{\parallel}
=\sum_{j<k}
\left[
v_{m,n}\left(\left|(k-j)b+s_k-s_j\right|\right)
-v_{m,n}\left((k-j)b\right)
\right].
$$
- English: intrarow compatibility energy
- 한국어: 동일 row 내부 변위-호환성 에너지
- Physical meaning: exact same-row pair-energy change caused by nonuniform registry displacement of neighboring repeats
- Unit: energy
- Status: CANDIDATE EXACT ENERGY within the declared moving-row embedding

### $E_{\mathrm{rk}}[\{a_j,s_j\}]$
$$
E_{\mathrm{rk}}
=\sum_j U_0(a_j,s_j)+E_{\parallel}[\{s_j\}].
$$
- English: spatial registry-kink candidate energy
- 한국어: 공간 registry-kink 후보 에너지
- Physical meaning: cross-row periodic registry energy plus same-row compatibility energy
- Unit: energy
- Status: CANDIDATE TOTAL ENERGY; external work is separate

### $Q_{a,j}^{\mathrm{int}}$
$$
Q_{a,j}^{\mathrm{int}}=-\partial_aU_0(a_j,s_j).
$$
- English: local intrinsic normal generalized force
- 한국어: 국소 고유 normal 일반화힘
- Physical meaning: normal force on spacing coordinate $a_j$ produced by the local registry state
- Unit: force
- Status: CANDIDATE EXACT DERIVATIVE of $U_0$

### $a_{\mathrm{eq}}(s;Q_a)$
$$
\partial_aU_0(a_{\mathrm{eq}},s)=Q_a,
\qquad
\partial_a^2U_0(a_{\mathrm{eq}},s)>0.
$$
- English: stable normal equilibrium conditional on registry
- 한국어: registry 조건부 안정 normal 평형거리
- Physical meaning: normal spacing preferred by a fixed registry position and applied normal generalized force
- Unit: length
- Status: CANDIDATE DIAGNOSTIC

### $\Delta G_{\mathrm{kp}}$
- English: kink-pair nucleation activation barrier
- 한국어: kink-pair 핵생성 활성화 장벽
- Mathematical meaning: minimum energy-saddle excess for creating a spatially nonuniform registry transition path
- Unit: energy
- Status: OPEN CANDIDATE QUANTITY; must be computed from $E_{\mathrm{rk}}$, not replaced by an arbitrary $N\Delta G_s$
