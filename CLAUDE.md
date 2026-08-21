# CLAUDE.md

이 프로젝트의 목표·규칙·데이터 스펙을 요약한 파일. 세션 시작 시 이 문서를 먼저 참고할 것.

## 프로젝트 개요

- 목표: `binding_affinity`(화합물-단백질 결합 친화도)를 `molecular_weight`, `logp`, `protein_pi` 등으로 예측하는 회귀분석
- 대상 사용자: 신약 후보를 스크리닝하는 연구원 (실험 전 컴퓨터로 1차 거르는 용도)
- 성공 기준: 초보 연구원이 실험 전 스크리닝에 참고할 만한 수준
- 자세한 문제 정의: [`docs/problem_definition.md`](docs/problem_definition.md) (v0, 오전 작성 — 오후에 갱신될 수 있음)

## 데이터 스펙

- 원본: `data/raw/drug_discovery_virtual_screening.csv` — 2,000행 × 17컬럼
- 결측치: `logp`, `polar_surface_area`, `hydrophobicity` 각 3.0%. 나머지 0%
- **`active` 컬럼은 `binding_affinity`를 임곗값 7.0으로 이진화한 파생값 — 피처로 쓰면 데이터 누수. 반드시 제외**
- `compound_id`: 완전 유일한 식별자(2000/2000). 피처 아님
- `protein_id`: 400개 그룹, 식별자 아닌 그룹 변수(그룹당 평균 5건)
- `mw_ratio` ≈ `molecular_weight`/`protein_length`, `logp_pi_interaction` ≈ `logp`×`protein_pi` — 다른 컬럼의 파생값, 다중공선성(상관 0.8 이상) 주의
- 자세한 진단 근거: `docs/WORKFLOW.md`의 "STAGE ① 원본 탐색" 항목

## 작업 규칙

- `data/raw/`는 절대 수정하지 않는다. 모든 변환은 `data/processed/`에서.
- 컬럼 삭제/값 변경은 **실행 전에 방법(안)을 제시하고 확인받은 뒤** 진행한다.
- `random_state`는 전 과정에서 42로 고정한다.
- 한글 출력 시 인코딩 깨짐 방지를 위해 `PYTHONIOENCODING=utf-8 python ...`로 실행한다.
- 데이터 전처리·분석 절차는 [`docs/reference/SOP.md`](docs/reference/SOP.md)의 STAGE ①~⑤ 순서를 따른다.
- 진행 상황과 판단 근거는 매 작업 후 [`docs/WORKFLOW.md`](docs/WORKFLOW.md)에 누적 기록한다 (최상단 "현재 상태" 블록을 최신화).

## 폴더 구조 (현재)

```
./
├── CLAUDE.md              # 이 파일
├── README.md              # 포트폴리오용 최종 보고서 (영문)
├── requirements.txt
├── data/
│   ├── raw/                # 원본, 수정 금지
│   └── processed/          # 전처리 산출물 + train/test 분할
├── docs/
│   ├── reference/
│   │   └── SOP.md               # 범용 분석 절차 (참고용, 이 프로젝트 전용 아님)
│   ├── archive/
│   │   └── CLAUDE_CODE_WORKFLOW.md  # 예전 다른 프로젝트(자전거 대여) 기록 — 이번 프로젝트와 무관
│   ├── problem_definition.md    # 이 프로젝트의 목표·문제 정의
│   ├── WORKFLOW.md              # 이 프로젝트의 진행 기록(로그) + Research Decision Log
│   ├── CLAUDE_WORKFLOW.md       # Claude Code 활용 방식 정리 (영문, representative workflow)
│   ├── biological_context.md    # 바이오 R&D 맥락 설명 (영문)
│   └── next_experiment.md       # 제안된 후속 실험 (영문, 미실행)
├── output/
│   └── day6/
│       ├── figures/, models/    # 전체 산출물
│       └── final/               # 핵심 결과만 curation + FINAL_RESULTS.md + INDEX.md
└── scripts/
```
