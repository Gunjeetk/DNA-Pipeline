# WCSA DNA Sequencing Pipeline  
**NRT Researcher-a-thon Submission**

A research-oriented DNA sequencing and candidate alignment prototype developed for the **NRT Researcher-a-thon**. This repository demonstrates a **traditional DNA sequencing workflow enhanced with a Weighted Column-Sum Alignment (WCSA) scoring engine**, designed as a **hardware-friendly alternative to conventional alignment scoring**.

---

## Problem
Traditional DNA sequencing pipelines spend significant compute time in the **candidate evaluation and alignment stage**, especially after basecalling and candidate window generation. As read volumes increase, this stage becomes a major bottleneck for **real-time genomic analysis**.

---

## Innovation
This project introduces **WCSA (Weighted Column-Sum Alignment)**, a novel scoring mechanism that reformulates alignment-event detection into a **column-wise weighted sum arithmetic operation**.

Instead of relying on conventional dynamic programming or direct string-by-string comparison, WCSA converts aligned bases into numeric digits and computes:

```text
S = ref × 1 + read × α
```

Where:
- `ref` = encoded reference base digit
- `read` = encoded read base digit
- `α` = row-weight multiplier (`α = 5`)

This enables:
- simplified match / mismatch decoding
- lightweight arithmetic scoring
- hardware-friendly accelerator design
- efficient candidate evaluation
- future FPGA / ASIC mapping opportunities

---

## Why It Matters
WCSA is designed as a **compute-efficient DNA alignment scoring prototype** with strong relevance to:

- nanopore-style sequencing pipelines
- hardware acceleration research
- near-memory genomic compute
- edge and handheld sequencers
- publishable architecture + systems research

This work is especially suitable for **algorithm–hardware co-design exploration**.

---

## Key Features
- End-to-end toy DNA sequencing workflow
- Synthetic nanopore-like signal generation
- Lightweight basecalling
- Candidate window generation
- Traditional scoring baseline
- **WCSA scoring engine**
- Runtime comparison metrics
- Explainability alignment trace
- Competition-ready reproducible notebook + Python script

---

## Repository Structure
```text
wcsa-dna-pipeline/
├── README.md
├── Requirements.txt
├── DNA_pipeline_wcsa_NRT_Hackhathon.ipynb
├── WCSA_vs_Traditional_DNA_pipeline.py
├── mapped_output.txt
└── pipeline_comparison_metrics.png
```

---

## Pipeline Overview
The current prototype follows this workflow:

```text
Raw Signal Acquisition
        ↓
Signal Transfer
        ↓
Basecalling
        ↓
Read Quality Control
        ↓
Candidate Window Generation
        ↓
Traditional Scoring
        ↓
WCSA Scoring Engine
        ↓
Best-Hit Selection
        ↓
Final Mapped Read Output
```

---

## WCSA Core Logic
The WCSA engine performs four internal stages:

1. **Base Digit Encoder**  
   `A,C,G,T,- → 1,2,3,4,0`

2. **Column-Sum Generator**  
   `S = ref + α × read`

3. **Alignment Event Decoder**  
   Detects:
   - match
   - mismatch
   - gap

4. **Window Score Aggregator**  
   Computes cumulative candidate score.

---

## Installation
Clone the repository:

```bash
git clone https://github.com/<your-username>/wcsa-dna-pipeline.git
cd wcsa-dna-pipeline
```

Install dependencies:

```bash
pip install -r Requirements.txt
```

---

## Quick Start
Run the complete WCSA vs Traditional comparison pipeline:

```bash
python WCSA_vs_Traditional_DNA_pipeline.py
```

---

## Example Outputs
This repository currently includes:

- `mapped_output.txt` → best-hit alignment output
- `pipeline_comparison_metrics.png` → runtime and scoring comparison
- notebook explainability traces
- WCSA column-by-column arithmetic validation

These outputs are intentionally stored in the **repository root for quick hackathon review**.

---

## Evaluation Goals
This repository evaluates WCSA against a traditional scoring baseline using:

- execution time
- scoring complexity
- candidate scoring behavior
- explainability
- hardware suitability
- scalability potential

---

## Competition Context
Developed for the **NRT Researcher-a-thon** as a **research prototype for compute-efficient DNA alignment scoring**.

The repository is structured for:
- quick judge onboarding
- easy execution
- reproducible notebook validation
- clean research communication
- future publication extension

---

## Future Work
- larger synthetic reference genomes
- minimizer-based candidate generation
- indel-aware robustness testing
- FPGA / ASIC WCSA mapping
- near-memory accelerator simulation
- handheld sequencer deployment feasibility

---

## Author
**Gunjeet Kaur**  
PhD Researcher – Electrical & Computer Engineering  
DNA sequencing acceleration • WCSA • architecture research

---



---

## License

MIT License

