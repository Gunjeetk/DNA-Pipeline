# DNA-Pipeline
WCSA DNA Sequencing Pipeline
# WCSA DNA Sequencing Pipeline

A research-oriented DNA sequencing and candidate alignment pipeline developed for the **NRT Researcher-a-thon**. This repository demonstrates an enhanced sequencing workflow with a **Weighted Column-Sum Alignment (WCSA) scoring engine** integrated into the mapping stage to explore a **hardware-friendly alternative to conventional alignment scoring**.

---

## Problem

Traditional DNA sequencing analysis pipelines spend significant compute time in the **candidate evaluation and alignment stage**, especially after basecalling and candidate window generation. As read volumes grow, this stage becomes a major bottleneck for efficient real-time genomic analysis.

---

## Innovation

This project introduces **WCSA (Weighted Column-Sum Alignment)**, a novel scoring mechanism that reformulates alignment-event detection into a **column-wise weighted sum operation**.

Instead of relying on conventional string-by-string comparison or dynamic programming at every step, WCSA converts aligned bases into numeric digits and computes:

```text
S = ref × 1 + read × α
```

where:

* `ref` = encoded reference base digit
* `read` = encoded read base digit
* `α` = row-weight multiplier (default 5)

This enables:

* simplified match/mismatch/gap decoding
* hardware-friendly arithmetic blocks
* efficient candidate scoring
* future accelerator integration opportunities

---

## Why It Matters

WCSA is designed as a **research prototype for compute-efficient DNA alignment scoring**, with strong relevance to:

* nanopore-style streaming pipelines
* hardware acceleration research
* near-memory genomic compute architectures
* low-power portable sequencers
* future handheld and edge sequencing devices

The approach is particularly suitable for **publishable architecture and systems research**, where algorithm–hardware co-design is central.

---

## Key Features

* End-to-end toy DNA sequencing pipeline
* Synthetic nanopore-style signal generation
* Lightweight basecalling stage
* Candidate window generation
* **WCSA scoring engine**
* Best-hit selection
* Baseline vs WCSA comparison
* Timing visualization and result plots
* Competition-ready reproducible workflow

---

## Repository Structure

wcsa-dna-pipeline/
├── README.md
├── Requirements.txt
├── DNA_pipeline_wcsa_NRT_Hackhathon.ipynb
├── WCSA_vs_Traditional_DNA_pipeline.py
├── mapped_output.txt
├── pipeline_comparison_metrics.png
```

---

## Pipeline Overview

The current prototype follows this system flow:

```text
Raw Signal Acquisition
        ↓
Signal Transfer
        ↓
Basecalling
        ↓
Read Quality Control
        ↓
Reference Indexing
        ↓
Candidate Window Generation
        ↓
WCSA Scoring Engine
        ↓
Best-Hit Selection
        ↓
Final Mapped Read Output
```

---

## WCSA Core Scoring Logic

The WCSA engine performs four internal stages:

1. **Base Digit Encoder**
   `A,C,G,T,- → 1,2,3,4,0`

2. **Column-Sum Generator**
   `S = ref + α × read`

3. **Alignment Event Decoder**
   Detects:

   * match
   * mismatch
   * gap

4. **Window Score Aggregator**
   Computes cumulative alignment score across the candidate window.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/<your-username>/wcsa-dna-pipeline.git
cd wcsa-dna-pipeline
```

Install dependencies:

```bash
pip install -r requirements.txt
```

---

## Quick Start

Run the full WCSA pipeline:

```bash
python scripts/run_pipeline.py
```

Run baseline comparison:

```bash
python scripts/compare_wcsa_vs_baseline.py
```

---

## Example Output

Generated outputs include:

* mapped read summaries
* candidate alignment scores
* stage-wise timing tables
* runtime comparison plots
* WCSA score distributions
* pipeline architecture figures

All generated files are saved under:

```text
results/
```

---

## Evaluation Goals

This repository is designed to evaluate WCSA against a traditional candidate scoring pipeline using:

* execution time
* scoring complexity
* candidate evaluation behavior
* hardware suitability
* scalability potential

Future comparisons may extend toward:

* GENPIP-inspired pipelines
* minimap2-style seed-and-extend baselines
* accelerator-aware genomic workflows

---

## Competition Context

Developed for the **NRT Researcher-a-thon** as a research prototype demonstrating a novel **DNA alignment scoring architecture**.

The repository is structured for:

* quick judge onboarding
* reproducible execution
* clean research communication
* future publication extension

---

## Future Work

* larger synthetic reference genomes
* minimizer-based candidate generation
* streaming chunk-level processing
* ONT-style signal chunk integration
* FPGA / ASIC scoring block mapping
* near-memory accelerator simulation
* handheld device deployment feasibility study

---

## Author

**Gunjeet Kaur 
**
PhD Researcher – Electrical & Computer Engineering
DNA sequencing acceleration • WCSA • architecture research

---

## License

MIT License

