# PrimerForge Development Notes

## Purpose

PrimerForge is an intelligent primer design platform focused on:

- Species-specific PCR
- eDNA assays
- Conservation genetics
- Population-wide primer robustness

The goal is to go beyond Primer3 by automatically evaluating and ranking primer pairs using biological evidence.

---

# Architecture

Configuration
↓

Download

↓

Cache

↓

Sequence Manager

↓

Merge

↓

MAFFT Alignment

↓

SNP Discovery

↓

Conserved Regions

↓

Primer3

↓

Primer Pair Generation

↓

Scoring

↓

Validation

↓

Reports

---

# Design Principles

- Modular architecture
- One responsibility per module
- Automated testing
- Preserve backwards compatibility
- Never rewrite stable code unnecessarily
- Complete one sprint before starting another

---

# Scoring Strategy

Current scoring includes:

- Thermodynamics
- Conservation
- SNP robustness
- Product size
- GC balance
- Tm balance
- Multiplex compatibility
- Specificity (future)

---

# Coding Rules

- Every new feature must include tests.
- All existing tests must continue passing.
- No duplicated logic.
- Small, well-defined commits.
- Tag important milestones.

---

# Current Version

v0.9-alpha

24 / 24 tests passing

---

# Next Sprint

Species-specific specificity engine

Automatically:

- Download target taxa
- Download related taxa
- Build BLAST database
- Evaluate every primer
- Reject cross-amplifying primers
- Keep only species-specific primers

---

# Publication Goal

PrimerForge v1.0

Publication-quality software for:

- PCR primer design
- eDNA primer design
- Species identification
- Biodiversity monitoring

Target publication:

Bioinformatics

or

BMC Bioinformatics

or

Nucleic Acids Research (Web Server/Database if expanded)

---

Last Updated

Sprint 3
Version v0.9-alpha
