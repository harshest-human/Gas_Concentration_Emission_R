# Manuscript 3 literature knowledge base

This directory stores compact, auditable project memory for Manuscript 3. Original PDFs remain in `D:\Literature`. Extracted full text is a local search cache and is not version-controlled.

## Working sequence

```text
source registry
  -> extracted local text
  -> verified paper note
  -> cross-paper synthesis
  -> claim audit
  -> manuscript text
```

## Directory map

- `PROJECT_STATE.md`: current scientific framing, active files, decisions and next task.
- `00_source_registry/`: citation keys, authoritative PDF paths, DOI and verification state.
- `01_scope/`: research question, study facts and analysis provenance.
- `02_source_cache/`: local extracted text. This directory is ignored by Git.
- `03_paper_evidence/`: one concise verified note per central paper.
- `04_synthesis/`: evidence organised by manuscript section or scientific theme.
- `05_claim_audit/`: manuscript claims linked to sources and verification status.
- `06_prompts/`: reusable questions for NotebookLM, Gemini or another discovery tool.
- `07_search_gaps/`: missing sources and unresolved claims.

## Evidence rules

1. Search verified paper notes before extracted text or PDFs.
2. Treat the original PDF as authoritative.
3. Verify numerical values, sampling dimensions, instrument details, statistical tests and limitations in the primary source.
4. Record the supporting page, table or figure for every important claim.
5. Mark unavailable information as `NOT REPORTED`.
6. Compare results from different barns conditionally.
7. Do not commit extracted full papers or copied copyrighted text.

## Manuscript location

The active complete drafting file is outside this knowledge base:

```text
../Manuscript_3_Mapping_high_resolution_concentration_Latex_draft/
  LaTex_drafts/20260813_Manuscript3_draft.tex
```

Draft prose should be edited there after its supporting claims are entered in the evidence maps and claim register.
