# Record & Thinking: Model behaviour on ReAct tool-use

**Date**: 2026-02-24  
**Context**: Gene Agent (ReAct + MCP tools). Same codebase, same prompts, same query; only the Ollama model was changed.

---

## Finding

**Under identical setup (code, skills, prompt, query), different Ollama models can succeed or fail differently on the same task.**

In particular, for the **gene_proteomics** skill and the query *"Give me a full analysis of TP53"*:

- **qwen3:30b-a3b**: Calls `get_uniprot_function` once with `organism: 'human'`, then on the next turn outputs the full "Protein Function & Disease Report" — **succeeds** (2 tool calls total for the skill).
- **qwen3-coder-next**: Keeps calling `get_uniprot_function` in a loop (alternating `organism: 'Homo sapiens'` and `organism: 'human'`), never produces a final text report in the same skill run — **fails** (hits max iterations, 6 tool calls, no final answer).

So: **instruction-following and tool-use discipline** (e.g. “call once, then output the report”) vary by model; some models complete the flow directly, others get stuck in repeated tool calls without returning a final answer.

---

## Summary table

| Model               | gene_genomics (TP53) | gene_proteomics (TP53) | Outcome                          |
|---------------------|----------------------|-------------------------|----------------------------------|
| qwen3-coder-next    | ✅ 1 call → report   | ❌ 5 rounds of tool calls, no report | Fails (max iterations)          |
| qwen3:30b-a3b       | ✅ 1 call → report   | ✅ 1 call → report       | Succeeds (2 tool calls total)    |

---

## Logs (qwen3-coder-next — failing run)

```
✅ Loaded 2 skills: ['gene_genomics', 'gene_proteomics']
🤖 Ollama model: qwen3-coder-next

============================================================
🧬 Query: Give me a full analysis of TP53
🤖 Model: qwen3-coder-next
============================================================
🔌 MCP connected · tools discovered: ['search_gene_ncbi', 'get_uniprot_function']
🗺️  Router selected: ['gene_genomics', 'gene_proteomics']

────────────────────────────────────────────────────────────
📖 Running skill: gene_genomics
────────────────────────────────────────────────────────────

[Round 1] Reasoning...

🔧 Call tool: search_gene_ncbi({'gene_name': 'TP53', 'organism': 'human'})
   ↳ Returned 702 characters.

[Round 2] Reasoning...

📋 [gene_genomics] Result:
# Genomic Identity Card: TP53
...

────────────────────────────────────────────────────────────
📖 Running skill: gene_proteomics
────────────────────────────────────────────────────────────

[Round 1] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'Homo sapiens'})
   ↳ Returned 716 characters.

[Round 2] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'human'})
   ↳ Returned 716 characters.

[Round 3] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'Homo sapiens'})
   ↳ Returned 716 characters.

[Round 4] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'human'})
   ↳ Returned 716 characters.

[Round 5] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'Homo sapiens'})
   ↳ Returned 716 characters.

✅ Done, 6 tool call(s)
```

---

## Logs (qwen3:30b-a3b — succeeding run)

```
✅ Loaded 2 skills: ['gene_genomics', 'gene_proteomics']
🤖 Ollama model: qwen3:30b-a3b

============================================================
🧬 Query: Give me a full analysis of TP53
🤖 Model: qwen3:30b-a3b
============================================================
🔌 MCP connected · tools discovered: ['search_gene_ncbi', 'get_uniprot_function']
🗺️  Router selected: ['gene_genomics', 'gene_proteomics']

────────────────────────────────────────────────────────────
📖 Running skill: gene_genomics
────────────────────────────────────────────────────────────

[Round 1] Reasoning...

🔧 Call tool: search_gene_ncbi({'gene_name': 'TP53', 'organism': 'human'})
   ↳ Returned 702 characters.

[Round 2] Reasoning...

📋 [gene_genomics] Result:
# Genomic Identity Card: TP53

- **Official Symbol**: TP53
- **Full Name**: tumor protein p53
- **NCBI Gene ID**: 7157
- **Organism**: human
- **Location**: 17p13.1
- **Aliases**: BCC7, BMFS5, LFS1, P53, TRP53

## Genomic Summary
TP53 encodes a tumor suppressor protein with transcriptional activation, DNA binding, and oligomerization domains...

────────────────────────────────────────────────────────────
📖 Running skill: gene_proteomics
────────────────────────────────────────────────────────────

[Round 1] Reasoning...

🔧 Call tool: get_uniprot_function({'gene_name': 'TP53', 'organism': 'human'})
   ↳ Returned 716 characters.

[Round 2] Reasoning...

📋 [gene_proteomics] Result:
# Protein Function & Disease Report: TP53

- **Protein Name**: Cellular tumor antigen p53
- **UniProt Accession**: P04637
- **Function**: Multifunctional transcription factor that induces cell cycle arrest, DNA repair, or apoptosis in response to cellular stress. Acts as a tumor suppressor by regulating key pathways that prevent malignant transformation.

## Disease Associations
N/A

## Research Insight
TP53 is the most frequently mutated tumor suppressor gene in human cancers, with over 50% of all malignancies harboring p53 pathway disruptions, making it a central focus for cancer therapeutics and biomarker development.

## Resource Link
- [UniProtKB](https://www.uniprot.org/uniprotkb/P04637)

✅ Done, 2 tool call(s)
```

---

## Takeaway

For ReAct-style tool use, **model choice matters**: the same pipeline can work with one model and loop/fail with another. Worth recording which models behave well under your prompts and tool contracts (e.g. “one call then report”) for future reference and for fallbacks (e.g. reuse of prior tool result when the same tool is called again).
