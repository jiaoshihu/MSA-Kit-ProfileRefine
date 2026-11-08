# Profile-Based Star Refinement (experimental V4)

Based on MSA-Kit V2. After adaptive center-star alignment and insertion-slot adjustment, `profile_iterative_refine` performs bounded deterministic profile-profile DP on three bipartitions, for at most two iterations. Only strictly improving cross-group sum-of-pairs moves are accepted, preserving all ungapped input sequences and original row order. Alignment scoring: identical +2, mismatch -1, residue-gap -2, N=0, T/U equivalent. Limits: 3-64 sequences, width <=1200, <=700,000 DP cells per split. These limits can cause the procedure to skip refinement.

Build with `bash benchmark/build.sh`, test with `python3 benchmark/test_adapter.py`, and run with `bash benchmark/run.sh examples/tiny.fasta output.fasta 1 dna`. Consult `PROFILE_REFINEMENT.zh-CN.md` for details. This is experimental and not a complete reimplementation of MAFFT. Better internal SP does not guarantee improvements in benchmark Q/TC.
