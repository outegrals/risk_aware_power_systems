# Masked-consensus revision checks

Run from the repository root:

```sh
python3 scripts/validate_masked_consensus.py
```

The script writes `output/privacy/masked_consensus_validation.json`. It is a reproducible mathematical validation model, not a production cryptographic implementation. Seeds and masks are intentionally public test data.

The protocol uses public positive integer bounds B, resolution Delta, and prime p > nB; privately exchanged pairwise masks; real-valued average consensus of masked residues; and certified rounding of the masked integer sum before reduction modulo p. The utility participates with a zero initial value. The private operator mask graph is separate from the public star used for averaging.

The three-operator examples use exact rational arithmetic. One includes modular wraparound. Checks with 10, 50, and 100 operators use binary64 and validate numerical correctness on the saved inputs; they do not certify accumulated roundoff or benchmark deployed performance. All iteration budgets satisfy the ideal-arithmetic spectral bound, verified by an integer inequality.

Small exhaustive examples compare the complete public masked-input distribution and the larger view of a utility colluding with one operator. They illustrate the theorem's honest-mask-connectivity condition and its failure when the honest graph disconnects. The manuscript supplies the general proof; finite examples are not a replacement for it.

The guarantee assumes a fixed honest-but-curious coalition, a fixed roster and graph, private authenticated mask channels, and fresh independent uniform masks. It permits the forecast total and uncapped ERS quantities to be disclosed. It does not cover malicious participants, dropouts, compromised endpoints, later meter/settlement releases, or the long manuscript's saturated comparisons. Exact target bids are for encoded forecasts, before any additional market bid rounding.

Production secure-channel establishment, secret cryptographic randomness, and their byte or timing overhead are not implemented here. The pre-existing MATLAB experiments have not been converted into secure implementations.

Revised review PDFs:

- `output/pdf/main_privacy_revision.pdf`
- `output/pdf/main_LONG_privacy_revision.pdf`

Original submitted PDFs and the earlier beta-revision PDFs are retained.
