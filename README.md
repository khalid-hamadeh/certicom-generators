# The recycled numbers inside Certicom's curve generators

[![reproduce](https://github.com/khalid-hamadeh/certicom-generators/actions/workflows/reproduce.yml/badge.svg)](https://github.com/khalid-hamadeh/certicom-generators/actions/workflows/reproduce.yml)

Every standard elliptic curve comes with a published starting point, its generator *G*. For the curves of SEC 2 (2000) and the binary-field examples of ANSI X9.62 (1998), most of them from Certicom, no document says how the generators were chosen. Bitcoin's curve, secp256k1, is one of them.

Many curves have a cofactor *h*: the full group of points is *h* times larger than the subgroup in use. Dividing each published generator by *h* (or by 2 when *h* = 1), and keeping every answer, gives a structured "hidden" value on 25 of 69 standard curves. Seven were known before. On the other 18, hidden values are shared between sibling curves, reused from earlier curves, or too uneven to be random. For example, the generators of c2tnb191v1, v2 and v3 look unrelated, but divided by their cofactors (2, 4 and 6) all three give the same 47 hex digits; only the last digit differs.

## Run it

```
python3 reproduce.py
```

No dependencies; under a minute on a laptop. Curve parameters come from OpenSSL's `crypto/ec/ec_curve.c`, and every generator is checked to lie on its curve with order *n* before use. The script runs four blind tests over every curve (short hidden values, shared digit runs between curves, runs shared with another curve's published generator, and uneven digits), then checks each specific claim in the write-up. It exits non-zero if any claim fails. The badge above runs it on every push.

The last lines of a run:

```
CLAIMS IN THE NOTE
   PASS  sect113r1: G = 2P with x(P) = 3
   PASS  sect113r2: G = 2P with x(P) = 2
   ...
   PASS  c2tnb191v1 G = 2P, x(P) = 3.B.4
   PASS  c2tnb191v2 G = 4P, x(P) = 3.B.d
   PASS  c2tnb191v3 G = 6P, x(P) = 3.B.0
   ...
   PASS  secp256k1 (Bitcoin) G = 2P, x(P) = 0x3b78ce...9c63

21 of 21 claims pass.
```

The full output is in `expected-output.txt`.

## Files

- `reproduce.py`: the curves, the division, the blind tests and the claim checks
- `expected-output.txt`: a complete run
- `CITATION.cff`: citation metadata

## Write-up

The patterns curve by curve, a family tree of which curves share or reuse a value, Bitcoin's 2009 source, and what was known before: [The recycled numbers inside Certicom's curve generators](https://khalidhamadeh.com/writing/certicom-generators).

## Notes

What this establishes: the shared and reused hidden values are exact computations, and unrelated curves share at most 6 digits where these share 23 to 56. That one routine starting from stored values made them is an inference. Where the stored values came from, including the 38-digit core behind secp256k1's generator (noticed on two curves in 2015 and on all four secp…k1 curves by 2022), is still unknown.

None of this affects security. A generator's origin cannot weaken the discrete-log problem on its curve.

The patterns were found by an AI agent (Anthropic's Claude) running blind tests over the standards, and the main results were rebuilt independently by a second model (OpenAI's Codex). The checks above do not depend on trusting either. Corrections and earlier sources are welcome as issues and will be credited.

Khalid Hamadeh, October 2026. MIT licence.
