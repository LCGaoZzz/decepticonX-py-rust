# Third-party backend notice

This repository contains an MIT-licensed orchestration and adapter layer. It
does **not** vendor, redistribute, or install the optional algorithm backends,
their reference data, or their native binaries.

The workflow can interoperate at runtime with independently obtained software:

| Backend | Upstream used by this project | License/status observed at the pinned revision |
|---|---|---|
| Reference builders | [LCGaoZzz/decepticon-fast](https://github.com/LCGaoZzz/decepticon-fast) `a98b66c7da5fad81ba6a1eaddb9c826daac31fad` | Code declares MIT; bundled BayesPrism-derived annotation tables are identified by its NOTICE as GPL-2 or GPL-3 |
| CIBERSORT / CIBERSORT-ABS | [LCGaoZzz/python-cibersort-rs](https://github.com/LCGaoZzz/python-cibersort-rs) `v0.1.1` / `cd957af06357d2636ba18af1135bf70b488dfbbb` | GPL-3.0-or-later |
| EPIC | `LCGaoZzz/epic-py-rust` `2e661cccb5c9d7c8327bd65253e3495a1d3c1991` | Private repository; LICR academic/non-commercial terms; no redistribution without permission |
| DeconRNASeq | [LCGaoZzz/deconrnaseq-py](https://github.com/LCGaoZzz/deconrnaseq-py) `a2ec1b9ebe341628e945a6d25c6068cad975f6b5` | GPL-2.0-only as declared by that repository |
| MuSiC | [LCGaoZzz/music-py](https://github.com/LCGaoZzz/music-py) `a4464c84f87fbeb61f0fa7032625bad21add52dd` | GPL-3.0-or-later |

The original [DECEPTICONx](https://github.com/Hao-Zou-lab/DECEPTICONx) and
[DECEPTICON](https://github.com/Hao-Zou-lab/DECEPTICON) projects declare MIT
licenses. Their scientific methods and the original method authors should be
cited in publications.

Upstream MIT notices retained for attribution:

- DECEPTICON: Copyright (c) 2023 Fulan Deng.
- DECEPTICONx: Copyright (c) 2025 Hao-Zou-lab.

Permission is hereby granted, free of charge, to any person obtaining a copy
of the software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished
to do so, subject to the following conditions: the above copyright notice and
this permission notice shall be included in all copies or substantial portions
of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.

These license summaries are provided to explain why backends are plugins and
are not legal advice. Users and distributors remain responsible for checking
the applicable upstream terms, especially for EPIC and for any combined
distribution.
