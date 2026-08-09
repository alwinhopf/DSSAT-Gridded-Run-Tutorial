# Bundled executable status

The three binaries in this directory are legacy development artifacts. The
pipeline does **not** auto-select them: it uses `DSSAT_EXE`, `DSSAT_BASE`, or a
sibling `../DSSAT48` installation. These binaries are not self-contained and
will fail unless the matching DSSAT support installation (`DSSATPRO`,
`MODEL.ERR`, `StandardData`, genotype files, and related data) is available.

Their original source commit, build flags, and licence provenance have not been
recorded. Do not publish them as release assets or copy them into another
project until they are replaced by documented DSSAT-CSM-OS builds or removed
after confirming they are unnecessary. The hashes below identify the current
artifacts for corruption checks only; they are not provenance:

| Platform | File | SHA-256 |
|---|---|---|
| Linux x86-64 | `linux/dscsm048` | `58864b80adde7e10d50065fff29bdfb6a48af7836ce23bf0a2cf83d653eb4d85` |
| macOS arm64 | `macos/dscsm048` | `d1df7f41e93f88288cfe817ad94b9a3af67a41a11f74bdbce93207035f3c75b9` |
| Windows x86-64 | `windows/DSCSM048.EXE` | `7b75b4df50f61733b04b70922c1c26dc048108cae9e934535b5af94ccd9cec60` |

For supported setup instructions, use the repository README's **Installing
DSSAT** section and point `DSSAT_EXE` at that verified installation.
