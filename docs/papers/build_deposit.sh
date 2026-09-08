#!/bin/bash
# Rebuild the ChemRxiv deposit PDF from its two sources.
#
# deposit.md is DERIVED - never edit it. The two sources are
# screening_decomposition_preprint.md and screening_decomposition_supplementary.md;
# editing the derived copy is how a corrected paper ships unchanged.
# The standalone '---' rules are dropped because pandoc renders them as page-wide lines
# that read as section breaks in the PDF.
set -euo pipefail
cd "$(dirname "$0")"
{ grep -v '^---$' screening_decomposition_preprint.md
  printf '\n\\newpage\n\n'
  cat screening_decomposition_supplementary.md
} > chemrxiv/deposit.md
pandoc chemrxiv/deposit.md -o chemrxiv/screening_decomposition_preprint.pdf \
  --pdf-engine=xelatex -H chemrxiv/header.tex \
  -V geometry:margin=2.5cm -V fontsize=10pt -V colorlinks=true \
  --resource-path=.:figures
echo "built chemrxiv/screening_decomposition_preprint.pdf"
