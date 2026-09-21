# Repository publication and manuscript-link update

Date: 2026-09-21.

The initial project snapshot was successfully pushed to
`git@github.com:maris205/riemann_clock.git` over SSH on branch `main`.
Initial commit: `573346bdf091266bdd2617d824176dcc98f9fb9a`.

Both manuscripts now identify the public project at
<https://github.com/maris205/riemann_clock> in their data/code availability
paragraphs. The two project READMEs provide the same URL and SSH clone command.
The recovery guide explains which large third-party downloads are excluded
from Git and how to restore and check them. The original local files are retained.

## Scoped checks for this revision

- Both papers compiled successfully: spectroscopy 20 pages; original 35 pages.
- PDF read-integrity preflight passed for both rebuilt files.
- Actual PDF URI annotations point to the intended repository. Upstream
  Riemann-data and Zenodo links remain intact.
- Scoped visual inspection covered spectroscopy pp. 15–16 and original pp. 32–33:
  repository URLs and their surrounding paragraphs render without clipping.
- The focused delivery audit passed 72/72 grouped checks. The original paper's
  saved-output delivery audit was also rerun; its separate report records the count.
- No scientific/model source, stored numerical result, figure or prediction
  protocol was changed. The earlier scientific review's main.tex hash identifies
  the pre-link-update snapshot; this revision changes only that file's data/code
  availability paragraph.
- An independent read-only inventory confirmed that both papers' recursively
  referenced TeX inputs, bibliographies and figures are tracked, including the
  prospective appendix figure. Source manifests and upstream licenses are retained.

The `spectroscopy_delivery_verification.json` and
`available_data_delivery_verification.json` reports identify the rebuilt PDFs and
current source snapshots. Earlier reports remain historical records of their
stated analysis stage. This publication check adds no new physical evidence or
external peer review.
