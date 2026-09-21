## Material Passport

- Origin Skill: academic-research-suite / experiment-agent
- Origin Mode: validate (bounded pre-fit selection review)
- Origin Date: 2026-09-21
- Verification Status: VERIFIED for raw-source reconstruction and frozen-window counts; morphology assessment is conditional.
- Version Label: j0530_selection_review_v1

# J053007−250329: independent pre-fit selection review

The fixed four-line selection is acceptable for exploratory shared-gas profile fitting, with material saturation, continuum, and blend limitations. This review made no changes to masks, model configuration, or data and performed no profile fits.

## Source and numerical checks

The archived original FITS was opened directly from `data/raw/archive_expansion/J053007-250329_Final_Spectrum.tar.gz`. Its SHA256 agrees with the source manifest. The full coadd wavelength grid was reconstructed from FITS header values; normalized flux, expected fluctuations, and status exactly reproduce the processed full-coadd arrays. Atomic isotope-weighted reference wavelengths were separately parsed from the archived VPFIT atomic list. The coadd SHA256 agrees with `selection_frozen.json`.

There are exactly **120 valid native pixels in each of 1608, 1611, 2374 and 2382**, hence **480 total**, inside the same frozen **[−110, +130] km/s** window. All selected source indices are contiguous with raw status 1, finite flux, and positive finite expected fluctuation. Native spacing is approximately 1.999993 km/s. No residual-based rejection or extra mask was introduced. The corresponding source-index lists and checks are saved in `results/archive_complete/J053007-250329/selection_review.json`.

## Morphology and the weak 1611 line

I inspected the saved full ±500 km/s and enlarged pre-fit plots against the full coadd arrays. The isolated feature near −77 km/s, the principal complex from about −60 to +30 km/s, and the structure around +80 km/s recur in 1608, 2374 and 2382. Their common velocity extent supports a shared absorber model.

Fe II 1611 is useful here despite being weak. Its deepest flux is **0.84699 at +2.92 km/s**, matching the principal 1608 trough at **+2.96 km/s**. The mean flux over [−20,+20] km/s is approximately **0.916**, against nearby edge levels approximately **0.99**. Its archived total oscillator strength is approximately 0.00138, much lower than 1608 (0.05770), 2374 (0.03130), or 2382 (0.32000). It therefore supplies a weak, unsaturated constraint on the main column distribution that the black strong-line core cannot supply alone. This is not proof of contamination-free absorption or a precise independent velocity measurement. Faint −77 and +80 structures are not individually established in 1611.

## Material caveats

- **Saturation:** 2382 has 34 of 120 window pixels below normalized flux 0.05 (33 below three expected-fluctuation errors); the corresponding counts below 0.05 are 7 for 1608 and 6 for 2374. Multiple narrow gas decompositions can reproduce such cores. Fitting success must be assessed before interpreting relative shifts.
- **Continuum:** in [−110,−95] km/s, weighted flux means for 2374 and 2382 are approximately 0.964 and 0.965; their red-edge means over [+120,+130] are approximately 0.989 and 0.973. The 1611 blue/red edge means are approximately 0.995/0.988. Thus “continuum at both ends” should mean approximately continuum-like, with per-line continuum freedom; it must not imply known unity normalization.
- **Nearby separate features:** the pronounced additional 1608 red absorption starts beyond approximately +145 km/s, outside the +130 boundary. The [+130,+145] interval averages approximately 1.002 for 1608 and gives a modest guard against that feature. Other obvious absorption near +280 km/s in 2382 and below −300 km/s in 1611 is outside the frozen window. Those positions alone do not establish blend identities, and weak wings or unresolved interlopers remain possible.
- **Selection:** the window follows observed morphology. The declaration that it was frozen before all relative-shift fits is recorded in the supplied selection; a present-day file hash cannot independently prove chronology. Treat this as an exploratory selection fixed for the current comparison, not a blind preregistration.
- **Statistical meaning:** band means and their diagonal-error standard errors are descriptive checks. This review neither measures covariance/calibration uncertainties nor assigns detection significance to them. It provides no α or cosmic-time inference.

## Fallacy scan

11/11 statistical-fallacy categories were considered. Selection/survivorship, look-elsewhere, and researcher-choice concerns are applicable: this target follows an earlier archive screen, and the window and extra weak line were chosen from morphology. They are addressed here by recording the choice, retaining all source-valid pixels, and limiting interpretation. Simpson aggregation, ecological inference, collider adjustment, prevalence claims, regression-to-the-mean inference, causal attribution, and reverse-causality claims are not used in this bounded review. No confirmatory or population-general conclusion is drawn.

## Reproduction verdict

**38/38 numerical checks passed.** This verifies the selected inputs and bookkeeping, not the adequacy or uniqueness of any later physical fit. See the JSON companion for exact hashes, source pixels, per-line counts, and diagnostic bands.
