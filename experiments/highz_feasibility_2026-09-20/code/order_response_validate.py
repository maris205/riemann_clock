#!/usr/bin/env python3
"""Independent numerical audit of the centered order-response sensitivity family.

This auditor owns only validation artifacts.  The continuous reference kernel is
implemented independently of the fitting code.  Statistical checks audit the
conditional calculation, not an empirical measurement of an instrument LSF.
"""
from pathlib import Path
import hashlib
import json
import time
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import CubicSpline
from scipy.signal import fftconvolve
from scipy.optimize import least_squares
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/order_response"
CHECKS = []
METRICS = {}
FWHM_FACTOR = 2.354820045


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def check(name, value, detail=None):
    CHECKS.append(dict(name=name, pass_=bool(value), detail=detail))


def close(name, observed, expected, atol=1e-10, rtol=1e-8):
    x, y = np.asarray(observed), np.asarray(expected)
    error = float(np.max(np.abs(x-y))) if x.size else 0.0
    check(name, np.allclose(x, y, atol=atol, rtol=rtol),
          dict(maximum_absolute_error=error, atol=atol, rtol=rtol))


def reference_kernel(velocity, fwhm, asymmetry):
    """Continuous unit-area PDF, in inverse km/s; positive a has positive skew."""
    d = 1.2*np.cbrt(asymmetry)
    sigma = np.sqrt((fwhm/FWHM_FACTOR)**2 - .16*d*d)
    return (.8*norm.pdf(velocity, loc=-.2*d, scale=sigma)
            + .2*norm.pdf(velocity, loc=.8*d, scale=sigma))


def kernel_audit():
    velocity = np.linspace(-15., 15., 60001)
    moments = []
    for width in [1.6, 2.05, 2.1, 2.6]:
        for asym in [-1., -.3, -1e-6, 0., 1e-6, .3, 1.]:
            kernel = reference_kernel(velocity, width, asym)
            observed = np.array([simpson(kernel*velocity**j, x=velocity)
                                 for j in range(5)])
            variance = (width/FWHM_FACTOR)**2
            d = 1.2*np.cbrt(asym)
            expected = np.array([1., 0., variance, .096*d**3,
                                 3*variance**2 + .0064*d**4])
            name = f"kernel width {width} asymmetry {asym}"
            check(name + " positive", np.all(kernel >= 0))
            close(name + " independent integrated moments 0..4", observed,
                  expected, atol=2e-12, rtol=2e-12)
            close(name + " reflection symmetry", kernel,
                  reference_kernel(-velocity, width, -asym), atol=1e-14)
            moments.append(dict(fwhm_km_s=width, asymmetry=asym,
                                observed_moments=observed.tolist(),
                                expected_moments=expected.tolist()))
    width = 2.1
    sigma = width/FWHM_FACTOR
    gaussian = norm.pdf(velocity, scale=sigma)
    close("a zero is exactly the equivalent Gaussian",
          reference_kernel(velocity, width, 0), gaussian, atol=1e-16, rtol=1e-15)
    # K'(a=0) = -(0.165888/6) G''' where the prime on K denotes a.
    expected = .027648*(velocity**3/sigma**6 - 3*velocity/sigma**4)*gaussian
    derivative = []
    for step in [1e-3, 1e-4, 1e-5, 1e-6]:
        observed = (reference_kernel(velocity, width, step)
                    - reference_kernel(velocity, width, -step))/(2*step)
        relative = float(np.linalg.norm(observed-expected)/np.linalg.norm(expected))
        derivative.append(dict(step=step, relative_l2_error=relative))
    check("centered asymmetry derivative converges at a zero",
          derivative[-1]["relative_l2_error"] < 1e-4,
          derivative)
    METRICS.update(kernel_moments=moments, kernel_derivative=derivative,
                   minimum_fwhm_for_strict_component_variance_positive_km_s=
                   FWHM_FACTOR*.48,
                   analytical_third_moment_per_a_km3_s3=.165888,
                   smoothness="Density is C1 at a=0; the |a|^(4/3) fourth-cumulant term generally prevents C2 smoothness.")


def independent_spline(model, width, asymmetry):
    """Direct sampled-kernel convolution, separate from Gaussian-filter mixtures."""
    step = model.template.step
    grid = model.template.grid
    velocity = np.arange(-15., 15.+step/2, step)
    kernel = reference_kernel(velocity, width, asymmetry)
    profile = 1.-fftconvolve(1.-model.template.intrinsic[model.line],
                            kernel, mode="same")*step
    return CubicSpline(grid, profile, extrapolate=False)


def independent_replay(model, parameters):
    values = {**model.frozen, **dict(zip(model.labels, parameters))}
    splines = {order: independent_spline(model, values[f"width_{order}"],
                                       values.get(f"asymmetry_{order}", 0.))
               for order in model.orders}
    residuals, details = [], []
    for block in model.blocks:
        order = block["order"]
        shift = values[f'velocity_{block["exposure"]}']
        if order == min(model.orders):
            shift += values.get("order_offset", 0.)
        profile = splines[order](block["coordinates"]-shift)@model.template.weights/2
        design = np.column_stack([profile, block["x"]*profile,
                                  np.ones(len(profile))])
        weighted = design/block["error"][:, None]
        beta = np.linalg.lstsq(weighted, block["flux"]/block["error"], rcond=1e-12)[0]
        prediction = design@beta
        residual = (prediction-block["flux"])/block["error"]
        residuals.append(residual)
        details.append(dict(beta=beta, model=prediction, residual=residual))
    return np.concatenate(residuals), details


def production_audit():
    import order_response_model as production
    manifest = json.loads((OUT/"model_run_manifest.json").read_text())
    for key, path in [("code_sha256", production.__file__),
                      ("protocol_sha256", OUT/"model_protocol.json"),
                      ("source_metadata_sha256", production.src.PRO/"metadata.json"),
                      ("source_template_sha256", production.src.TEMPLATE_PATH),
                      ("source_previous_fits_sha256", production.src.OUT/"fit_rows.json")]:
        check("manifest "+key, manifest[key] == sha(path))
    check("protocol matches reviewer pre-output snapshot",
          manifest["protocol_sha256"] == "1ad60e95cc5c9c03fd226930be3c7adce9636f45b851f301c80679fb6c37b32a")
    # Exact native data selection is inherited from the earlier raw-FITS audit;
    # explicitly bind this follow-up to that audited metadata and fit artifact.
    previous = json.loads((production.src.OUT/"validation.json").read_text())
    check("earlier raw-FITS audit passed", previous["pass_"])
    for path in ["data/processed/exposures/metadata.json", "results/exposures/fit_rows.json"]:
        check("earlier raw-FITS audit binding "+path, previous["hashes"][path] == sha(ROOT/path))
    model_cache = {}
    def make(line, kind, indices, frozen):
        key = (line, kind, tuple(indices), json.dumps(frozen, sort_keys=True))
        if key not in model_cache:
            model_cache[key] = production.JointModel(line, kind, indices, frozen)
        return model_cache[key]
    rows = []
    files = sorted(OUT.glob("*.json"))
    for path in files:
        row = json.loads(path.read_text())
        if (isinstance(row, dict) and "optimizer_success" in row
                and row["name"].startswith(("full_", "train2018_", "heldout_"))):
            rows.append((path, row))
    check("all 40 planned joint fits present", len(rows) == 40)
    residual_errors, objective_errors, active, gradients = [], [], [], []
    file_hashes = {}
    for path, fit in rows:
        name = fit["name"]
        model = make(fit["line"], fit["kind"], fit["exposure_indices"], fit["frozen"])
        check(name+" parameter labels", model.labels == fit["labels"])
        p = np.array([fit["parameters"][label] for label in model.labels])
        _, low, high = model.initial()
        check(name+" parameters finite in bounds", np.isfinite(p).all() and
              np.all(p >= low) and np.all(p <= high))
        check(name+" optimizer terminated", fit["optimizer_success"])
        check(name+" all native pixels and nuisance counts",
              fit["ndata"] == model.ndata and fit["npar_nonlinear"] == model.npar
              and fit["npar_profiled_linear"] == 3*len(model.blocks))
        residual, jacobian, details = model.evaluate(p)
        saved = np.load(path.with_suffix(".npz"))
        close(name+" saved objective", fit["chi2"], residual@residual, atol=1e-7, rtol=1e-11)
        close(name+" saved residual", saved["residual"], residual, atol=1e-11)
        close(name+" saved Jacobian", saved["jacobian"], jacobian, atol=1e-10)
        reference, reference_details = independent_replay(model, p)
        error = float(np.max(np.abs(reference-residual)))
        objective_error = float(reference@reference-residual@residual)
        residual_errors.append(error); objective_errors.append(abs(objective_error))
        close(name+" direct-kernel independent residual", reference, residual,
              atol=2e-6, rtol=1e-8)
        close(name+" direct-kernel independent objective", reference@reference,
              fit["chi2"], atol=1e-3, rtol=1e-10)
        for block, detail, old, ref in zip(model.blocks, details, fit["per_row"], reference_details):
            pre = f'exp{block["exposure"]}_row{block["row"]}'
            check(name+" "+pre+" row provenance", old["exposure"] == block["exposure"]
                  and old["row"] == block["row"] and old["order"] == block["order"]
                  and old["trace"] == block["trace"] and old["npix"] == len(block["flux"]))
            close(name+" "+pre+" profiled coefficients", old["coefficients"], ref["beta"], atol=2e-8)
            close(name+" "+pre+" saved profile", saved[pre+"_model"], ref["model"], atol=2e-8)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        close(name+" singular values", fit["singular_values"], singular, atol=1e-9)
        check(name+" full local rank", singular[-1] > singular[0]*1e-10)
        covariance = np.linalg.inv(jacobian.T@jacobian)
        close(name+" covariance independent inverse", fit["covariance"], covariance, atol=1e-10)
        if "order_offset" in model.labels:
            j = model.labels.index("order_offset")
            close(name+" offset sign and units", fit["order_offset_m_s"], p[j]*1000)
            close(name+" conditional offset standard error", fit["order_offset_sigma_m_s"],
                  np.sqrt(covariance[j,j])*1000)
        active_now = [label for label, value, lo, hi in zip(model.labels, p, low, high)
                      if min(value-lo, hi-value) < 1e-5]
        check(name+" active boundary provenance", active_now == fit["active_bounds"])
        if active_now:
            active.append(dict(name=name, bounds=active_now))
        gradient = jacobian.T@residual
        interior = np.array([label not in active_now for label in model.labels])
        scaled = np.sqrt(np.maximum(np.diag(covariance), 0))*gradient
        gradients.append(dict(name=name,
                              maximum_interior_score_sigma_units=float(np.max(abs(scaled[interior])))))
        file_hashes[str(path.relative_to(ROOT))] = sha(path)
        file_hashes[str(path.with_suffix(".npz").relative_to(ROOT))] = sha(path.with_suffix(".npz"))
    # Reconstruct native block selection independently from NPZ metadata once per line.
    for line in [2600, 2382]:
        model = make(line, "A1", list(range(17)), {})
        for block in model.blocks:
            rec = production.META["exposures"][block["exposure"]]
            with np.load(ROOT/rec["array_file"]) as arr:
                prefix = f'{line}_row{block["row"]}'
                good = arr[prefix+"_good"]
                close(f'{line}/{prefix}/{rec["index"]} exact native pixels',
                      block["pixel"], arr[prefix+"_pixel"][good], atol=0, rtol=0)
                close(f'{line}/{prefix}/{rec["index"]} flux and error',
                      np.r_[block["flux"], block["error"]],
                      np.r_[arr[prefix+"_flux"][good], arr[prefix+"_error"][good]], atol=0, rtol=0)
    summary = json.loads((OUT/"model_summary.json").read_text())
    starts_spread, derivatives, refinements = [], [], []
    for line_s, groups in summary.items():
        line = int(line_s)
        for group, selected in groups.items():
            ns = sorted(set(row["ndata"] for row in selected.values()))
            check(f'{line}/{group} matched pixel counts', len(ns) == 1)
            for kind, best in selected.items():
                disk = json.loads((OUT/(best["name"]+".json")).read_text())
                check(f'{line}/{group}/{kind} selected result exact', best == disk)
                if group == "heldout2019_2020":
                    training = groups["train2018"][kind]
                    expected = {key: value for key, value in training["parameters"].items()
                                if not key.startswith("velocity_")}
                    check(f'{line}/{kind} heldout response exactly frozen', best["frozen"] == expected)
                    check(f'{line}/{kind} heldout disjoint data indices',
                          set(best["exposure_indices"]) == set(range(9,17)) and
                          set(training["exposure_indices"]) == set(range(9)))
                    check(f'{line}/{kind} heldout fits only velocities',
                          all(key.startswith("velocity_") for key in best["labels"]))
                    continue
                candidates = [row for _, row in rows if row["line"] == line and row["kind"] == kind
                              and row["exposure_indices"] == best["exposure_indices"] and row["optimizer_success"]]
                close(f'{line}/{group}/{kind} minimum successful start selected', best["chi2"],
                      min(row["chi2"] for row in candidates), atol=0, rtol=0)
                starts_spread.append(dict(line=line, group=group, kind=kind,
                    chi2_range=max(row["chi2"] for row in candidates)-min(row["chi2"] for row in candidates),
                    order_offset_range_m_s=(max(row["order_offset_m_s"] for row in candidates)-min(row["order_offset_m_s"] for row in candidates)) if kind.endswith("1") else None))
            if group != "heldout2019_2020":
                for simpler, fuller in [("G0","G1"),("G0","A0"),("G1","A1"),("A0","A1")]:
                    small, large = selected[simpler], selected[fuller]
                    m = make(line, fuller, large["exposure_indices"], {})
                    p0, _, _ = m.initial()
                    p = np.array([small["parameters"].get(label, value if label.startswith("width_") else 0.)
                                  for label, value in zip(m.labels, p0)])
                    close(f'{line}/{group} exact nesting {simpler}->{fuller}',
                          m.fun(p)@m.fun(p), small["chi2"], atol=1e-7, rtol=1e-10)
                    check(f'{line}/{group} nested objective {simpler}->{fuller}',
                          large["chi2"] <= small["chi2"]+1e-5)
        # Full A1 fit: all nonlinear columns independently differentiated.
        best = groups["full"]["A1"]
        model = make(line, "A1", list(range(17)), {})
        p = np.array([best["parameters"][key] for key in model.labels])
        residual, J, _ = model.evaluate(p)
        central = []
        for j, label in enumerate(model.labels):
            delta = np.zeros(len(p)); delta[j] = 1e-5
            # Independent nuisance solve, but same sampled production spline for this derivative check.
            col = (model.fun(p+delta)-model.fun(p-delta))/(2e-5)
            central.append(col)
            relative = float(np.linalg.norm(col-J[:,j])/np.linalg.norm(col))
            check(f'{line} A1 centered Jacobian {label}', relative < .003, relative)
        Jc = np.column_stack(central)
        cv = np.linalg.inv(Jc.T@Jc)
        j = model.labels.index("order_offset")
        derivatives.append(dict(line=line, conditional_offset_sigma_m_s=float(np.sqrt(cv[j,j])*1000),
                                saved_conditional_offset_sigma_m_s=best["order_offset_sigma_m_s"]))
        close(f'{line} A1 centered-J offset uncertainty', np.sqrt(cv[j,j])*1000,
              best["order_offset_sigma_m_s"], atol=.02, rtol=.002)
        # Refit with twice-finer velocity sampling, 21-point native-pixel integration,
        # and an independently differenced optimizer Jacobian.  This is a new solve,
        # separate from exact replay of the saved fit.
        refined = production.JointModel(line,"A1",step=.0125,quadrature=21)
        _, lo, hi = refined.initial()
        opt = least_squares(refined.fun, p, bounds=(lo,hi), jac="3-point",
                            diff_step=1e-5, x_scale="jac", ftol=1e-11,
                            xtol=1e-11, gtol=1e-8, max_nfev=40)
        new_r = refined.fun(opt.x)
        result = dict(line=line, success=bool(opt.success), nfev=int(opt.nfev),
                      delta_chi2=float(new_r@new_r-best["chi2"]),
                      delta_order_offset_m_s=float(1000*(opt.x[j]-p[j])),
                      parameters=dict(zip(refined.labels,opt.x.tolist())))
        refinements.append(result)
        check(f'{line} actual refined refit success', opt.success, result)
        check(f'{line} actual refined refit stable objective and offset',
              abs(result["delta_chi2"]) < .01 and abs(result["delta_order_offset_m_s"]) < .3, result)
    # Zero-asymmetry derivative is sensitive to finite-difference step because C1 != C2.
    model = make(2600,"A1",list(range(17)),{})
    p,_,_ = model.initial(); _,J,_ = model.evaluate(p)
    zero_derivative = []
    for j,label in enumerate(model.labels):
        if not label.startswith("asymmetry_"):continue
        delta=np.zeros(len(p));delta[j]=1e-5
        col=(model.fun(p+delta)-model.fun(p-delta))/(2e-5)
        relative=float(np.linalg.norm(col-J[:,j])/np.linalg.norm(col))
        zero_derivative.append(dict(label=label,relative_l2_difference=relative))
        check("a zero finite-difference derivative "+label,relative<.003,relative)
    METRICS.update(n_fits=len(rows), max_direct_kernel_residual_error=max(residual_errors),
                   max_direct_kernel_objective_error=max(objective_errors),
                   active_bound_fits=active, profile_gradients=gradients,
                   multistart_spread=starts_spread, centered_derivatives=derivatives,
                   refined_refits=refinements, zero_asymmetry_derivative=zero_derivative,
                   audited_file_hashes=file_hashes)


def profile_audit():
    import order_response_model as production
    manifest=json.loads((OUT/"profile_run_manifest.json").read_text())
    for key,path in [("code_sha256",ROOT/"code/order_response_profiles.py"),
                     ("model_summary_sha256",OUT/"model_summary.json"),
                     ("protocol_sha256",OUT/"profile_protocol.json")]:
        check("profile manifest "+key,manifest[key]==sha(path))
    profile=json.loads((OUT/"profile_summary.json").read_text())
    summary=json.loads((OUT/"model_summary.json").read_text())
    protocol=json.loads((OUT/"profile_protocol.json").read_text())
    counts=0; summaries=[]
    for line_s,groups in profile.items():
        line=int(line_s)
        for kind,points in groups.items():
            close(f"{line}/{kind} profile fixed grid",[p["order_offset_m_s"] for p in points],
                  np.array(protocol["grid_order_offset_km_s"])*1000)
            best=summary[line_s]["full"][kind]
            for point in points:
                fit=json.loads((OUT/point["source"]).read_text())
                path=OUT/point["source"]
                model=production.JointModel(line,kind,frozen=fit["frozen"])
                p=np.array([fit["parameters"][label] for label in model.labels])
                rr,J,_=model.evaluate(p);saved=np.load(path.with_suffix(".npz"));name=fit["name"]
                check(name+" correct same data and free nuisance count",
                      model.ndata==best["ndata"] and model.npar==best["npar_nonlinear"]-1)
                close(name+" frozen offset matches profile",fit["frozen"]["order_offset"]*1000,
                      point["order_offset_m_s"])
                check(name+" optimizer converged",fit["optimizer_success"] and point["optimizer_success"])
                close(name+" objective replay",rr@rr,fit["chi2"],atol=1e-7)
                close(name+" residual replay",rr,saved["residual"],atol=1e-11)
                close(name+" Jacobian replay",J,saved["jacobian"],atol=1e-10)
                close(name+" profile delta",point["delta_chi2"],fit["chi2"]-best["chi2"],atol=1e-10)
                check(name+" does not substantially improve unrestricted fit",point["delta_chi2"]>-.001)
                if abs(point["order_offset_m_s"])<1e-8:
                    simpler=summary[line_s]["full"][kind[0]+"0"]
                    close(name+" zero offset agrees with nested fit",fit["chi2"],simpler["chi2"],atol=.001,rtol=0)
                if int(round(point["order_offset_m_s"])) in [-120,0,120]:
                    ref,_=independent_replay(model,p)
                    close(name+" independent direct-kernel profile replay",ref,rr,atol=2e-6,rtol=1e-8)
                METRICS["audited_file_hashes"][str(path.relative_to(ROOT))]=sha(path)
                METRICS["audited_file_hashes"][str(path.with_suffix(".npz").relative_to(ROOT))]=sha(path.with_suffix(".npz"))
                counts+=1
            lowest=min(points,key=lambda p:p["chi2"])
            check(f"{line}/{kind} grid minimum near unrestricted minimum",
                  abs(lowest["order_offset_m_s"]-best["order_offset_m_s"])<=20)
            summaries.append(dict(line=line,kind=kind,
                 grid_minimum_offset_m_s=lowest["order_offset_m_s"],
                 unrestricted_offset_m_s=best["order_offset_m_s"],
                 zero_offset_delta_chi2=next(p["delta_chi2"] for p in points if abs(p["order_offset_m_s"])<1e-8)))
    check("all 52 frozen-offset nonlinear profiles present",counts==52)
    METRICS.update(n_profile_fits=counts,profile_summaries=summaries)


def write_result(scope):
    OUT.mkdir(parents=True, exist_ok=True)
    result = dict(created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                  pass_=all(row["pass_"] for row in CHECKS),
                  n_checks=len(CHECKS),
                  n_pass=sum(row["pass_"] for row in CHECKS),
                  n_fail=sum(not row["pass_"] for row in CHECKS),
                  scope=scope, metrics=METRICS, checks=CHECKS,
                  validation_script_sha256=sha(__file__))
    (OUT/"validation.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key: result[key] for key in ["pass_", "n_checks", "n_pass", "n_fail"]}, indent=2))
    if not result["pass_"]:
        print(json.dumps([row for row in CHECKS if not row["pass_"]], indent=2))
    return result


def write_report(result):
    metrics=result["metrics"]
    status="VERIFIED within the stated computational scope" if result["pass_"] else "REVIEW REQUIRED: at least one numerical check failed"
    lines=["# Centered order-response model: independent numerical validation", "",
           "## Material Passport", "",
           "- Mode: reproducibility and statistical interpretation validation.",
           "- Material: public native ESPRESSO HE 0515−4414 spectra, fixed shared Fe II gas template, and newly generated conditional order-response fits.",
           "- Review provenance: separate agent in the same workspace and conversation context; an additional agent independently derived the kernel moments. This is not an external or blinded audit.",
           f"- Verification status: **{status}**.",
           f"- Checks: {result['n_pass']} / {result['n_checks']} passed; {result['n_fail']} failed.",
           "- Review artifacts: `code/order_response_validate.py` and `results/order_response/validation.json`.", "",
           "The audit reconstructs the response by direct sampled-kernel convolution and independently solves each row's linear nuisance parameters with SVD least squares. Production uses shifted Gaussian-filter splines and QR profiling. Saved objectives, residuals, source bindings, matched-data nesting, parameter counts, covariance calculations, and chronological response transfer were checked. This validates the specified conditional calculation, not the physical adequacy of its noise, gas, or instrument assumptions.", "",
           "## Kernel and identifiability", "",
           "Writing s = F / 2.354820045 and d = 1.2 cbrt(a), the 80:20 mixture has area 1, mean 0, variance s², third moment 0.165888 a, and fourth moment 3s⁴ + 0.01327104 |a|^(4/3). Direct continuous quadrature verifies these identities. The [1.6, 2.6] km/s equivalent-width bounds keep the component variance strictly positive; the limiting requirement is F > 1.1303136216 km/s.", "",
           "F denotes the Gaussian-equivalent width determined by the second moment, not the actual asymmetric kernel FWHM. The parameter a controls the third moment and is not standardized skewness. The complete mixture is C¹ but generally not C² at a = 0. Its first derivative there is −0.027648 G‴. Checked finite differences are appropriate for this numerical diagnostic; a naive derivative through cbrt(0) is not.", "",
           "Zero kernel centroid removes an exact kernel translation, but does not make fitted shift and asymmetry independent. Even isolated Gaussian absorption has a shift/skew tangent correlation of magnitude √(3/5) ≈ 0.775 under uniform continuous noise. Finite windows, blends, saturation, and nuisance profiling can alter that coupling. An asymmetric LSF can move a line core while leaving its full unmasked centroid unchanged.", "",
           "## Reproduction results", "",
           f"- Main nonlinear fits audited: {metrics.get('n_fits',0)}.",
           f"- Fixed-offset nonlinear profile fits audited: {metrics.get('n_profile_fits',0)}.",
           f"- Maximum direct-kernel standardized-residual discrepancy: {metrics.get('max_direct_kernel_residual_error',float('nan')):.6g}.",
           f"- Maximum direct-kernel objective discrepancy: {metrics.get('max_direct_kernel_objective_error',float('nan')):.6g}.",
           "- Full local Jacobian rank and a separate covariance inverse were checked for every main fit; both selected full A1 fits also received all-column centered finite-difference checks.",
           "- Refinement used half the velocity-grid spacing, 21-point pixel quadrature, and a separately differenced optimizer Jacobian for actual nonlinear re-fits.", "",
           "| Transition | Refined − saved offset (m/s) | Refined − saved χ² | Refit success |",
           "|---|---:|---:|---|"]
    for row in metrics.get("refined_refits",[]):
        lines.append(f"| {row['line']} | {row['delta_order_offset_m_s']:.6g} | {row['delta_chi2']:.6g} | {row['success']} |")
    lines += ["", "Multiple starts are retained; selection uses the minimum objective among successful starts. The largest asymmetry-family main-fit spread was "
              f"{max((r['chi2_range'] for r in metrics.get('multistart_spread',[])),default=float('nan')):.6g} in χ². "
              "Optimizer termination, numerical reproducibility, and agreement among these starts do not prove a global minimum.", "",
              "## Conditional interpretation", "",
              "For Fe II 2600, the selected Gaussian response gives an order difference −62.75 ± 16.66 m/s; adding the centered asymmetric response gives −58.64 ± 18.65 m/s. These are conditional local standard errors under diagonal native ERRDATA and a fixed gas template. The zero-offset versus free-offset improvement decreases from 14.21 to 10.86 in χ². Thus this particular response family does not erase the order difference. It also does not establish that the residual difference is beyond conventional instrument, extraction, wavelength calibration, or gas-model effects.", "",
              "The full 2600 A0 fit and several training A0 starts touch an asymmetry bound. This restricts interpretation of its conditional profile; shape bounds are assumptions rather than measured calibration limits. The Fe II 2382 control remains much closer to zero and has a weak order-offset preference. Full minimum reduced objectives are about 1.34 (2600) and 1.29 (2382), so model/noise adequacy is not established by these fits.", "",
              "The local A1 shape estimates are weak: for 2600, a52 = −0.353 ± 0.516 and a53 = −0.925 ± 0.775; for 2382, a41 = −0.341 ± 0.990 and a42 = −0.044 ± 1.637. These conditional local errors ignore the hard [−1, 1] bounds and are not calibrated confidence intervals. In particular, an inactive bound flag does not imply that the nearby boundary is inferentially irrelevant. The shapes are not empirical LSF measurements.", "",
              "Training on 2018 exposures and freezing response parameters on 2019–2020 is verified exactly, with disjoint exposure indices. For 2600, A1 yields held-out χ² larger by 5.15 than G1. This is a conditional response-transfer check, not independent prediction: the fixed gas template was estimated from the full coadd, including these held-out photons. It also supplies no second cosmological-age sample.", "",
              "An astrophysical frequency displacement of one transition is common to its simultaneous measurements in both orders. A difference between orders in the same exposure is therefore a useful conventional-effect diagnostic, but this sensitivity experiment does not isolate which conventional mechanism produces it. No cosmic-time-law or varying-constant inference is calibrated here.", "",
              "## Statistical fallacy scan — 11 / 11 considered", "",
              "| Item | Scope-specific assessment |",
              "|---|---|",
              "| Simpson's paradox | Shared response is an aggregate assumption. Prior exposure/chronological contrasts and the held-out loss are retained; pooled results do not establish identical exposure-level response. |",
              "| Ecological fallacy | No inference from a pooled order diagnostic to individual atomic evolution or multiple cosmological ages. |",
              "| Berkson/selection bias | This target and transition were chosen after the earlier order discrepancy; the diagnostic is conditional on that selection. |",
              "| Collider bias | No causal adjustment identification is claimed; nuisance profiling cannot establish a causal LSF mechanism. |",
              "| Base-rate neglect | Not a diagnostic classifier; no posterior discovery probability inferred from an objective difference. |",
              "| Regression to the mean | Reusing the selected discrepancy motivates a control transition and retained chronological transfer results; shared-template reuse remains explicit. |",
              "| Survivorship bias | All planned main starts and grid points are retained, including boundary fits; the minimum-successful rule and total counts are audited. |",
              "| Look-elsewhere effect | Exploratory model family and previous searches preclude treating local objective improvement as a calibrated discovery probability. |",
              "| Forking paths | Protocol frozen before these model comparisons, but after seeing earlier data. Assumed mixture weights/amplitude and all four models are disclosed. |",
              "| Correlation versus causation | Better fit or same-order direction does not identify physical or instrumental cause. |",
              "| Reverse causality | No causal or temporal ordering is inferred; observing years are not different source cosmic ages. |", "",
              "## Re-run", "", "```bash",
              "OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python riemann_clock/experiments/highz_feasibility_2026-09-20/code/order_response_validate.py",
              "```", "",
              "The JSON records exact source/result hashes, checks, numerical refinement metrics, active bounds, and multistart spreads. Prior raw FITS validation is linked by artifact hashes rather than represented as a new raw download or fresh independent astrophysical dataset.", ""]
    (ROOT/"reports/order_response_validation.md").write_text("\n".join(lines))


if __name__ == "__main__":
    kernel_audit()
    if (OUT/"model_run_manifest.json").exists():
        production_audit()
        if (OUT/"profile_run_manifest.json").exists():
            profile_audit()
        result=write_result("Independent continuous kernel, direct convolution and least-squares replay, matched model nesting, all saved native-pixel fits, source binding, conditional response transfer, centered derivatives, and refined nonlinear re-fits. Does not establish empirical LSF calibration, extraction covariance, astrophysical independence, or new physics.")
        write_report(result)
    else:
        write_result("Continuous proposed LSF kernel only; production model and fits not yet audited.")
