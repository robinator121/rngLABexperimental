"""Tables of the B2 numerical-limits investigation, from its result file only.

    PYTHONPATH=src .venv/bin/python research/v020-nuclear-b2-numerics/report_tables.py

Reads results/numerics.json and writes results/tables.md. Computes nothing new: it formats recorded
values (maxima over rows are taken where a table would otherwise be too long, and are labelled).
"""
from __future__ import annotations

import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def fmt(value) -> str:
    if value is None:
        return "-"
    if isinstance(value, str):
        return value
    value = float(value)
    if value == 0:
        return "0"
    return f"{value:.3g}" if 1e-3 <= abs(value) < 1e4 else f"{value:.1e}"


def main() -> None:
    r = json.loads((HERE / "results" / "numerics.json").read_text())
    out = ["# B2 numerical limits: tables from `results/numerics.json`", "",
           f"Code `{r['code']['git_commit']}` (dirty: {r['code']['git_dirty']}); started {r['started_utc']}; "
           f"{r['elapsed_seconds']:.0f} s. Python {r['environment']['python']}, NumPy {r['environment']['numpy']}, "
           f"SciPy {r['environment']['scipy']}.", ""]

    out += ["## 1. Binomial distribution function against exact sums (60-digit decimal), n up to 1e9", "",
            "Largest absolute error over nine points from -6 to +6 standard deviations.", "",
            "| n | p | Boost F | Boost S | cephes `bdtr` (not used) | as consulted by the sampler |",
            "|---:|---:|---:|---:|---:|---:|"]
    for row in r["cdf_exact"]["rows"]:
        out.append(f"| {row['n']:.0e} | {row['p']} | {fmt(row['max_abs_error_F_boost'])} | "
                   f"{fmt(row['max_abs_error_S_boost'])} | {fmt(row['max_abs_error_cephes_bdtr'])} | "
                   f"{fmt(row['max_abs_error_sampler_tails'])} |")

    out += ["", "## 2. Large n with small means (reference: 50-digit pmf recurrence)", "",
            f"The sampler uses the recurrence table for means up to {r['cdf_small_mean']['table_mean_max']:g} "
            "and Boost above.", "",
            "| n | mean | regime | Boost F | Boost S | recurrence table | as consulted by the sampler |",
            "|---:|---:|---|---:|---:|---:|---:|"]
    for row in r["cdf_small_mean"]["rows"]:
        out.append(f"| {row['n']:.0e} | {row['mean']} | {row['regime']} | {fmt(row['max_abs_error_F_boost'])} | "
                   f"{fmt(row['max_abs_error_S_boost'])} | "
                   f"{fmt(row['max_abs_error_table']) if row['regime'] == 'table' else '-'} | "
                   f"{fmt(row['max_abs_error_sampler_tails'])} |")

    out += ["", "## 3. n from 1e12 to 2^53, p = 1/2: exact values (centre value 1/2 plus exact pmf sums)", "",
            "| n | z | exact F | Boost F error | Boost S error | error of the tail the sampler consults | Edgeworth error |",
            "|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["cdf_huge_n_exact"]["rows"]:
        out.append(f"| {row['n']:.4g} | {row['z']} | {row['exact_F']:.6f} | {fmt(row['boost_F_error'])} | "
                   f"{fmt(row['boost_S_error'])} | {fmt(row['sampler_consulted_error'])} | "
                   f"{fmt(row['edgeworth_error'])} |")
    out += ["", "Centre anchor and Edgeworth comparison (error of the reference is of order 1/n):", "",
            "| n | Boost F error at the centre | Boost S error at the centre | max abs(Boost F - Edgeworth), p = 0.5 | p = 0.1 | "
            "ms per evaluation |", "|---:|---:|---:|---:|---:|---:|"]
    for row in r["cdf_large_n"]["rows"]:
        gap = row["max_abs_boost_minus_edgeworth"]
        out.append(f"| {row['n']:.4g} | {fmt(row['centre_error_F'])} | {fmt(row['centre_error_S'])} | "
                   f"{fmt(gap['0.5'])} | {fmt(gap['0.1'])} | {1e3 * row['seconds_per_evaluation']:.2f} |")

    inv, tails, traced = r["inversion"], r["inversion_tails"], r["cdf_finiteness"]["traced"]
    out += ["", "## 4. Inversion against exact arithmetic", "",
            f"- Grid of uniforms, n up to 55, five values of p: {inv['grid_mismatches']} mismatches in "
            f"{inv['grid_draws']} draws against exact rational inversion at the cell midpoint.",
            f"- Uniforms 1e-12 either side of every distribution-function step: {inv['near_step_mismatches']} "
            f"mismatches in {inv['near_step_draws']} draws.",
            f"- Extreme uniforms (0, 2^-53, 1e-15, ..., 1 - 2^-53), n up to 1e9, both regimes: "
            f"{tails['total_mismatches']} mismatches in {sum(len(row['u']) for row in tails['rows'])} draws "
            "against exact decimal quantiles.", "",
            f"- Uniforms adjacent to every distribution-function step (n = 1..100 at four p, and three Boost-regime "
            f"cases): {r['table_boundaries']['adjacent_mismatches']} of {r['table_boundaries']['adjacent_draws']} "
            f"draws differ from exact inversion (the evaluated F puts the step on the other side of the midpoint); "
            f"{r['table_boundaries']['far_mismatches']} of {r['table_boundaries']['far_draws']} draws eight cells away do.",
            f"- Exact symmetry anchors (odd n, p = 1/2, n up to PART_MAX - 1): all exact: "
            f"{r['inversion_anchor']['anchors_ok']}; differences between brackets of half-width 1, 4 and 64: "
            f"{r['inversion_anchor']['bracket_differences']}.",
            f"- Boost's F and S at every point the bisection can reach up to PART_MAX: "
            f"{r['cdf_finiteness']['nonfinite_within_domain']} non-finite of "
            f"{r['cdf_finiteness']['points_within_domain']} points; beyond PART_MAX (near 2^53): "
            f"{r['cdf_finiteness']['nonfinite_beyond_domain']} non-finite (PART_MAX = "
            f"{r['cdf_finiteness']['part_max']}).",
            f"- Every evaluation made by {traced['draws']} traced draws (n from 1e6 to PART_MAX, extreme and central "
            f"uniforms, narrow brackets that force the fallbacks): {traced['evaluations']['F']} of F and "
            f"{traced['evaluations']['S']} of S, {traced['nonfinite']['F'] + traced['nonfinite']['S']} non-finite; "
            f"standardised points reached from {traced['z_min']:.3g} to {traced['z_max']:.3g}; draws refused: "
            f"{traced['refused_draws']}.",
            f"- Mass beyond the recurrence table: at most 10^{r['table_truncation']['max_log10_omitted_mass']:.1f}.", "",
            "## 5. The discretisation of the midpoint rule (exact F; exact decimal evaluation)", "",
            "| n | p | mean | z at u = 0 | z at u = 1 - 2^-53 | relative variance bound (with tails) | exact relative "
            "variance error | old bound without tails | variance added by inverting at u itself |",
            "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for row in r["midpoint_convention"]["rows"]:
        out.append(f"| {row['n']:.4g} | {row['p']:.3g} | {row['mean']:.3g} | {fmt(row['z_at_u_0'])} | "
                   f"{fmt(row['z_at_u_max'])} | {fmt(row['relative_variance_bound'])} | "
                   f"{fmt(row['relative_variance_error_exact'])} | {fmt(row['old_bound_without_tails'])} | "
                   f"{fmt(row['left_endpoint_share'])} |")
    out.append("")
    out.append(f"The bound holds in every row: {r['midpoint_convention']['bound_holds_everywhere']}.")
    out.append("")

    rel = [row["relative_error"] for row in r["small_p"]["rows"] if row["relative_error"] is not None]
    out += ["", "## 6. Very small probabilities", "",
            f"P(X >= 1) = betainc(1, n, p) against 80-digit decimal for n in 1, 1e6, 1e12, 2^53 and p from 1e-300 "
            f"to 1e-5 ({len(rel)} cases): largest relative error {fmt(max(rel))}.", "",
            "## 7. Interval laws (decimal closed-form Bateman)", "",
            "Close decay constants (four-member chain, half-lives 1 + i x spacing, interval 1):", "",
            "| relative spacing | digits used | max abs(decimal - matrix exponential) | float64 closed form (Layer A `bateman_chain`) | "
            "brackets hold |", "|---:|---:|---:|---|---|"]
    for row in r["paths_close"]["rows"]:
        closed = row["float_closed_form_deepest"]
        closed = "refused" if isinstance(closed, str) else f"{closed:.12g} (decimal: {row['deepest_path_probability']:.12g})"
        out.append(f"| {row['relative_spacing']:.0e} | {row['digits_used']} | "
                   f"{fmt(row['max_abs_decimal_minus_expm'])} | {closed} | {row['bounds_hold']} |")
    extreme = r["paths_extreme"]["rows"]
    out += ["", f"Extreme rate ratios (half-lives from 1e-12 to 1e12, intervals from 1e-15 to 1e15, {len(extreme)} "
                f"cases): largest abs(decimal - matrix exponential) {fmt(max(x['max_abs_decimal_minus_expm'] for x in extreme))}; "
                f"largest abs(sum - 1) {fmt(max(abs(x['sum_minus_one']) for x in extreme))}; most digits used "
                f"{max(x['digits_used'] for x in extreme)}; brackets hold in all: {all(x['bounds_hold'] for x in extreme)}.",
            "", "Tiny intervals (S2 from Xa-400, deepest path of four decays):", "",
            "| interval | path probability (decimal) | relative error of the float64 matrix exponential |",
            "|---:|---:|---:|"]
    for row in r["paths_tiny"]["rows"]:
        out.append(f"| {row['delta']:.0e} | {row['deepest_path_probability']:.4e} | "
                   f"{fmt(row['expm_relative_error_on_it'])} |")
    stiff = r["expm_stiff"]["rows"]
    out += ["", f"Stiff, close-rate chain S5: largest abs(float64 matrix exponential - 150-digit decimal) "
                f"{fmt(max(x['max_abs_expm_minus_decimal'] for x in stiff))}; largest abs(default-precision decimal - "
                f"150-digit decimal) {fmt(max(x['decimal_default_vs_150_digits'] for x in stiff))}.",
            "", "## 8. numpy's binomial at large n (the calibration's null; not used by the B2 sampler)", "",
            f"Each row: {r['numpy_binomial']['rows'][0]['draws']:,} draws of numpy's `binomial(n, p)`; mean z is "
            "(sample mean - np) / its standard error; the variance ratio is the sample variance over np(1 - p), "
            f"with standard error {r['numpy_binomial']['variance_ratio_standard_error']:.4f} under the binomial law.", "",
            "| n | p | expected count | mean z | variance ratio | (ratio - 1) / se |",
            "|---:|---:|---:|---:|---:|---:|"]
    se = r["numpy_binomial"]["variance_ratio_standard_error"]
    for row in r["numpy_binomial"]["rows"]:
        out.append(f"| {row['n']:,} | {row['p']:.6g} | {row['expected_count']:.4g} | {row['mean_z']:.2f} | "
                   f"{row['variance_ratio']:.4f} | {(row['variance_ratio'] - 1) / se:.1f} |")
    call = r["numpy_binomial"]["smallest_call_of_the_big_null"]
    out += ["", f"Smallest expected count among the BIG null's binomial calls with n above 1e13: {call['expected']:.4g} "
                f"(n = {call['n']:.6g}, p = {call['p']:.12g}; interval {call['interval']}, {call['nuclide']}, "
                f"{call['call']} category {call['category']}, initial size {call['size']}). The last row above is "
                "that call's n and p.",
            "", "## 9. Populations from 1e6 to 2^62 (S2, six intervals, one run each; descriptive)", "",
            "| initial nuclei | parts | seconds | words read | conserved | max abs z against exact means (54 values) |",
            "|---:|---:|---:|---:|---|---:|"]
    for row in r["populations"]["rows"]:
        out.append(f"| {row['n0']:.4g} | {row['parts']} | {row['seconds']:.2f} | {row['words_read']} | "
                   f"{row['conserved']} | {row['max_abs_z_against_exact_means']:.2f} |")
    out += ["", f"A total population above 2^63 - 1 is refused: {r['populations']['population_above_int64_refused']}.", ""]
    (HERE / "results" / "tables.md").write_text("\n".join(out))
    print(f"wrote {HERE / 'results' / 'tables.md'}")


if __name__ == "__main__":
    main()
