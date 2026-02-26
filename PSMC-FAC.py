#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import argparse
import math

import numpy as np
import pandas as pd

from scipy.spatial.distance import directed_hausdorff
import similaritymeasures
import matplotlib.pyplot as plt


NUM_PSMC_RE = re.compile(r"\.(\d+)\.psmc$")  # coverage nomenclature


# =====================================================
# 0) logspace and build years
# =====================================================

def logspace(limit, start, n):
    '''
    Returns a vector of logarithmically spaced integers. limit and start are +1
    '''
    result = [start]
    n = n + 1
    if n > 1:
        ratio = (float(limit) / result[-1]) ** (1.0 / (n - len(result)))
    while len(result) < n:
        next_value = result[-1] * ratio
        if next_value - result[-1] >= 1:
            result.append(next_value)
        else:
            result.append(result[-1] + 1)
            ratio = (float(limit) / result[-1]) ** (1.0 / (n - len(result)))
    return np.array(list(map(lambda x: round(x) - 1, result)), dtype=np.uint64)


def build_v_years(tmin, tmax, n_timepoints):
    '''
    Build the custom vector of years
    '''
    return logspace(limit=int(tmax) + 1, start=int(tmin) + 1, n=int(n_timepoints) - 1)


# =====================================================
# Coverage table (optional)
# =====================================================

def load_coverage_table(path):
    '''
    User can specify a table with coverages
    a tsv file (tab separated)

    Mandatory columns:
    - psmc_file
    - coverage

    coverage can be:
    - numeric
    - baseline
    '''
    if path is None:
        return None

    df = pd.read_csv(path, sep="\t")

    if "psmc_file" not in df.columns or "coverage" not in df.columns:
        raise ValueError("Coverage table must contain columns: psmc_file and coverage")

    mapping = {}

    for _, row in df.iterrows():
        fname = str(row["psmc_file"]).strip()
        cov = row["coverage"]

        if pd.isna(cov):
            raise ValueError("Coverage is missing for file: %s" % fname)

        cov_str = str(cov).strip().lower()

        if cov_str == "baseline":
            mapping[fname] = "baseline"
        else:
            try:
                mapping[fname] = float(cov)
            except:
                raise ValueError("Coverage must be numeric or 'baseline' for file: %s" % fname)

    return mapping


def list_all_psmc_filenames(root):
    found = []
    for dirpath, _, filenames in os.walk(root):
        for f in filenames:
            if f.endswith(".psmc"):
                found.append(f)
    return sorted(found)


# =====================================================
# 1) Parse PSMC and coverage. FNR correction
# =====================================================

def PSMC_last_iteration(psmc_path):
    '''
    psmc does more than one iteration of the EM algorithm. Just get the last one
    '''
    with open(psmc_path) as n:
        last = n.read().split("//")[-2]
    return last


def PSMC_to_FNR_dataframe(psmc_path, mu, g, FNR_min, FNR_max, svalue, v_years):
    '''
    We parse the file
    '''
    svalue = int(svalue)
    FNR_max = float(FNR_max)
    FNR_min = float(FNR_min)
    g = float(g)
    mu = float(mu)

    psmc_last = PSMC_last_iteration(psmc_path=psmc_path)
    all_file = psmc_last.split("\n")[1:-1]

    times = [line.split("\t")[2] for line in all_file if line.startswith("RS")]
    lambdas = [line.split("\t")[3] for line in all_file if line.startswith("RS")]
    theta_0 = float([line.split("\t")[1] for line in all_file if line.startswith("TR")][0])

    FNR_vector = np.arange(FNR_min, FNR_max + 0.01, 0.01)
    FNR_vector = [round(FNR, 2) for FNR in FNR_vector]

    FNR_psmc_dict = {}

    for FNR in FNR_vector:
        # We consider each FNR
        theta_0_FNR = theta_0 / (1 - FNR) ## theta correction
        N0 = theta_0_FNR / (4 * float(mu) * svalue)
        generations = [2 * N0 * float(t) for t in times]
        years = [gen * g for gen in generations]
        N = [N0 * float(lmd) for lmd in lambdas]

        psmc_in = pd.DataFrame({"years": years, "N": N})

        nans = np.empty(len(v_years))
        nans[:] = np.nan
        x_df = pd.DataFrame({"years": v_years, "N": nans})

        psmc_in = pd.concat([psmc_in, x_df])
        psmc_in = psmc_in.sort_values(by=["years"])
        psmc_in = psmc_in.fillna(method="ffill")

        psmc_sampled = psmc_in[psmc_in.years.isin(v_years)]
        psmc_sampled = psmc_sampled.drop_duplicates(subset=["years"])

        this_psmc_FNR = list(psmc_sampled["N"])
        FNR_psmc_dict[FNR] = this_psmc_FNR

    FNR_psmc_dict["years"] = list(psmc_sampled["years"])
    df_FNR_psmc = pd.DataFrame(FNR_psmc_dict)
    cols = ["years"] + [c for c in df_FNR_psmc.columns if c != "years"]
    df_FNR_psmc = df_FNR_psmc[cols]
    return df_FNR_psmc


def write_psmc_tsv(psmc_path, df, outdir):
    '''
    Get a tsv with all the fnr corrections
    '''
    name = os.path.basename(psmc_path).split(".psmc")[0] + ".tsv"
    outpath = os.path.join(outdir, name)
    df.to_csv(outpath, sep="\t", index=False)
    return outpath


# =====================================================
# Sample + coverage extraction (coverage nomenclature)
# =====================================================

def sample_from_filename(filename):
    m = NUM_PSMC_RE.search(filename)
    if m:
        return filename[:m.start()]
    if filename.endswith(".psmc"):
        return filename[:-5]
    raise ValueError("Not a .psmc: " + filename)


def coverage_from_filename(filename):
    m = NUM_PSMC_RE.search(filename)
    if m:
        return float(m.group(1))
    if filename.endswith(".psmc"):
        return "baseline"
    raise ValueError("Not a .psmc: " + filename)


def get_base_set(base_files, coverage_mapping):
    '''
    Decide references from:
    1) coverage table, if provided
    2) base_files argument otherwise
    '''
    if coverage_mapping is not None:
        base_set = set()
        for fname, cov in coverage_mapping.items():
            if cov == "baseline":
                base_set.add(fname)
        return base_set

    return set(base_files)


def get_coverage(filename, coverage_mapping, base_set):
    '''
    Returns:
    - "baseline" for references
    - float for downsampled files
    '''
    if coverage_mapping is None:
        if filename in base_set:
            return "baseline"
        return coverage_from_filename(filename)

    if filename not in coverage_mapping:
        raise ValueError("File %s not found in coverage_table" % filename)

    table_value = coverage_mapping[filename]

    if filename in base_set:
        if table_value != "baseline":
            raise ValueError("Reference file %s must be tagged as 'baseline' in coverage_table" % filename)
        return "baseline"

    if table_value == "baseline":
        raise ValueError("Non-reference file %s cannot be tagged as 'baseline' in coverage_table" % filename)

    return float(table_value)


# =====================================================
# 2) H anf F Distances + SSE
# =====================================================

def hausdorff_symmetric(A, B):
    d1 = directed_hausdorff(A, B)[0]
    d2 = directed_hausdorff(B, A)[0]
    return float(max(d1, d2))


def frechet_dist(A, B):
    return float(similaritymeasures.frechet_dist(A, B))


def compute_sse(expected, predicted):
    return float(np.sum((predicted - expected) ** 2))


def compute_sse_log10(expected, predicted):
    mask = (expected > 0) & (predicted > 0)
    if not np.any(mask):
        return np.nan
    return float(np.sum((np.log10(predicted[mask]) - np.log10(expected[mask])) ** 2))


def coords_distance_transform(years, ne, mode):
    years = years.astype(float)
    ne = ne.astype(float)

    if mode == "linear": # no scaling
        x = years
        y = ne
    elif mode == "loglog":
        eps = 1e-300 # avoid zeros
        x = np.log10(np.maximum(years, eps))
        y = np.log10(np.maximum(ne, eps))
    elif mode == "psmc": # the same as used in PSMC plots. preferred for empirists. 
        eps = 1e-300 # avoid zeros
        x = np.log10(np.maximum(years, eps))
        y = ne * 1e4
    else:
        raise ValueError("Unknown distance_scale: " + str(mode))

    return np.column_stack([x, y])


def coords_raw_from_df(df, fnr_value): # just get the coordinates from the dataframe 
    years = df["years"].to_numpy(dtype=float)
    ne = df[fnr_value].to_numpy(dtype=float)
    return years, ne


# =====================================================
# 3) Curves (polynomials) + plot
# =====================================================

def get_optimal_rows(df, metric):
    '''
    Separates by sample and coverage and gets the row where the selected metric is min
    '''
    idx = df.groupby(["Sample", "Coverage"])[metric].idxmin()
    return df.loc[idx].sort_values(["Sample", "Coverage"]).reset_index(drop=True)


def fit_poly2(x, y):
    '''
    Fits degree 2 polynom
    '''
    coeffs = np.polyfit(x, y, 2)
    poly = np.poly1d(coeffs)
    return coeffs, poly


def write_curves_and_plot(big_df, out_prefix, out_pdf, out_txt):
    # Get tables by H and F
    df_h = get_optimal_rows(big_df, "Hausdorff_Distance")
    df_f = get_optimal_rows(big_df, "Frechet_Distance")
    
    # Provide the csv

    df_h.to_csv(out_prefix + "_hausdorff.csv", index=False)
    df_f.to_csv(out_prefix + "_frechet.csv", index=False)

    samples = sorted(set(big_df["Sample"].unique()))

    # check if makes sense to do curves

    any_three_points = False
    for sample in samples:
        sdf = df_h[df_h["Sample"] == sample]
        if len(sdf) >= 3: # only if there are more than 3 samples
            any_three_points = True
            break

    if not any_three_points:
        with open(out_txt, "w") as f:
            f.write("Curves skipped: no sample has >=3 coverage points.\n")
        return

    ncols = 4
    nrows = int(math.ceil(len(samples) / float(ncols))) if samples else 1

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(4*ncols, 3*nrows), sharex=True, sharey=True)
    axes = np.array(axes).flatten()

    txt_lines = []

    for ax, sample in zip(axes, samples):
        txt_lines.append("Sample: %s\n" % sample)

        sdf_h = df_h[df_h["Sample"] == sample].sort_values("Coverage")
        xh = sdf_h["Coverage"].values.astype(float)
        yh = sdf_h["FNR_Value"].values.astype(float)

        if len(xh) < 3:
            ax.scatter(xh, yh, s=15, label="Hausdorff (1 point)")
            txt_lines.append("  Approach: Hausdorff\n")
            txt_lines.append("  Curves: skipped (needs >= 3 coverages; found 1)\n\n")
        else:
            ch, ph = fit_poly2(xh, yh)
            xfit = np.linspace(xh.min(), xh.max(), 200)
            ax.scatter(xh, yh, s=15)
            ax.plot(xfit, ph(xfit), label="Hausdorff")
            txt_lines.append("  Approach: Hausdorff\n")
            txt_lines.append("  Equation: y=%.2e x^2 + %.2e x + %.2e\n\n" % (ch[0], ch[1], ch[2]))

        sdf_f = df_f[df_f["Sample"] == sample].sort_values("Coverage")
        xf = sdf_f["Coverage"].values.astype(float)
        yf = sdf_f["FNR_Value"].values.astype(float)

        if len(xf) < 3:
            ax.scatter(xf, yf, s=15, label="Frechet (1 point)")
            txt_lines.append("  Approach: Frechet\n")
            txt_lines.append("  Curves: skipped (needs >= 3 coverages; found 1)\n\n")
        else:
            cf, pf = fit_poly2(xf, yf)
            xfit2 = np.linspace(xf.min(), xf.max(), 200)
            ax.scatter(xf, yf, s=15)
            ax.plot(xfit2, pf(xfit2), linestyle="--", label="Frechet")
            txt_lines.append("  Approach: Frechet\n")
            txt_lines.append("  Equation: y=%.2e x^2 + %.2e x + %.2e\n\n" % (cf[0], cf[1], cf[2]))

        ax.set_title(sample, fontsize=9)
        ax.legend(fontsize=7, frameon=False)
        ax.grid(alpha=0.3)

    for ax in axes[len(samples):]:
        ax.axis("off")

    fig.text(0.5, 0.04, "Coverage", ha="center")
    fig.text(0.04, 0.5, "Best FNR", va="center", rotation="vertical")
    fig.tight_layout(rect=[0.03, 0.06, 0.97, 0.97])

    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)

    with open(out_txt, "w") as f:
        f.writelines(txt_lines)


# =====================================================
# 5) Trajectory plots
# =====================================================

def plot_trajectories_by_sample(root, base_files, coverage_mapping, mu, g, svalue, v_years,
                               FNR_min, FNR_max, big_df, outdir_plots):
    grp = big_df.groupby(["Sample", "Coverage"])
    idx_h = grp["Hausdorff_Distance"].idxmin()
    idx_f = grp["Frechet_Distance"].idxmin()
    best_h = big_df.loc[idx_h].set_index(["Sample", "Coverage"])
    best_f = big_df.loc[idx_f].set_index(["Sample", "Coverage"])

    base_set = get_base_set(base_files, coverage_mapping)

    baseline_cache = {}
    ds_cache = {}

    for dirpath, _, filenames in os.walk(root):
        psmcs = [f for f in filenames if f.endswith(".psmc")]
        if not psmcs:
            continue

        baselines_here = [f for f in psmcs if f in base_set]
        if not baselines_here:
            continue

        if coverage_mapping is None:
            downs_here = [f for f in psmcs if NUM_PSMC_RE.search(f)]
        else:
            downs_here = [f for f in psmcs if f in coverage_mapping]

        for base in sorted(baselines_here):
            base_path = os.path.join(dirpath, base)
            sample = sample_from_filename(base)

            key = (dirpath, sample)
            if key not in baseline_cache:
                df_base = PSMC_to_FNR_dataframe(base_path, mu, g, 0.0, 0.0, svalue, v_years)
                baseline_cache[key] = df_base
            else:
                df_base = baseline_cache[key]

            years_base, ne_base = coords_raw_from_df(df_base, 0.0)

            matching = []
            for f in downs_here:
                cov = get_coverage(f, coverage_mapping, base_set)
                if cov == "baseline":
                    continue
                matching.append((f, cov))

            if not matching:
                continue

            matching = sorted(matching, key=lambda x: float(x[1]))

            pdf_path = os.path.join(outdir_plots, "%s.trajectories.pdf" % sample)

            n = len(matching)
            ncols = 3
            nrows = n

            fig, axes = plt.subplots(
                nrows=nrows,
                ncols=ncols,
                figsize=(4 * ncols, 3 * nrows),
                sharex=True,
                squeeze=False
            )

            for row_idx, (fname, cov) in enumerate(matching):
                ds_path = os.path.join(dirpath, fname)

                if ds_path not in ds_cache:
                    df_ds = PSMC_to_FNR_dataframe(ds_path, mu, g, FNR_min, FNR_max, svalue, v_years)
                    ds_cache[ds_path] = df_ds
                else:
                    df_ds = ds_cache[ds_path]

                years_ds = df_ds["years"].to_numpy(dtype=float)
                ne_unc = df_ds[0.0].to_numpy(dtype=float) if 0.0 in df_ds.columns else None

                k = (sample, float(cov))
                fnr_h = float(best_h.loc[k]["FNR_Value"]) if k in best_h.index else None
                fnr_f = float(best_f.loc[k]["FNR_Value"]) if k in best_f.index else None

                ax0 = axes[row_idx, 0]
                ax1 = axes[row_idx, 1]
                ax2 = axes[row_idx, 2]

                # Panel 1: no correction
                ax0.plot(years_base, ne_base, label="Baseline")
                if ne_unc is not None:
                    ax0.plot(years_ds, ne_unc, label="Downsample (FNR=0)")
                ax0.set_xscale("log")
                ax0.set_title("%s | cov=%s | no correction" % (sample, str(cov)), fontsize=9)
                ax0.grid(alpha=0.3)
                ax0.legend(fontsize=6, frameon=False)

                # Panel 2: best Hausdorff correction
                ax1.plot(years_base, ne_base, label="Baseline")
                if fnr_h is not None and fnr_h in df_ds.columns:
                    ax1.plot(years_ds, df_ds[fnr_h].to_numpy(dtype=float),
                             label="Corrected (best H, FNR=%.2f)" % fnr_h)
                ax1.set_xscale("log")
                ax1.set_title("%s | cov=%s | Hausdorff" % (sample, str(cov)), fontsize=9)
                ax1.grid(alpha=0.3)
                ax1.legend(fontsize=6, frameon=False)

                # Panel 3: best Frechet correction
                ax2.plot(years_base, ne_base, label="Baseline")
                if fnr_f is not None and fnr_f in df_ds.columns:
                    ax2.plot(years_ds, df_ds[fnr_f].to_numpy(dtype=float),
                             label="Corrected (best F, FNR=%.2f)" % fnr_f)
                ax2.set_xscale("log")
                ax2.set_title("%s | cov=%s | Frechet" % (sample, str(cov)), fontsize=9)
                ax2.grid(alpha=0.3)
                ax2.legend(fontsize=6, frameon=False)

            fig.tight_layout()
            fig.savefig(pdf_path, format="pdf")
            plt.close(fig)


# =====================================================
# MAIN PIPELINE
# =====================================================

def run_pipeline(args):
    root = os.path.abspath(args.root)
    outdir = os.path.abspath(args.outdir)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    outdir_psmc_tsv = os.path.join(outdir, "1_psmc_tsv")
    outdir_plots = os.path.join(outdir, "5_trajectory_plots")

    if args.write_psmc_tsv and not os.path.exists(outdir_psmc_tsv):
        os.makedirs(outdir_psmc_tsv)

    if not args.no_trajectory_plots and not os.path.exists(outdir_plots):
        os.makedirs(outdir_plots)

    coverage_mapping = load_coverage_table(args.coverage_table)

    if coverage_mapping is None and args.base_files is None:
        raise SystemExit("ERROR: you must provide either --coverage_table or --base_files.")

    base_files = args.base_files if args.base_files is not None else []
    base_set = get_base_set(base_files, coverage_mapping)

    if len(base_set) == 0:
        raise SystemExit("ERROR: no baseline files were identified.")

    all_psmc = list_all_psmc_filenames(root)

    if coverage_mapping is not None:
        table_names = set(coverage_mapping.keys())
        found_names = set(all_psmc)

        missing_in_table = sorted(list(found_names - table_names))
        extra_in_table = sorted(list(table_names - found_names))

        if missing_in_table:
            raise SystemExit(
                "ERROR: --coverage_table is missing %d PSMC files found under --root.\nMissing:\n%s"
                % (len(missing_in_table), "\n".join(missing_in_table))
            )

        if extra_in_table:
            raise SystemExit(
                "ERROR: --coverage_table contains %d PSMC files that were not found under --root.\nExtra:\n%s"
                % (len(extra_in_table), "\n".join(extra_in_table))
            )

    v_years = build_v_years(tmin=args.tmin, tmax=args.tmax, n_timepoints=args.n_timepoints)

    big_rows = []

    for dirpath, _, filenames in os.walk(root):
        psmcs = [f for f in filenames if f.endswith(".psmc")]
        if not psmcs:
            continue

        baselines_here = [f for f in psmcs if f in base_set]
        if not baselines_here:
            continue

        if coverage_mapping is None:
            downs_here = [f for f in psmcs if NUM_PSMC_RE.search(f)]
        else:
            downs_here = [f for f in psmcs if f in coverage_mapping]

        for base in sorted(baselines_here):
            base_path = os.path.join(dirpath, base)
            sample = sample_from_filename(base)

            df_base = PSMC_to_FNR_dataframe(base_path, args.mu, args.g, 0.0, 0.0, args.svalue, v_years)
            if args.write_psmc_tsv:
                write_psmc_tsv(base_path, df_base, outdir_psmc_tsv)

            years_base, ne_base = coords_raw_from_df(df_base, 0.0)
            base_xy = coords_distance_transform(years_base, ne_base, args.distance_scale)

            matching = []
            for f in downs_here:
                cov = get_coverage(f, coverage_mapping, base_set)
                if cov == "baseline":
                    continue
                matching.append((f, cov))
            matching = sorted(matching, key=lambda x: float(x[1]))

            for fname, cov in matching:
                ds_path = os.path.join(dirpath, fname)

                df_ds = PSMC_to_FNR_dataframe(ds_path, args.mu, args.g, args.FNR_min, args.FNR_max, args.svalue, v_years)
                if args.write_psmc_tsv:
                    write_psmc_tsv(ds_path, df_ds, outdir_psmc_tsv)

                for col in df_ds.columns:
                    if col == "years":
                        continue
                    fnr_val = float(col)

                    years_ds, ne_ds = coords_raw_from_df(df_ds, col)
                    ds_xy = coords_distance_transform(years_ds, ne_ds, args.distance_scale)

                    hd = hausdorff_symmetric(ds_xy, base_xy)
                    fr = frechet_dist(ds_xy, base_xy)
                    e = compute_sse(ne_base, ne_ds)
                    elog = compute_sse_log10(ne_base, ne_ds)

                    big_rows.append({
                        "Sample": sample,
                        "Coverage": float(cov),
                        "FNR_Value": round(fnr_val, 2),
                        "Hausdorff_Distance": hd,
                        "Frechet_Distance": fr,
                        "SSE": e,
                        "SSE_log10": elog,
                        "Distance_Scale": args.distance_scale
                    })

    big_df = pd.DataFrame(big_rows)
    if big_df.empty:
        raise SystemExit("No results. Check inputs and naming consistency.")

    big_tsv = os.path.join(outdir, "2_big_distances.tsv")
    big_df.to_csv(big_tsv, sep="\t", index=False)

    small_rows = []
    for (sample, cov), grp in big_df.groupby(["Sample", "Coverage"]):
        grp = grp.reset_index(drop=True)

        i_h = grp["Hausdorff_Distance"].astype(float).idxmin()
        best_h = grp.loc[i_h]
        i_f = grp["Frechet_Distance"].astype(float).idxmin()
        best_f = grp.loc[i_f]

        small_rows.append({
            "Sample": sample,
            "Coverage": float(cov),

            "Best_FNR_Hausdorff": float(best_h["FNR_Value"]),
            "Min_Hausdorff": float(best_h["Hausdorff_Distance"]),
            "SSE_at_Best_Hausdorff": float(best_h["SSE"]),
            "SSE_log10_at_Best_Hausdorff": float(best_h["SSE_log10"]) if pd.notna(best_h["SSE_log10"]) else np.nan,

            "Best_FNR_Frechet": float(best_f["FNR_Value"]),
            "Min_Frechet": float(best_f["Frechet_Distance"]),
            "SSE_at_Best_Frechet": float(best_f["SSE"]),
            "SSE_log10_at_Best_Frechet": float(best_f["SSE_log10"]) if pd.notna(best_f["SSE_log10"]) else np.nan,

            "Distance_Scale": args.distance_scale
        })

    small_df = pd.DataFrame(small_rows).sort_values(["Sample", "Coverage"]).reset_index(drop=True)
    small_tsv = os.path.join(outdir, "3_small_best.tsv")
    small_df.to_csv(small_tsv, sep="\t", index=False)

    curves_pdf = None
    curves_txt = None
    if not args.no_curves:
        out_prefix = os.path.join(outdir, "4_optimal_fnr")
        curves_pdf = os.path.join(outdir, "4_polynomials.pdf")
        curves_txt = os.path.join(outdir, "4_polynomials.txt")
        write_curves_and_plot(big_df, out_prefix, curves_pdf, curves_txt)

    if not args.no_trajectory_plots:
        plot_trajectories_by_sample(
            root=root,
            base_files=base_files,
            coverage_mapping=coverage_mapping,
            mu=args.mu, g=args.g, svalue=args.svalue,
            v_years=v_years,
            FNR_min=args.FNR_min, FNR_max=args.FNR_max,
            big_df=big_df,
            outdir_plots=outdir_plots
        )

    print("DONE")
    print("Wrote BIG:", big_tsv)
    print("Wrote SMALL:", small_tsv)
    if args.write_psmc_tsv:
        print("Wrote per-PSMC TSVs in:", outdir_psmc_tsv)
    if not args.no_curves:
        print("Curves outputs (if not skipped due to insufficient points):", curves_pdf, curves_txt)
    if not args.no_trajectory_plots:
        print("Wrote trajectory plots in:", outdir_plots)


def main():
    p = argparse.ArgumentParser(description="Unified pipeline (1->5) with distance-scale and optional baseline definition from table or input.")
    p.add_argument("--root", required=True, help="Root directory containing .psmc files")
    p.add_argument("--base_files", nargs="+", default=None, help="Exact baseline filenames (optional if --coverage_table defines baseline)")
    p.add_argument("--outdir", required=True, help="Output directory")

    p.add_argument("--mu", type=float, required=True)
    p.add_argument("--g", type=float, required=True)
    p.add_argument("--svalue", type=int, default=100)

    p.add_argument("--tmin", type=float, default=5e4)
    p.add_argument("--tmax", type=float, default=1.5e6)
    p.add_argument("--n_timepoints", type=int, default=60)

    p.add_argument("--FNR_min", type=float, default=0.0)
    p.add_argument("--FNR_max", type=float, default=0.99)

    p.add_argument("--distance_scale", choices=["psmc", "linear", "loglog"], default="psmc",
                   help="Scale used ONLY for Hausdorff/Frechet. Default: log10(time) and Ne*1e4 (PSMC-like).")

    p.add_argument("--coverage_table", default=None,
                   help="Optional TSV with columns: psmc_file and coverage. If provided, it may define baseline rows with 'baseline'.")

    p.add_argument("--no_curves", action="store_true", help="Do not generate polynomials/curve plots.")
    p.add_argument("--no_trajectory_plots", action="store_true", help="Do not generate trajectory plots (baseline vs corrected/uncorrected).")

    p.add_argument("--write_psmc_tsv", action="store_true",
                   help="Write step-1 per-PSMC TSVs (FNR trajectories sampled on the custom vector).")

    args = p.parse_args()
    run_pipeline(args)


if __name__ == "__main__":
    main()