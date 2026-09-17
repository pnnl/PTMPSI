#!/usr/bin/env python3

import argparse
from pathlib import Path
import warnings

import numpy as np
import pandas as pd

from sklearn.model_selection import StratifiedKFold, cross_validate, cross_val_predict
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, classification_report


def parse_args():
	ap = argparse.ArgumentParser(
		description="Combo1 v2: classify Viruses vs Bacteria from CP12 HMM per-position features"
	)
	ap.add_argument("--accession_long", required=True, help="*_accession_per_column_step2v6.txt")
	ap.add_argument("--out_dir", required=True, help="Output directory")
	ap.add_argument(
		"--feature_set",
		choices=["chem", "aa", "both"],
		default="chem",
		help="Use ChemClass, AA, or both as features"
	)
	ap.add_argument(
		"--use_gaps",
		action="store_true",
		help="If set, keep gap states as a category, otherwise drop gap positions for that accession"
	)
	ap.add_argument("--kfold", type=int, default=5, help="Number of CV folds")
	ap.add_argument("--seed", type=int, default=1, help="Random seed")
	ap.add_argument(
		"--silence_warnings",
		action="store_true",
		help="If set, silence DeprecationWarnings for cleaner logs"
	)
	return ap.parse_args()


def add_domain_and_label(df: pd.DataFrame) -> pd.DataFrame:
	df = df.copy()

	df["Domain"] = df["LineageCleaned"].astype(str).str.split(":", n=1).str[0]
	df = df[df["Domain"].isin(["Viruses", "Bacteria"])].copy()

	df["y"] = df["Domain"].map({"Bacteria": 0, "Viruses": 1}).astype(int)

	return df


def pivot_features(df: pd.DataFrame, feature_cols, use_gaps: bool = False) -> pd.DataFrame:
	tmp = df.copy()

	if not use_gaps:
		tmp = tmp[tmp["IsGap"] == 0].copy()

	wide_parts = []

	for feat in feature_cols:
		w = tmp.pivot_table(
			index="Accession",
			columns="HMM_col",
			values=feat,
			aggfunc="first"
		)
		w.columns = [f"{feat}_{int(c)}" for c in w.columns]
		wide_parts.append(w)

	X = pd.concat(wide_parts, axis=1)
	X = X.fillna("MISSING")

	return X


def get_ohe_feature_names(enc: OneHotEncoder, input_features):
	if hasattr(enc, "get_feature_names_out"):
		return enc.get_feature_names_out(input_features)
	return enc.get_feature_names(input_features)


def main():
	args = parse_args()

	if args.silence_warnings:
		warnings.filterwarnings("ignore", category=DeprecationWarning)

	out_dir = Path(args.out_dir).resolve()
	out_dir.mkdir(parents=True, exist_ok=True)

	df = pd.read_csv(args.accession_long, sep="\t", dtype=str)

	needed = ["Accession", "HMM_col", "AA", "ChemClass", "IsGap", "LineageCleaned"]
	for c in needed:
		if c not in df.columns:
			raise ValueError(f"Missing required column in accession_long: {c}")

	df["HMM_col"] = df["HMM_col"].astype(int)
	df["IsGap"] = df["IsGap"].astype(int)

	df2 = add_domain_and_label(df)

	acc_labels = (
		df2[["Accession", "y", "Domain", "LineageCleaned"]]
		.drop_duplicates("Accession")
		.set_index("Accession")
	)

	if args.feature_set == "chem":
		feature_cols = ["ChemClass"]
	elif args.feature_set == "aa":
		feature_cols = ["AA"]
	else:
		feature_cols = ["ChemClass", "AA"]

	X = pivot_features(df2, feature_cols=feature_cols, use_gaps=args.use_gaps)

	X = X.loc[X.index.intersection(acc_labels.index)].copy()
	acc_labels = acc_labels.loc[X.index].copy()
	y = acc_labels["y"].astype(int).values

	counts = pd.Series(y).value_counts().to_dict()
	print("Class balance (0=Bacteria, 1=Viruses):", counts)

	cat_cols = list(X.columns)

	pre = ColumnTransformer(
		transformers=[
			("cat", OneHotEncoder(handle_unknown="ignore"), cat_cols)
		],
		remainder="drop"
	)

	clf = LogisticRegression(
		C=1.0,
		solver="liblinear",
		class_weight="balanced",
		max_iter=2000,
		random_state=args.seed
	)

	pipe = Pipeline([
		("pre", pre),
		("clf", clf)
	])

	cv = StratifiedKFold(
		n_splits=args.kfold,
		shuffle=True,
		random_state=args.seed
	)

	scoring = {
		"auc": "roc_auc",
		"ap": "average_precision",
		"acc": "accuracy",
		"f1": "f1"
	}

	res = cross_validate(
		pipe,
		X,
		y,
		cv=cv,
		scoring=scoring,
		return_train_score=False
	)

	summary = {k: float(np.mean(v)) for k, v in res.items() if k.startswith("test_")}
	summary_sd = {k + "_sd": float(np.std(v)) for k, v in res.items() if k.startswith("test_")}

	print("CV mean scores:", summary)
	print("CV sd scores:", summary_sd)

	# Out-of-fold prediction probabilities for ROC and PR curves
	y_prob_cv = cross_val_predict(
		pipe,
		X,
		y,
		cv=cv,
		method="predict_proba"
	)[:, 1]

	y_pred_cv = (y_prob_cv >= 0.5).astype(int)

	cv_pred_df = pd.DataFrame({
		"Accession": X.index,
		"y_true": y,
		"true_label": np.where(y == 1, "Viruses", "Bacteria"),
		"prob_viral": y_prob_cv,
		"y_pred": y_pred_cv,
		"pred_label": np.where(y_pred_cv == 1, "Viruses", "Bacteria"),
		"Domain": acc_labels["Domain"].values,
		"LineageCleaned": acc_labels["LineageCleaned"].values
	})

	cv_pred_path = out_dir / "combo1_cv_predictions.tsv"
	cv_pred_df.to_csv(cv_pred_path, sep="\t", index=False)
	print(f"Wrote: {cv_pred_path}")

	# CV confusion matrix
	cm_cv = confusion_matrix(y, y_pred_cv)
	rep_cv = classification_report(
		y,
		y_pred_cv,
		target_names=["Bacteria", "Viruses"],
		digits=3
	)

	cm_cv_path = out_dir / "combo1_confusion_matrix_cv.tsv"
	pd.DataFrame(
		cm_cv,
		index=["Bacteria_true", "Viruses_true"],
		columns=["Bacteria_pred", "Viruses_pred"]
	).to_csv(cm_cv_path, sep="\t")
	print(f"Wrote: {cm_cv_path}")

	rep_cv_path = out_dir / "combo1_classification_report_cv.txt"
	with rep_cv_path.open("w", encoding="utf-8") as f:
		f.write(rep_cv + "\n")
	print(f"Wrote: {rep_cv_path}")

	# Fit on full dataset for coefficients and full-fit diagnostics
	pipe.fit(X, y)

	enc = pipe.named_steps["pre"].named_transformers_["cat"]
	feat_names = get_ohe_feature_names(enc, cat_cols)

	coefs = pipe.named_steps["clf"].coef_.ravel()

	coef_df = (
		pd.DataFrame({
			"feature": feat_names,
			"coef": coefs,
			"abs_coef": np.abs(coefs)
		})
		.sort_values("abs_coef", ascending=False)
	)

	coef_path = out_dir / "combo1_logreg_coefficients.tsv"
	coef_df.to_csv(coef_path, sep="\t", index=False)
	print(f"Wrote: {coef_path}")

	yhat_full = pipe.predict(X)

	cm_full = confusion_matrix(y, yhat_full)
	rep_full = classification_report(
		y,
		yhat_full,
		target_names=["Bacteria", "Viruses"],
		digits=3
	)

	cm_full_path = out_dir / "combo1_confusion_matrix_fullfit.tsv"
	pd.DataFrame(
		cm_full,
		index=["Bacteria_true", "Viruses_true"],
		columns=["Bacteria_pred", "Viruses_pred"]
	).to_csv(cm_full_path, sep="\t")
	print(f"Wrote: {cm_full_path}")

	rep_full_path = out_dir / "combo1_classification_report_fullfit.txt"
	with rep_full_path.open("w", encoding="utf-8") as f:
		f.write(rep_full + "\n")
	print(f"Wrote: {rep_full_path}")

	run_path = out_dir / "combo1_run_summary.tsv"

	out = {
		"n_accessions": int(X.shape[0]),
		"n_features_raw": int(X.shape[1]),
		"feature_set": args.feature_set,
		"use_gaps": int(args.use_gaps),
		"kfold": args.kfold,
		"seed": args.seed,
		"test_auc_mean": summary.get("test_auc", np.nan),
		"test_auc_sd": summary_sd.get("test_auc_sd", np.nan),
		"test_ap_mean": summary.get("test_ap", np.nan),
		"test_ap_sd": summary_sd.get("test_ap_sd", np.nan),
		"test_acc_mean": summary.get("test_acc", np.nan),
		"test_acc_sd": summary_sd.get("test_acc_sd", np.nan),
		"test_f1_mean": summary.get("test_f1", np.nan),
		"test_f1_sd": summary_sd.get("test_f1_sd", np.nan),
		"cv_predictions_file": str(cv_pred_path),
		"cv_confusion_matrix_file": str(cm_cv_path),
		"fullfit_confusion_matrix_file": str(cm_full_path),
		"coefficients_file": str(coef_path)
	}

	pd.DataFrame([out]).to_csv(run_path, sep="\t", index=False)
	print(f"Wrote: {run_path}")


if __name__ == "__main__":
	main()