from flask import Flask, render_template
import json
import os
import pandas as pd

app = Flask(__name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
SNP_FILE = os.path.join(BASE_DIR, "snpdict.json")
KETO_FILE = os.path.join(BASE_DIR, "ketoDict.json")

COMPLEMENT = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}


def load_json_file(path):
    with open(path, "r", encoding="utf-8") as f:
        return json.load(f)


def normalize_genotype(genotype):
    if genotype is None:
        return ""

    cleaned = str(genotype).upper().replace("/", "").replace(" ", "")
    cleaned = "".join(ch for ch in cleaned if ch in {"A", "C", "G", "T"})

    if not cleaned:
        return ""

    return "".join(sorted(cleaned))


def flip_genotype(genotype):
    normalized = normalize_genotype(genotype)
    if not normalized:
        return ""

    flipped = "".join(COMPLEMENT.get(base, base) for base in normalized)
    return "".join(sorted(flipped))


def normalize_reference_options(reference_genotype):
    if reference_genotype is None:
        return set()

    raw = str(reference_genotype).upper()
    parts = [part.strip() for part in raw.split("/")]

    options = set()
    for part in parts:
        norm = normalize_genotype(part)
        if norm:
            options.add(norm)

    if not options:
        norm = normalize_genotype(raw)
        if norm:
            options.add(norm)

    return options


def evaluate_match(user_genotype, reference_genotype, orientation):
    user_norm = normalize_genotype(user_genotype)
    flipped_norm = flip_genotype(user_genotype)
    ref_options = normalize_reference_options(reference_genotype)

    is_minus = str(orientation).strip().lower() == "minus"

    direct_match = user_norm in ref_options if user_norm else False
    flipped_match = flipped_norm in ref_options if is_minus and flipped_norm else False

    if direct_match:
        match_status = "Direct match"
    elif flipped_match:
        match_status = "Match after strand flip"
    else:
        match_status = "No match"

    return {
        "user_genotype_normalized": user_norm,
        "flipped_genotype": flipped_norm if is_minus else "",
        "reference_options": ", ".join(sorted(ref_options)),
        "match_status": match_status
    }


def join_snp_data(snp_dict, keto_list):
    df_snp = pd.DataFrame(snp_dict.items(), columns=["rsid", "user_genotype"])
    df_keto = pd.DataFrame(keto_list)

    if "rsid" not in df_keto.columns:
        raise ValueError("ketoDict.json must contain an 'rsid' field.")

    df_merged = df_snp.merge(df_keto, on="rsid", how="inner")

    rows = []
    merged_dict = {}

    for row in df_merged.to_dict(orient="records"):
        comparison = evaluate_match(
            user_genotype=row.get("user_genotype"),
            reference_genotype=row.get("genotype"),
            orientation=row.get("orientation", "plus")
        )

        enriched = {
            **row,
            **comparison
        }
        rows.append(enriched)

        rsid = row["rsid"]
        merged_dict[rsid] = {
            "user_genotype": row.get("user_genotype"),
            "user_genotype_normalized": comparison["user_genotype_normalized"],
            "gene": row.get("gene"),
            "category": row.get("category"),
            "orientation": row.get("orientation"),
            "reference_genotype": row.get("genotype"),
            "reference_options": comparison["reference_options"],
            "flipped_genotype": comparison["flipped_genotype"],
            "match_status": comparison["match_status"],
            "evidence_level": row.get("evidence_level"),
            "actionability": row.get("actionability"),
            "nutrition_relevance": row.get("nutrition_relevance"),
            "keto_effect": row.get("keto_effect"),
            "description": row.get("description"),
            "source": row.get("source"),
        }

    categories = sorted(
        {row.get("category") for row in rows if row.get("category")}
    )
    evidence_levels = sorted(
        {row.get("evidence_level") for row in rows if row.get("evidence_level")}
    )
    actionability_levels = sorted(
        {row.get("actionability") for row in rows if row.get("actionability")}
    )

    return (
        merged_dict,
        rows,
        len(snp_dict),
        len(keto_list),
        len(rows),
        categories,
        evidence_levels,
        actionability_levels
    )


@app.route("/")
def index():
    try:
        snp_dict = load_json_file(SNP_FILE)
        keto_list = load_json_file(KETO_FILE)

        (
            merged_dict,
            rows,
            snp_count,
            keto_count,
            common_count,
            categories,
            evidence_levels,
            actionability_levels
        ) = join_snp_data(snp_dict, keto_list)

        return render_template(
            "index.html",
            rows=rows,
            merged_dict=merged_dict,
            snp_count=snp_count,
            keto_count=keto_count,
            common_count=common_count,
            categories=categories,
            evidence_levels=evidence_levels,
            actionability_levels=actionability_levels,
            error=None,
        )

    except Exception as e:
        return render_template(
            "index.html",
            rows=[],
            merged_dict={},
            snp_count=0,
            keto_count=0,
            common_count=0,
            categories=[],
            evidence_levels=[],
            actionability_levels=[],
            error=str(e),
        )


if __name__ == "__main__":
    app.run(debug=True)