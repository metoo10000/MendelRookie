from pathlib import Path

path = Path("nutrition_stage/prepare_c_n08.py")
text = path.read_text(encoding="utf-8")
old = '''    merged = pieces[0]
    for piece in pieces[1:]:
        merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")
'''
new = '''    merged = pieces[0]
    for piece in pieces[1:]:
        if not piece.empty and "SEQN" in piece.columns:
            merged = merged.merge(piece.drop_duplicates("SEQN"), on="SEQN", how="left")
'''
if old not in text:
    raise SystemExit("Frozen merge target not found; refusing broader modification")
path.write_text(text.replace(old, new, 1), encoding="utf-8")
print("Applied optional-module merge guard; no cohort, exposure, outcome, covariate, or gate changed")
