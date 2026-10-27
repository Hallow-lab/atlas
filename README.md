# Molecule Atlas

A redesigned Streamlit molecule explorer with a molecule detail view containing:

- 2D structure image (PubChem PUG REST; avoids RDKit's native 2D drawing dependencies)
- Interactive 3D view powered by 3Dmol.js
- Ball and Stick, Space Filling, Wedge and Dash, Stick, and Wireframe representations
- Clear fallback when Cartoon is selected for a small molecule
- Molecular formula, molecular weight, exact mass, cLogP, TPSA, H-bond donors/acceptors, rotatable bonds, ring count, heavy atoms, charge, and fraction C sp3
- Canonical isomeric SMILES
- SDF and property CSV downloads
- Starter collection, PubChem name search, and direct SMILES input

## Run locally

```bash
python -m pip install -r requirements.txt
streamlit run app.py
```

## Deploy to Streamlit Community Cloud

1. Upload `app.py`, `requirements.txt`, and `packages.txt` to a GitHub repository.
2. Create a new Streamlit app using `app.py` as the entry point.
3. The app requires internet access for PubChem name lookups, 2D PNGs, and the 3Dmol.js CDN. Starter library records and RDKit calculations are local.

## Notes

- The app does not claim to preload every molecule ever described. PubChem name search gives access to a much larger public compound database.
- 3D conformers are computational estimates, not experimental coordinates.
- If PubChem does not resolve a name, try a more specific name or enter a valid SMILES string.
- The `quote_from_bytes()` issue is avoided by URL-encoding string values with `urllib.parse.quote()`.
