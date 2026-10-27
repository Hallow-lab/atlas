import io, json, re, urllib.parse
import requests
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Lipinski, Crippen

st.set_page_config(page_title="MoleculeVoyager", page_icon="🧬", layout="wide", initial_sidebar_state="expanded")

# Molecules reflected in the supplied chemistry reports, plus useful starter examples.
STARTER_MOLECULES = [
    {"name":"Acetylsalicylic acid (Aspirin)","smiles":"CC(=O)Oc1ccccc1C(=O)O","category":"Pharmaceuticals"},
    {"name":"Acetic acid", "smiles":"CC(=O)O", "category":"Small molecules"},
    {"name":"Ethanol", "smiles":"CCO", "category":"Small molecules"},
    {"name":"Glucose", "smiles":"C(C1C(C(C(C(O1)O)O)O)O)O", "category":"Biomolecules"},
    {"name":"Caffeine", "smiles":"Cn1c(=O)c2c(ncn2C)n(C)c1=O", "category":"Pharmaceuticals"},
    {"name":"Paracetamol (Acetaminophen)", "smiles":"CC(=O)Nc1ccc(O)cc1", "category":"Pharmaceuticals"},
    {"name":"Ibuprofen", "smiles":"CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "category":"Pharmaceuticals"},
    {"name":"Adenine", "smiles":"C1=NC2=NC=NC(=C2N1)N", "category":"Biomolecules"},
    {"name":"Dopamine", "smiles":"NCCc1ccc(O)c(O)c1", "category":"Biomolecules"},
    {"name":"Lactic acid", "smiles":"C[C@H](O)C(=O)O", "category":"Small molecules"},
    {"name":"Benzene", "smiles":"c1ccccc1", "category":"Organic compounds"},
    {"name":"Urea", "smiles":"NC(=O)N", "category":"Small molecules"},
]

HTTP = requests.Session()
HTTP.headers.update({"User-Agent":"MoleculeVoyager/1.0 (educational cheminformatics app)"})

@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_search(term):
    # quote() accepts strings; quote_from_bytes errors are avoided by never passing encoded bytes.
    safe_term = urllib.parse.quote(str(term).strip(), safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{safe_term}/property/CanonicalSMILES,IsomericSMILES,MolecularFormula,MolecularWeight,IUPACName/JSON"
    response = HTTP.get(url, timeout=20)
    if response.status_code != 200:
        raise ValueError("No matching PubChem compound was found. Try a more specific name or use SMILES input.")
    items = response.json().get("PropertyTable", {}).get("Properties", [])
    if not items:
        raise ValueError("PubChem returned no compound records.")
    return items[0]

@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_image(smiles):
    # PubChem returns a 2D PNG. This avoids RDKit's optional native image dependencies.
    encoded = urllib.parse.quote(str(smiles), safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/PNG?image_size=large"
    response = HTTP.get(url, timeout=25)
    response.raise_for_status()
    if not response.content.startswith(b"\x89PNG"):
        raise ValueError("PubChem did not return a PNG image.")
    return response.content

def molecule_from_smiles(smiles):
    mol = Chem.MolFromSmiles(str(smiles).strip())
    if mol is None:
        raise ValueError("This SMILES string is not valid or is not supported by RDKit.")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol

def make_3d_block(mol, include_h=True):
    work = Chem.AddHs(Chem.Mol(mol)) if include_h else Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 2026
    status = AllChem.EmbedMolecule(work, params)
    if status != 0:
        raise ValueError("Could not generate a 3D conformer for this structure.")
    if AllChem.MMFFHasAllMoleculeParams(work):
        AllChem.MMFFOptimizeMolecule(work, maxIters=250)
    elif AllChem.UFFHasAllMoleculeParams(work):
        AllChem.UFFOptimizeMolecule(work, maxIters=250)
    return Chem.MolToMolBlock(work)

def get_properties(mol):
    return {
        "Molecular formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular weight (g/mol)": round(Descriptors.MolWt(mol), 4),
        "Exact molecular weight (Da)": round(Descriptors.ExactMolWt(mol), 4),
        "cLogP": round(Crippen.MolLogP(mol), 3),
        "Topological polar surface area (Å²)": round(Descriptors.TPSA(mol), 3),
        "H-bond donors": Lipinski.NumHDonors(mol),
        "H-bond acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable bonds": Lipinski.NumRotatableBonds(mol),
        "Ring count": rdMolDescriptors.CalcNumRings(mol),
        "Heavy atoms": mol.GetNumHeavyAtoms(),
        "Formal charge": Chem.GetFormalCharge(mol),
        "Fraction C sp³": round(rdMolDescriptors.CalcFractionCSP3(mol), 3),
    }

def get_atom_data(mol):
    """Return calculated structural properties for each atom."""
    atoms = []
    for atom in mol.GetAtoms():
        atoms.append({
            "Atom index": atom.GetIdx(),
            "Element": atom.GetSymbol(),
            "Atomic number": atom.GetAtomicNum(),
            "Isotope": atom.GetIsotope() or "Most abundant",
            "Formal charge": atom.GetFormalCharge(),
            "Aromatic": atom.GetIsAromatic(),
            "Hybridization": str(atom.GetHybridization()),
            "Degree": atom.GetDegree(),
            "Total valence": atom.GetTotalValence(),
            "Implicit H": atom.GetNumImplicitHs(),
            "Explicit H": atom.GetNumExplicitHs(),
            "Chiral tag": str(atom.GetChiralTag()),
        })
    return atoms


def get_bond_data(mol):
    """Return calculated structural properties for each bond."""
    bonds = []
    for bond in mol.GetBonds():
        bonds.append({
            "Bond index": bond.GetIdx(),
            "Atom 1": bond.GetBeginAtomIdx(),
            "Element 1": bond.GetBeginAtom().GetSymbol(),
            "Atom 2": bond.GetEndAtomIdx(),
            "Element 2": bond.GetEndAtom().GetSymbol(),
            "Bond type": str(bond.GetBondType()),
            "Bond order": bond.GetBondTypeAsDouble(),
            "Aromatic": bond.GetIsAromatic(),
            "Conjugated": bond.GetIsConjugated(),
            "In ring": bond.IsInRing(),
            "Stereo": str(bond.GetStereo()),
        })
    return bonds


def get_element_composition(mol):
    """Count atoms by element and calculate their percentage of all atoms."""
    from collections import Counter
    counts = Counter(atom.GetSymbol() for atom in mol.GetAtoms())
    total = sum(counts.values())
    return [
        {"Element": element, "Atom count": count,
         "Percentage (%)": round(100 * count / total, 2) if total else 0}
        for element, count in sorted(counts.items())
    ]


# Common element reference: general uses, plus a few biological/chemical roles.
# Uses are examples, not exhaustive lists or safety guidance.
ELEMENT_REFERENCE = {
    "H": {"name": "Hydrogen", "family": "Nonmetal", "uses": "Ammonia production, hydrogenation, fuel cells and water chemistry.", "role": "Present in water and most organic molecules."},
    "C": {"name": "Carbon", "family": "Nonmetal", "uses": "Organic chemicals, polymers, graphite, activated carbon and carbon fibre.", "role": "Forms the backbone of most organic molecules."},
    "N": {"name": "Nitrogen", "family": "Nonmetal", "uses": "Fertilizers, ammonia, pharmaceuticals and inert atmospheres.", "role": "Found in amino acids, proteins, nucleic acids and many drugs."},
    "O": {"name": "Oxygen", "family": "Nonmetal", "uses": "Medical oxygen, steelmaking, oxidation processes and water treatment.", "role": "Present in water and many functional groups."},
    "F": {"name": "Fluorine", "family": "Halogen", "uses": "Fluoropolymers and selected medicines and specialty chemicals.", "role": "Fluorine substitution can change a molecule's stability and properties."},
    "P": {"name": "Phosphorus", "family": "Nonmetal", "uses": "Fertilizers, detergents and specialty chemicals.", "role": "Important in ATP, DNA/RNA and phospholipids."},
    "S": {"name": "Sulfur", "family": "Nonmetal", "uses": "Sulfuric acid, rubber vulcanization and chemical synthesis.", "role": "Found in cysteine, methionine and many medicines."},
    "Cl": {"name": "Chlorine", "family": "Halogen", "uses": "Water disinfection, PVC manufacture and chemical synthesis.", "role": "Chloride is an important electrolyte; organochlorine compounds vary widely."},
    "Br": {"name": "Bromine", "family": "Halogen", "uses": "Specialty chemicals, flame-retardant chemistry and synthesis.", "role": "Present in some pharmaceuticals and bioactive compounds."},
    "I": {"name": "Iodine", "family": "Halogen", "uses": "Antiseptics, contrast agents and iodized salt.", "role": "Required for thyroid hormone synthesis in appropriate dietary amounts."},
    "Na": {"name": "Sodium", "family": "Alkali metal", "uses": "Sodium compounds in food, glass, chemical manufacture and industry.", "role": "Sodium ions are central to nerve impulses and fluid balance."},
    "K": {"name": "Potassium", "family": "Alkali metal", "uses": "Fertilizers, soaps and industrial chemicals.", "role": "Potassium ions are important for membrane potentials and muscle function."},
    "Ca": {"name": "Calcium", "family": "Alkaline earth metal", "uses": "Cement, lime, supplements and materials.", "role": "Calcium ions support bones, signaling and muscle contraction."},
    "Mg": {"name": "Magnesium", "family": "Alkaline earth metal", "uses": "Lightweight alloys, refractories and chemical reagents.", "role": "Magnesium is a cofactor in many enzymes and is central to chlorophyll."},
    "Fe": {"name": "Iron", "family": "Transition metal", "uses": "Steel production, catalysts and magnetic materials.", "role": "Iron is used in hemoglobin and many enzymes."},
    "Zn": {"name": "Zinc", "family": "Transition metal", "uses": "Galvanizing, brass, batteries and selected topical products.", "role": "Zinc is a cofactor in many enzymes and proteins."},
    "Cu": {"name": "Copper", "family": "Transition metal", "uses": "Electrical wiring, plumbing, alloys and catalysts.", "role": "Copper is a cofactor in several enzymes."},
    "Si": {"name": "Silicon", "family": "Metalloid", "uses": "Semiconductors, glass, ceramics and silicones.", "role": "Silicon compounds are widely used in materials and devices."},
    "B": {"name": "Boron", "family": "Metalloid", "uses": "Borosilicate glass, detergents and semiconductor materials.", "role": "Boron compounds have specialized chemical and material applications."},
    "Se": {"name": "Selenium", "family": "Nonmetal", "uses": "Glass, electronics and selected industrial applications.", "role": "Selenium is a trace nutrient in specific selenoproteins."},
    "Li": {"name": "Lithium", "family": "Alkali metal", "uses": "Rechargeable batteries, ceramics and specialty alloys.", "role": "Lithium salts have specific medical uses under professional supervision."},
    "Al": {"name": "Aluminium", "family": "Post-transition metal", "uses": "Transport, packaging, construction and alloys.", "role": "Aluminium compounds have applications in materials and some formulations."},
    "Ag": {"name": "Silver", "family": "Transition metal", "uses": "Electronics, mirrors, jewellery and selected antimicrobial materials.", "role": "Silver compounds have limited specialized medical/material uses."},
    "Au": {"name": "Gold", "family": "Transition metal", "uses": "Electronics, jewellery and specialized catalysts.", "role": "Gold compounds have some specialized research and medical applications."},
}

def element_reference(symbol):
    """Return reference information for an element symbol."""
    info = ELEMENT_REFERENCE.get(symbol, {})
    try:
        atom = Chem.Atom(symbol)
        return {
            "Element": symbol,
            "Element name": info.get("name", symbol),
            "Atomic number": atom.GetAtomicNum(),
            "Atomic mass (approx.)": round(atom.GetMass(), 4),
            "Chemical family": info.get("family", "See periodic table reference"),
            "Common uses": info.get("uses", "Uses vary by element and compound; consult a chemistry reference."),
            "Role in biology/chemistry": info.get("role", "Depends on its chemical form and context."),
        }
    except Exception:
        return {"Element": symbol, "Note": "Element reference unavailable."}


def viewer_html(molblock, style):
    # Cartoon is intended for protein/polymer backbones, not small molecules.
    styles = {
        "Ball and Stick": {"stick": {"radius": 0.16}, "sphere": {"scale": 0.27}},
        "Space Filling": {"sphere": {"scale": 0.85}},
        "Wedge and Dash": {"stick": {"radius": 0.12}, "sphere": {"scale": 0.22}},
        "Stick": {"stick": {"radius": 0.20}},
        "Wireframe": {"line": {"linewidth": 2}},
    }

    sty = styles.get(style, styles["Ball and Stick"])
    block_json = json.dumps(molblock)
    style_json = json.dumps(sty)

    return f"""
<div id="molecule-viewer"
     style="height:480px;width:100%;background:#ffffff;border-radius:16px;overflow:hidden;">
</div>
<div style="font:12px sans-serif;color:#64748b;padding:6px 2px;">
    Drag to rotate · Scroll to zoom · 3D coordinates are computationally generated
</div>

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
(function() {{
    const el = document.getElementById("molecule-viewer");

    if (!window.$3Dmol) {{
        el.innerHTML = '<div style="padding:24px;color:#334155;background:#ffffff;font-family:sans-serif;">3D viewer could not load. Check your internet connection and reload.</div>';
        return;
    }}

    const isDark = window.matchMedia &&
        window.matchMedia("(prefers-color-scheme: dark)").matches;

    const viewer = window.$3Dmol.createViewer(el, {{
        backgroundColor: isDark ? "#0b1220" : "#ffffff",
        antialias: true
    }});

    viewer.addModel({block_json}, "mol");
    viewer.setStyle({{}}, {style_json});
    viewer.zoomTo();
    viewer.render();

    window.addEventListener("resize", () => viewer.resize());
}})();
</script>
"""

def resolve_record(label, source, typed):
    if source == "SMILES":
        mol = molecule_from_smiles(typed)
        return {"name":"User-entered molecule", "smiles":Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), "mol":mol}
    if source == "Starter library":
        record = next((x for x in STARTER_MOLECULES if x["name"] == label), None)
        if not record: raise ValueError("Select a molecule from the library.")
        mol = molecule_from_smiles(record["smiles"])
        return {**record, "smiles":Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), "mol":mol}
    props = pubchem_search(typed)
    smiles = props.get("ConnectivitySMILES") or props.get("SMILES") or props.get("CanonicalSMILES") or props.get("ConnectivitySMILES")
    if not smiles: raise ValueError("PubChem record did not include a usable SMILES string.")
    mol = molecule_from_smiles(smiles)
    return {"name":typed, "smiles":Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True), "mol":mol, "cid":props.get("CID"), "iupac":props.get("IUPACName")}

# Visual identity
st.markdown("""<style>
@import url('https://fonts.googleapis.com/css2?family=DM+Sans:wght@400;500;600;700&family=Space+Grotesk:wght@500;600;700&display=swap');
html, body, [class*="css"] {font-family:'DM Sans',sans-serif;}
.stApp {background:linear-gradient(180deg,#f6f8fc 0%,#ffffff 38%);}
.block-container {max-width:1450px;padding-top:1.7rem;padding-bottom:3rem;}
.hero {
    padding: 26px 30px;
    border-radius: 22px;
    background: linear-gradient(
        120deg,
        #111b35 0%,
        #203e69 62%,
        #2b7b84 100%
    );
    color: #ffffff;
    margin-bottom: 22px;
}

.hero h1,
.hero h1 *,
.hero h2,
.hero h2 *,
.hero p,
.hero p *,
.hero span {
    color: #ffffff !important;
    -webkit-text-fill-color: #ffffff !important;
}
.hero h1 {
    font-family: 'Space Grotesk', sans-serif;
    font-size: 2.45rem;
    letter-spacing: -1.4px;
    margin: 0 0 7px 0;
    color: #ffffff !important;
    -webkit-text-fill-color: #ffffff !important;
    line-height: 1.2;
    font-weight: 700;
}
.hero p {
    color: #d7e6f7 !important;
    -webkit-text-fill-color: #d7e6f7 !important;
    margin: 0;
    font-size: 1rem;
    line-height: 1.6;
}
.hero .eyebrow {
    color: #77e0d1 !important;
    -webkit-text-fill-color: #77e0d1 !important;
    font-size: .72rem;
    font-weight: 700;
    letter-spacing: 2px;
    text-transform: uppercase;
    margin-bottom: 10px;
}

.panel {
    background: var(--secondary-background-color, #ffffff);
    color: var(--text-color, #17253f);
    border: 1px solid rgba(128, 128, 128, .25);
    border-radius: 18px;
    padding: 18px 20px;
    box-shadow: 0 8px 28px rgba(28, 47, 79, .08);
}

.panel * {
    color: inherit;
}

.small-label {
    font-size: .72rem;
    letter-spacing: 1px;
    text-transform: uppercase;
    color: var(--text-color, #718096);
    opacity: .75;
    font-weight: 700;
}

.big-value {
    font-family: 'Space Grotesk', sans-serif;
    font-size: 1.55rem;
    color: var(--text-color, #17253f);
    font-weight: 700;
}

.smiles {
    font-family: ui-monospace, SFMono-Regular, Consolas, monospace;
    background: var(--secondary-background-color, #f3f6fb);
    padding: 13px;
    border-radius: 10px;
    overflow-wrap: anywhere;
    color: var(--text-color, #233e63);
    border: 1px solid rgba(128, 128, 128, .2);
}

section[data-testid="stSidebar"] {
    background: var(--secondary-background-color, #101a2e);
}

section[data-testid="stSidebar"] * {
    color: var(--text-color, #edf4ff);
}
/* Consistent high-contrast text and surfaces throughout Streamlit widgets */
.stApp, .stApp p, .stApp label, .stApp li, .stApp span,
.stApp div[data-testid="stMarkdownContainer"], .stApp [data-testid="stWidgetLabel"] {
    color: #17253f;
}
.stApp { background: #f5f7fb; }
div[data-testid="stMetric"] {
    background: #ffffff;
    border: 1px solid #dce4ef;
    padding: 12px;
    border-radius: 12px;
}
div[data-testid="stMetric"] label,
div[data-testid="stMetric"] [data-testid="stMetricValue"],
div[data-testid="stMetric"] [data-testid="stMetricDelta"] {
    color: #17253f !important;
}
.stTabs [data-baseweb="tab"] { color: #334155 !important; }
.stTabs [aria-selected="true"] { color: #0f766e !important; }
[data-testid="stDataFrame"], [data-testid="stTable"] { color: #17253f; }
[data-testid="stAlert"] p { color: #17253f !important; }
section[data-testid="stSidebar"] {
    background: #101a2e !important;
}
section[data-testid="stSidebar"] * {
    color: #edf4ff !important;
}
section[data-testid="stSidebar"] input,
section[data-testid="stSidebar"] textarea,
section[data-testid="stSidebar"] [data-baseweb="select"] > div {
    color: #17253f !important;
    background-color: #ffffff !important;
}
.hero, .hero * { color: #ffffff !important; -webkit-text-fill-color: #ffffff !important; }
.hero .eyebrow { color: #77e0d1 !important; -webkit-text-fill-color: #77e0d1 !important; }
.hero p { color: #d7e6f7 !important; -webkit-text-fill-color: #d7e6f7 !important; }
</style>""", unsafe_allow_html=True)

st.markdown('''<div class="hero"><div class="eyebrow">CHEMISTRY · STRUCTURE · DISCOVERY</div><h1>🧬 MoleculeVoyager</h1><p>Explore molecular structures, 3D geometry, canonical SMILES, and physicochemical properties in one workspace.</p></div>''', unsafe_allow_html=True)

with st.sidebar:
    st.markdown("## Explore")
    source = st.radio("Molecule source", ["Starter library", "PubChem name search", "SMILES"], index=0)
    if source == "Starter library":
        categories = sorted(set(x["category"] for x in STARTER_MOLECULES))
        category = st.selectbox("Collection", ["All collections"] + categories)
        choices = [x for x in STARTER_MOLECULES if category == "All collections" or x["category"] == category]
        search = st.text_input("Filter library", placeholder="e.g. aspirin, caffeine")
        if search.strip(): choices = [x for x in choices if search.lower() in x["name"].lower() or search.lower() in x["smiles"].lower()]
        selected_name = st.selectbox("Open molecule", [x["name"] for x in choices]) if choices else None
        typed = selected_name or ""
    elif source == "PubChem name search":
        typed = st.text_input("Compound name", value="aspirin", placeholder="e.g. caffeine")
        selected_name = typed
    else:
        typed = st.text_area("SMILES string", value="CCO", height=100, help="SMILES is parsed locally by RDKit.")
        selected_name = "User-entered molecule"
    st.markdown("---")
    style = st.selectbox("3D representation", ["Ball and Stick", "Space Filling", "Wedge and Dash", "Stick", "Wireframe", "Cartoon"])
    include_h = st.checkbox("Include hydrogens in 3D", value=True)
    open_button = st.button("Open molecule", type="primary", use_container_width=True)
    st.caption("PubChem name search and 2D images require internet access. The starter library works without name lookup.")

if "atlas_record" not in st.session_state:
    st.session_state.atlas_record = None
if open_button:
    try:
        with st.spinner("Loading structure and calculating properties…"):
            record = resolve_record(selected_name, source, typed)
            record["molblock"] = make_3d_block(record["mol"], include_h)
            st.session_state.atlas_record = record
            st.session_state.atlas_source = source
            st.session_state.atlas_error = None
    except Exception as exc:
        st.session_state.atlas_error = str(exc)
if st.session_state.get("atlas_error"):
    st.error(st.session_state.atlas_error)

record = st.session_state.atlas_record
if record:
    mol = record["mol"]
    props = get_properties(mol)
    st.markdown(f"### {record['name']}")
    st.caption("Molecule detail page · Select a different entry in the sidebar and choose Open molecule to switch compounds.")
    metric_cols = st.columns(4)
    for col, (label, value) in zip(metric_cols, [("FORMULA",props["Molecular formula"]),("MOLECULAR WEIGHT",f"{props['Molecular weight (g/mol)']:.3f} g/mol"),("HEAVY ATOMS",str(props["Heavy atoms"])),("RINGS",str(props["Ring count"]))]):
        with col:
            st.markdown(f'<div class="panel"><div class="small-label">{label}</div><div class="big-value">{value}</div></div>', unsafe_allow_html=True)
    st.write("")
    tab2d, tab3d, tabprops, tabatoms, tabsmiles = st.tabs([
        "◫  2D Structure", "◉  3D View", "▤  Properties",
        "⚛  Atom Explorer", "⌘  SMILES & Export"
    ])
    with tab2d:
        left, right = st.columns([1.2,1])
        with left:
            st.markdown('<div class="panel">', unsafe_allow_html=True)
            st.markdown("#### 2D molecular structure")
            try:
                png = pubchem_image(record["smiles"])
                st.image(png, use_container_width=True)
            except Exception:
                st.info("PubChem 2D image is unavailable. The molecule can still be explored in 3D and its properties remain available.")
            st.markdown('</div>', unsafe_allow_html=True)
        with right:
            st.markdown('<div class="panel">', unsafe_allow_html=True)
            st.markdown("#### At a glance")
            st.write(f"**Molecular formula:** {props['Molecular formula']}")
            st.write(f"**Formal charge:** {props['Formal charge']}")
            st.write(f"**Hydrogen-bond donors:** {props['H-bond donors']}")
            st.write(f"**Hydrogen-bond acceptors:** {props['H-bond acceptors']}")
            st.write(f"**Rotatable bonds:** {props['Rotatable bonds']}")
            st.caption("Physicochemical descriptors are computed from the parsed structure using RDKit.")
            st.markdown('</div>', unsafe_allow_html=True)
    with tab3d:
        if style == "Cartoon":
            st.warning("Cartoon representation is intended for polymer/protein backbones and is not applicable to ordinary small molecules. Showing Ball and Stick instead.")
            viewer_style = "Ball and Stick"
        else: viewer_style = style
        components.html(viewer_html(record["molblock"], viewer_style), height=520, scrolling=False)
        st.caption(f"Displayed as **{viewer_style}**. Generated coordinates are a computational conformer, not an experimentally measured structure.")
    with tabprops:
        st.markdown("#### Physicochemical descriptors")
        rows = [{"Property":k,"Value":v} for k,v in props.items()]
        st.dataframe(rows, use_container_width=True, hide_index=True)
        st.markdown("#### Interpretation notes")
        st.write("- **cLogP** estimates lipophilicity; it is not a measured partition coefficient.")
        st.write("- **TPSA** estimates topological polar surface area and can support early permeability comparisons.")
        st.write("- Hydrogen-bond donor/acceptor counts and rotatable bonds are structural descriptors, not direct clinical outcomes.")
    with tabatoms:
        st.subheader("Individual Atom Explorer")
        st.caption("Inspect atom-level properties and connected bonds. Atom indices start at 0.")
        atom_data = get_atom_data(mol)
        bond_data = get_bond_data(mol)
        element_data = get_element_composition(mol)

        if atom_data:
            atom_labels = [
                f"Atom {a['Atom index']} — {a['Element']} (Z={a['Atomic number']})"
                for a in atom_data
            ]
            selected_label = st.selectbox(
                "Select an individual atom", atom_labels,
                key="atom_explorer_selection"
            )
            selected_index = atom_labels.index(selected_label)
            selected_atom = atom_data[selected_index]

            st.markdown("#### Selected atom")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Element", selected_atom["Element"])
            c2.metric("Atomic number", selected_atom["Atomic number"])
            c3.metric("Formal charge", selected_atom["Formal charge"])
            c4.metric("Hybridization", selected_atom["Hybridization"])

            st.markdown("#### Complete atom properties")
            st.dataframe(
                [{"Property": k, "Value": v} for k, v in selected_atom.items()],
                use_container_width=True, hide_index=True
            )

            connected_bonds = [
                b for b in bond_data
                if b["Atom 1"] == selected_index or b["Atom 2"] == selected_index
            ]
            st.markdown("#### Bonds connected to selected atom")
            if connected_bonds:
                st.dataframe(connected_bonds, use_container_width=True, hide_index=True)
            else:
                st.info("This atom has no bonds.")

            st.markdown("#### Element composition")
            st.dataframe(element_data, use_container_width=True, hide_index=True)
            st.markdown("#### All atoms")
            st.dataframe(atom_data, use_container_width=True, hide_index=True)
            import pandas as pd
            st.download_button(
                "Download atom properties CSV",
                data=pd.DataFrame(atom_data).to_csv(index=False),
                file_name="atom_properties.csv", mime="text/csv"
            )
            st.markdown("#### All bonds")
            if bond_data:
                st.dataframe(bond_data, use_container_width=True, hide_index=True)
                st.download_button(
                    "Download bond properties CSV",
                    data=pd.DataFrame(bond_data).to_csv(index=False),
                    file_name="bond_properties.csv", mime="text/csv"
                )
            else:
                st.info("This molecule has no bonds.")
        else:
            st.info("No atoms are available for this molecule.")
    with tabuses:
        st.subheader("Elements, chemistry and common uses")
        st.caption("Element uses below are general examples. An element's behaviour depends on its chemical form, concentration and context.")
        composition = get_element_composition(mol)
        if composition:
            selected_symbol = st.selectbox(
                "Choose an element present in this molecule",
                [item["Element"] for item in composition],
                key="element_uses_selection"
            )
            info = element_reference(selected_symbol)
            st.markdown(f"### {info['Element name']} ({selected_symbol})")
            c1, c2, c3 = st.columns(3)
            c1.metric("Atomic number", info.get("Atomic number", "—"))
            c2.metric("Approx. atomic mass", info.get("Atomic mass (approx.)", "—"))
            c3.metric("Atoms in molecule", next(x["Atom count"] for x in composition if x["Element"] == selected_symbol))
            st.write(f"**Chemical family:** {info.get('Chemical family', '—')}")
            st.write(f"**Common uses:** {info.get('Common uses', '—')}")
            st.write(f"**Biological/chemical role:** {info.get('Role in biology/chemistry', '—')}")
            st.markdown("#### Composition of this molecule")
            st.dataframe(composition, use_container_width=True, hide_index=True)
            st.info("The element reference describes general uses of the element or its compounds; it does not mean the whole molecule has those same uses.")
            import pandas as pd
            st.download_button(
                "Download element composition CSV",
                data=pd.DataFrame(composition).to_csv(index=False),
                file_name="element_composition.csv", mime="text/csv"
            )
        else:
            st.info("No element composition is available for this molecule.")

    with tabsmiles:
        st.markdown("#### Canonical isomeric SMILES")
        st.markdown(f'<div class="smiles">{record["smiles"]}</div>', unsafe_allow_html=True)
        st.download_button("Download molecule as SDF", data=record["molblock"], file_name=re.sub(r"[^A-Za-z0-9_-]+","_",record["name"]).strip("_")+".sdf", mime="chemical/x-mdl-sdfile")
        st.download_button("Download properties as CSV", data=("Property,Value\n"+"\n".join(f'"{k}","{v}"' for k,v in props.items())), file_name="molecule_properties.csv", mime="text/csv")
        st.caption("Canonical SMILES is generated from the parsed structure. Stereochemistry is retained where encoded in the input.")
else:
    st.markdown("### Start exploring")
    a,b,c = st.columns(3)
    with a:
        st.markdown('<div class="panel"><div class="small-label">01 · Browse</div><h3>Open a molecule</h3><p>Choose a starter compound or search PubChem by compound name.</p></div>', unsafe_allow_html=True)
    with b:
        st.markdown('<div class="panel"><div class="small-label">02 · Visualize</div><h3>Compare structures</h3><p>Switch between 2D connectivity and interactive 3D representations.</p></div>', unsafe_allow_html=True)
    with c:
        st.markdown('<div class="panel"><div class="small-label">03 · Analyze</div><h3>Inspect properties</h3><p>Review molecular descriptors and export structure or property data.</p></div>', unsafe_allow_html=True)
   

st.markdown("---")
st.caption("MoleculeVoyager · Educational cheminformatics tool · RDKit + PubChem PUG REST + 3Dmol.js")
