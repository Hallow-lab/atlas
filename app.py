import io, json, re, urllib.parse
import requests
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, Lipinski, Crippen

st.set_page_config(page_title="Molecule Atlas", page_icon="⚛️", layout="wide", initial_sidebar_state="expanded")

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
HTTP.headers.update({"User-Agent":"MoleculeAtlas/1.0 (educational cheminformatics app)"})

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

def viewer_html(molblock, style):
    # Style definitions use 3Dmol.js representations. Cartoon is not applicable to ordinary small molecules.
    styles = {
        "Ball and Stick": "{stick:{radius:0.16}, sphere:{scale:0.27}}",
        "Space Filling": "{sphere:{scale:0.85}}",
        "Wedge and Dash": "{stick:{radius:0.12}, sphere:{scale:0.22}}",
        "Stick": "{stick:{radius:0.20}}",
        "Wireframe": "{line:{linewidth:2}}",
    }
    sty = styles.get(style, styles["Ball and Stick"])
    # JSON encoding keeps molblock safely quoted inside JavaScript.
    block_json = json.dumps(molblock)
    return f'''<div id="molecule-viewer" style="height:480px;width:100%;background:#0b1220;border-radius:16px;overflow:hidden"></div>
<div style="font:12px sans-serif;color:#64748b;padding:6px 2px">Drag to rotate · Scroll to zoom · 3D coordinates are computationally generated</div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
(function(){{
 const el=document.getElementById('molecule-viewer');
 if(!window.$3Dmol){{el.innerHTML='<div style="padding:24px;color:#fff">3D viewer could not load. Check your internet connection and reload.</div>';return;}}
 const viewer=$3Dmol.createViewer(el,{{backgroundColor:'#0b1220', antialias:true}});
 viewer.addModel({block_json},'mol');
 viewer.setStyle({{}},{sty});
 viewer.zoomTo(); viewer.render();
 window.addEventListener('resize',()=>viewer.resize());
}})();
</script>'''

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
.hero {padding:26px 30px;border-radius:22px;background:linear-gradient(120deg,#111b35 0%,#203e69 62%,#2b7b84 100%);color:white;margin-bottom:22px;}
.hero h1 {font-family:'Space Grotesk',sans-serif;font-size:2.45rem;letter-spacing:-1.4px;margin:0 0 7px 0;color:white;}
.hero p {color:#d7e6f7;margin:0;font-size:1rem;}
.eyebrow {font-size:.72rem;font-weight:700;letter-spacing:2px;text-transform:uppercase;color:#77e0d1;margin-bottom:10px;}
.panel {background:white;border:1px solid #e5eaf2;border-radius:18px;padding:18px 20px;box-shadow:0 8px 28px rgba(28,47,79,.045);}
.small-label {font-size:.72rem;letter-spacing:1px;text-transform:uppercase;color:#718096;font-weight:700;}
.big-value {font-family:'Space Grotesk',sans-serif;font-size:1.55rem;color:#17253f;font-weight:700;}
.smiles {font-family:ui-monospace,SFMono-Regular,Consolas,monospace;background:#f3f6fb;padding:13px;border-radius:10px;overflow-wrap:anywhere;color:#233e63;}
section[data-testid="stSidebar"] {background:#101a2e;}
section[data-testid="stSidebar"] * {color:#edf4ff;}
div[data-testid="stMetric"] {background:#f8fafc;border:1px solid #e8edf4;padding:12px;border-radius:12px;}
</style>""", unsafe_allow_html=True)

st.markdown('''<div class="hero"><div class="eyebrow">CHEMISTRY · STRUCTURE · DISCOVERY</div><h1>⚛ Molecule Atlas</h1><p>Explore molecular structures, 3D geometry, canonical SMILES, and physicochemical properties in one workspace.</p></div>''', unsafe_allow_html=True)

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
    tab2d, tab3d, tabprops, tabsmiles = st.tabs(["◫  2D Structure", "◉  3D View", "▤  Properties", "⌘  SMILES & Export"])
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
    st.info("The reference PDF contains example reports for aspirin, acetic acid, and ethanol; these are included in the starter library. No finite app can preload every known molecule, so PubChem search provides access to a much broader compound collection, subject to its availability and internet access.")

st.markdown("---")
st.caption("Molecule Atlas · Educational cheminformatics tool · RDKit + PubChem PUG REST + 3Dmol.js")
