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
   

st.markdown("---")
st.caption("MoleculeVoyager · Educational cheminformatics tool · RDKit + PubChem PUG REST + 3Dmol.js")
