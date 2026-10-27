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
