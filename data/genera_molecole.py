import requests
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
import os
import cairosvg
from PIL import Image
import warnings

# Suppress InsecureRequestWarning if verify=False is used for debugging
requests.packages.urllib3.disable_warnings(requests.packages.urllib3.exceptions.InsecureRequestWarning)

print("--- Script execution started ---")

def get_smiles_from_cas(cas_number):
    """
    Retrieves the canonical SMILES string of a molecule from PubChem by its CAS number
    using a two-step process: CAS to CID via xref/rn, then CID to SMILES.
    """
    cas_number = cas_number.strip()
    print(f"  Attempting to retrieve SMILES for CAS: '{cas_number}'")

    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/100.0.4896.88 Safari/537.36'
    }

    # Step 1: Get CID from CAS using xref/rn
    url_cas_to_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/xref/rn/{cas_number}/cids/JSON"
    try:
        print(f"    Trying CAS to CID lookup via xref/rn: {url_cas_to_cid}")
        response_cid = requests.get(url_cas_to_cid, timeout=10, headers=headers, verify=False)
        print(f"    CID lookup - Response status code: {response_cid.status_code}")
        print(f"    CID lookup - Response text (first 500 chars): {response_cid.text[:500]}")

        response_cid.raise_for_status() # Raise an exception for HTTP errors (4xx or 5xx)
        data_cid = response_cid.json()

        if 'IdentifierList' in data_cid and 'CID' in data_cid['IdentifierList'] and data_cid['IdentifierList']['CID']:
            cid = data_cid['IdentifierList']['CID'][0] # Take the first CID if multiple are returned
            print(f"    CID '{cid}' found for CAS '{cas_number}'.")

            # Step 2: Get SMILES from CID using TXT format for single property, using 'SMILES' property
            url_cid_to_smiles = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/SMILES/TXT" # Changed to 'SMILES'
            print(f"    Trying CID to SMILES lookup (TXT format, 'SMILES' property): {url_cid_to_smiles}")
            response_smiles = requests.get(url_cid_to_smiles, timeout=10, headers=headers, verify=False)
            print(f"    SMILES lookup - Response status code: {response_smiles.status_code}")
            print(f"    SMILES lookup - Response text (first 500 chars): {response_smiles.text[:500]}")

            response_smiles.raise_for_status()
            smiles = response_smiles.text.strip() # Get plain text and strip whitespace

            if smiles:
                print(f"    SMILES found for '{cas_number}'.")
                return smiles
            else:
                print(f"    SMILES property not found in TXT response for CID '{cid}' derived from CAS '{cas_number}'.")
        else:
            print(f"    CID not found in response for CAS '{cas_number}' via xref/rn.")

    except requests.exceptions.RequestException as e:
        print(f"    API request failed for CAS '{cas_number}': {e}")
    except (KeyError, IndexError):
        print(f"    CID property not found or unexpected JSON structure for CAS '{cas_number}'.")
    except Exception as e:
        print(f"    An unexpected error occurred for CAS '{cas_number}': {e}")

    print(f"  All attempts failed to retrieve SMILES for CAS '{cas_number}'.")
    return None

def get_molecules_from_file(filepath):
    """
    Reads molecule identifiers (e.g., CAS numbers) from a text file, one per line.
    """
    identifiers = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                identifier = line.strip()
                if identifier:
                    identifiers.append(identifier)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found. Please ensure it's in the correct directory.")
        return []
    except Exception as e:
        print(f"An error occurred while reading the file '{filepath}': {e}")
        return []
    return identifiers

def generate_molecule_png_white_transparent(identifier, output_folder="molecules_png_white_transparent"):
    """
    Generates a PNG image of a molecule with a white structure on a transparent background.
    """
    print(f"  [generate_molecule_png] Attempting to process identifier: {identifier}")

    try:
        os.makedirs(output_folder, exist_ok=True)
        print(f"  [generate_molecule_png] Ensured output folder exists: {output_folder}")

        smiles = get_smiles_from_cas(identifier)
        if not smiles:
            print(f"  [generate_molecule_png] No SMILES retrieved for '{identifier}', skipping image generation.")
            return

        print(f"  [generate_molecule_png] SMILES retrieved for '{identifier}': {smiles}")

        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                print(f"  [generate_molecule_png] Could not create RDKit molecule object from SMILES for CAS '{identifier}'.")
                return

            AllChem.Compute2DCoords(mol)
            print(f"  [generate_molecule_png] RDKit molecule created and 2D coords computed for '{identifier}'.")

            # --- Phase 1: SVG Drawing with RDKit (black structure on white background) ---
            drawer = rdMolDraw2D.MolDraw2DSVG(300, 300)
            drawer.drawOptions().bgColor = (1, 1, 1, 1)
            black = (0, 0, 0, 1)
            drawer.drawOptions().defaultColor = black
            drawer.drawOptions().atomPalette = {
                6: black, 7: black, 8: black, 9: black, 15: black,
                16: black, 17: black, 35: black, 53: black, 1: black,
            }
            drawer.drawOptions().bondPalette = {
                Chem.rdchem.BondType.SINGLE: black, Chem.rdchem.BondType.DOUBLE: black,
                Chem.rdchem.BondType.TRIPLE: black, Chem.rdchem.BondType.AROMATIC: black,
            }
            drawer.drawOptions().bondLineWidth = 2
            drawer.drawOptions().atomLabelFontSize = 18
            drawer.drawOptions().addStereoAnnotation = True
            drawer.drawOptions().useWedgeBond = True
            drawer.drawOptions().prepareMolsForDrawing = True
            drawer.drawOptions().addAtomIndices = False

            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            svg_data = drawer.GetDrawingText()
            print(f"  [generate_molecule_png] SVG data generated for '{identifier}'.")

            temp_output_filename = f"{identifier.replace(' ', '_')}_temp.png"
            temp_output_path_png = os.path.join(output_folder, temp_output_filename)
            print(f"  [generate_molecule_png] Temporary PNG path: {temp_output_path_png}")

            # Convert SVG to PNG with cairosvg
            cairosvg.svg2png(bytestring=svg_data.encode('utf-8'), write_to=temp_output_path_png,
                             output_width=300, output_height=300)
            print(f"  [generate_molecule_png] Temporary PNG created via cairosvg for '{identifier}'.")

            # --- Phase 2: Post-processing with Pillow for transparency and white colors ---
            final_output_filename = f"{identifier.replace(' ', '_')}.png"
            final_output_path_png = os.path.join(output_folder, final_output_filename)
            print(f"  [generate_molecule_png] Final PNG path: {final_output_path_png}")

            img = Image.open(temp_output_path_png).convert("RGBA")
            datas = img.getdata()

            new_data = []
            white_bg_threshold = 245
            black_mol_threshold = 10

            for item in datas:
                r, g, b, a = item
                if r > white_bg_threshold and g > white_bg_threshold and b > white_bg_threshold:
                    new_data.append((255, 255, 255, 0)) # Keep white but make transparent
                elif r < black_mol_threshold and g < black_mol_threshold and b < black_mol_threshold:
                    new_data.append((255, 255, 255, 255)) # Opaque white
                else:
                    new_data.append((255, 255, 255, a)) # Convert to white, keep original alpha
            
            img.putdata(new_data)
            img.save(final_output_path_png, "PNG")
            print(f"  [generate_molecule_png] Final PNG saved for '{identifier}'.")

            # Remove the temporary PNG file
            if os.path.exists(temp_output_path_png):
                os.remove(temp_output_path_png)
                print(f"  [generate_molecule_png] Removed temporary PNG for '{identifier}'.")
            else:
                print(f"  [generate_molecule_png] Warning: Temporary PNG for '{identifier}' not found for removal.")

            print(f"Structure for CAS '{identifier}' saved to '{final_output_path_png}' with white structure and transparent background.")
        except Exception as e:
            print(f"  [generate_molecule_png] Inner Error generating image for CAS '{identifier}': {e}")
    except Exception as e:
        print(f"  [generate_molecule_png] Outer Error processing CAS '{identifier}': {e}")

# --- Main execution block ---
molecules_filepath = "molecules.txt"

print(f"Attempting to read molecules from: {molecules_filepath}")
cas_numbers_to_generate = get_molecules_from_file(molecules_filepath)
print(f"Found CAS numbers: {cas_numbers_to_generate}")

if cas_numbers_to_generate:
    print(f"Starting image generation for {len(cas_numbers_to_generate)} molecules.")
    for cas_num in cas_numbers_to_generate:
        generate_molecule_png_white_transparent(cas_num)
    print("Finished processing all CAS numbers.")
else:
    print(f"No CAS numbers found in '{molecules_filepath}' or file not accessible. Please check the file and path.")

print("--- Script execution finished ---")