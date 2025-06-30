import requests
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import AllChem
import os
import cairosvg # Imports cairosvg for SVG to PNG conversion (still useful for the initial conversion, but not for final transparency)
from PIL import Image # Imports Pillow for image manipulation

def get_smiles_from_cas(cas_number):
    """
    Retrieves the canonical SMILES string of a molecule from PubChem by its CAS number.
    Includes multiple fallback mechanisms if direct lookup fails.
    """
    # Attempt 1: Direct SMILES lookup by CAS number
    url_direct_smiles = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cas/{cas_number}/property/CanonicalSMILES/JSON"
    try:
        response = requests.get(url_direct_smiles, timeout=10) # Added timeout
        response.raise_for_status() # Raise an exception for HTTP errors (e.g., 400, 404)
        data = response.json()
        
        # Check if PropertyTable and Properties exist and are not empty
        if 'PropertyTable' in data and 'Properties' in data['PropertyTable'] and data['PropertyTable']['Properties']:
            smiles = data['PropertyTable']['Properties'][0]['CanonicalSMILES']
            print(f"SMILES found via direct CAS lookup for '{cas_number}'.")
            return smiles
        else:
            print(f"Direct SMILES lookup for CAS '{cas_number}' returned no properties. Trying CID fallback...")
    except requests.exceptions.RequestException as e: # Catch all request exceptions (ConnectionError, HTTPError, Timeout, HTTPError)
        print(f"Direct SMILES lookup failed for CAS '{cas_number}': {e}. Trying CID fallback...")
    except (KeyError, IndexError):
        print(f"SMILES property not found or unexpected JSON structure in direct lookup for CAS '{cas_number}'. Trying CID fallback...")

    # Attempt 2: Fallback - Get CID by CAS number, then get SMILES by CID
    url_get_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cas/{cas_number}/cids/JSON"
    try:
        response_cid = requests.get(url_get_cid, timeout=10) # Added timeout
        response_cid.raise_for_status()
        data_cid = response_cid.json()
        
        # Check if IdentifierList and CID exist and are not empty
        if 'IdentifierList' in data_cid and 'CID' in data_cid['IdentifierList'] and data_cid['IdentifierList']['CID']:
            cid = data_cid['IdentifierList']['CID'][0] # Take the first CID if multiple are returned
            print(f"Found CID: {cid} for CAS '{cas_number}'. Now getting SMILES via CID...")

            # Now get SMILES using the retrieved CID
            url_smiles_from_cid = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON"
            response_smiles = requests.get(url_smiles_from_cid, timeout=10) # Added timeout
            response_smiles.raise_for_status()
            data_smiles = response_smiles.json()
            
            # Check if PropertyTable and Properties exist and are not empty
            if 'PropertyTable' in data_smiles and 'Properties' in data_smiles['PropertyTable'] and data_smiles['PropertyTable']['Properties']:
                smiles_from_cid = data_smiles['PropertyTable']['Properties'][0]['CanonicalSMILES']
                print(f"SMILES found via CID fallback for '{cas_number}'.")
                return smiles_from_cid
            else:
                print(f"SMILES property not found via CID '{cid}' for CAS '{cas_number}'. Trying name search fallback...")
        else:
            print(f"No CID found for CAS '{cas_number}'. Trying name search fallback...")
    except requests.exceptions.RequestException as e:
        print(f"Error in CID fallback for CAS '{cas_number}': {e}. Trying name search fallback...")
    except (KeyError, IndexError):
        print(f"CID or SMILES property structure unexpected in CID fallback for CAS '{cas_number}'. Trying name search fallback...")

    # Attempt 3: Fallback - Search by "name" (using CAS number as name)
    url_search_by_name = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{cas_number}/property/CanonicalSMILES/JSON"
    print(f"Attempting to find SMILES for CAS '{cas_number}' by treating it as a compound name...")
    try:
        response_name = requests.get(url_search_by_name, timeout=10)
        response_name.raise_for_status()
        data_name = response_name.json()

        if 'PropertyTable' in data_name and 'Properties' in data_name['PropertyTable'] and data_name['PropertyTable']['Properties']:
            smiles_from_name = data_name['PropertyTable']['Properties'][0]['CanonicalSMILES']
            print(f"SMILES found via name search fallback for '{cas_number}'.")
            return smiles_from_name
        else:
            print(f"SMILES property not found via name search for CAS '{cas_number}'.")
    except requests.exceptions.RequestException as e:
        print(f"Error in name search fallback for CAS '{cas_number}': {e}")
    except (KeyError, IndexError):
        print(f"SMILES property not found or unexpected JSON structure in name search fallback for CAS '{cas_number}'.")

    return None # If all attempts fail to retrieve SMILES


def get_molecules_from_file(filepath):
    """
    Reads molecule identifiers (e.g., CAS numbers) from a text file, one per line.

    Args:
        filepath (str): The path to the text file containing molecule identifiers.

    Returns:
        list: A list of molecule identifiers.
    """
    identifiers = []
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            for line in f:
                identifier = line.strip() # Remove leading/trailing whitespace
                if identifier: # Add non-empty lines
                    identifiers.append(identifier)
    except FileNotFoundError:
        print(f"Error: The file '{filepath}' was not found.")
        return []
    except Exception as e:
        print(f"An error occurred while reading the file '{filepath}': {e}")
        return []
    return identifiers

def generate_molecule_png_white_transparent(identifier, output_folder="molecules_png_white_transparent"):
    """
    Generates a PNG image of a molecule with a white structure on a transparent background.
    This method first generates a black structure on a white background, then post-processes it
    using Pillow to invert colors and add transparency.

    Args:
        identifier (str): The identifier of the molecule (e.g., CAS number).
        output_folder (str): The folder where the generated PNG images will be saved.
    """
    os.makedirs(output_folder, exist_ok=True)

    # Use the identifier (CAS number) to get the SMILES string
    smiles = get_smiles_from_cas(identifier)
    if smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                AllChem.Compute2DCoords(mol)

                # --- Phase 1: SVG Drawing with RDKit (black structure on white background) ---
                # The goal is to obtain a well-defined molecule in black on a white background
                # to facilitate subsequent post-processing with Pillow.
                drawer = rdMolDraw2D.MolDraw2DSVG(300, 300) # Set image dimensions

                # Set the RDKit drawing background to solid white for now.
                drawer.drawOptions().bgColor = (1, 1, 1, 1) # Opaque white

                # Set structure colors (atoms, bonds, text) to solid black.
                black = (0, 0, 0, 1) # Opaque black

                drawer.drawOptions().defaultColor = black
                drawer.drawOptions().atomPalette = {
                    6: black,  # Carbon (C)
                    7: black,  # Nitrogen (N)
                    8: black,  # Oxygen (O)
                    9: black,  # Fluorine (F)
                    15: black, # Phosphorus (P)
                    16: black, # Sulfur (S)
                    17: black, # Chlorine (Cl)
                    35: black, # Bromine (Br)
                    53: black, # Iodine (I)
                    1: black,  # Hydrogen (H)
                }
                drawer.drawOptions().bondPalette = {
                    Chem.rdchem.BondType.SINGLE: black,
                    Chem.rdchem.BondType.DOUBLE: black,
                    Chem.rdchem.BondType.TRIPLE: black,
                    Chem.rdchem.BondType.AROMATIC: black,
                }

                # Other drawing options (maintained for quality)
                drawer.drawOptions().bondLineWidth = 2
                drawer.drawOptions().atomLabelFontSize = 18
                drawer.drawOptions().addStereoAnnotation = True
                drawer.drawOptions().useWedgeBond = True
                drawer.drawOptions().prepareMolsForDrawing = True
                drawer.drawOptions().addAtomIndices = False

                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()

                svg_data = drawer.GetDrawingText()

                # Path for the temporary image (black on white)
                # Use the identifier (CAS number) for the filename
                temp_output_path_png = os.path.join(output_folder, f"{identifier.replace(' ', '_')}_temp.png")

                # Convert SVG to PNG with cairosvg (no transparency here, handled by Pillow later)
                cairosvg.svg2png(bytestring=svg_data.encode('utf-8'), write_to=temp_output_path_png,
                                 output_width=300, output_height=300) # No background_color for now

                # --- Phase 2: Post-processing with Pillow for transparency and white colors ---
                # Use the identifier (CAS number) for the final filename
                final_output_path_png = os.path.join(output_folder, f"{identifier.replace(' ', '_')}.png")

                img = Image.open(temp_output_path_png).convert("RGBA")
                datas = img.getdata()

                new_data = []
                # Define a tolerance for colors due to anti-aliasing
                white_bg_threshold = 245 # Almost white
                black_mol_threshold = 10 # Almost black

                for item in datas:
                    r, g, b, a = item
                    # If the pixel is close to white (background), make it transparent
                    if r > white_bg_threshold and g > white_bg_threshold and b > white_bg_threshold:
                        new_data.append((255, 255, 255, 0)) # Keep white but make transparent
                    # If the pixel is close to black (structure), make it opaque white
                    elif r < black_mol_threshold and g < black_mol_threshold and b < black_mol_threshold:
                        new_data.append((255, 255, 255, 255)) # Opaque white
                    else:
                        # For intermediate colors (anti-aliasing), convert to white while preserving alpha
                        # You can refine this or convert to pure opaque white if you prefer a sharper look
                        new_data.append((255, 255, 255, a)) # Convert to white, keep original alpha
                
                img.putdata(new_data)
                img.save(final_output_path_png, "PNG")

                # Remove the temporary PNG file
                os.remove(temp_output_path_png)

                print(f"Structure for CAS '{identifier}' saved to '{final_output_path_png}' with white structure and transparent background.")
            else:
                print(f"Could not create molecule object from SMILES for CAS '{identifier}'.")
        except Exception as e:
            print(f"Error generating image for CAS '{identifier}': {e}")
    else:
        print(f"Could not generate image for CAS '{identifier}' without a valid SMILES string.")

# --- Usage Example: Read molecules from a text file ---
# Specify the path to your text file here
# This file should contain CAS numbers, one per line.
molecules_filepath = "molecules.txt" # Change this to your actual file path

# Get the list of CAS numbers from the file
cas_numbers_to_generate = get_molecules_from_file(molecules_filepath)

if cas_numbers_to_generate:
    for cas_num in cas_numbers_to_generate:
        generate_molecule_png_white_transparent(cas_num)
else:
    print(f"No CAS numbers found in '{molecules_filepath}' or file not accessible. Please check the file and path.")
