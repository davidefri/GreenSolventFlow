import pandas as pd
import os
import re
import json

# Dizionario delle definizioni delle frasi H e EUH in INGLESE
h_phrases_definitions = {
    "H200": "Unstable explosive", "H201": "Explosive; mass explosion hazard",
    "H202": "Explosive; severe projection hazard", "H203": "Explosive; fire, blast or projection hazard",
    "H204": "Fire or projection hazard", "H205": "May mass explode in fire",
    "H220": "Extremely flammable gas", "H221": "Flammable gas",
    "H222": "Extremely flammable aerosol", "H223": "Flammable aerosol",
    "H224": "Extremely flammable liquid and vapour", "H225": "Highly flammable liquid and vapour",
    "H226": "Flammable liquid and vapour", "H228": "Flammable solid",
    "H229": "Pressurised container: May burst if heated", "H230": "May react explosively even in absence of air",
    "H231": "May react explosively even in absence of air at elevated pressure and/or temperature",
    "H232": "Catches fire spontaneously if exposed to air", "H240": "Heating may cause an explosion",
    "H241": "Heating may cause a fire or explosion", "H242": "Heating may cause a fire",
    "H250": "Catches fire spontaneously if exposed to air", "H251": "Self-heating: may catch fire",
    "H252": "Self-heating in large quantities; may catch fire",
    "H260": "In contact with water releases flammable gases which may ignite spontaneously",
    "H261": "In contact with water releases flammable gases", "H270": "May cause or intensify fire; oxidiser",
    "H271": "May cause fire or explosion; strong oxidiser", "H272": "May intensify fire; oxidiser",
    "H280": "Contains gas under pressure; may explode if heated",
    "H281": "Contains refrigerated gas; may cause cryogenic burns or injury",
    "H290": "May be corrosive to metals", "H300": "Fatal if swallowed",
    "H301": "Toxic if swallowed", "H302": "Harmful if swallowed",
    "H304": "May be fatal if swallowed and enters airways", "H310": "Fatal in contact with skin",
    "H311": "Toxic in contact with skin", "H312": "Harmful in contact with skin",
    "H314": "Causes severe skin burns and eye damage", "H315": "Causes skin irritation",
    "H317": "May cause an allergic skin reaction", "H318": "Causes serious eye damage",
    "H319": "Causes serious eye irritation", "H330": "Fatal if inhaled",
    "H331": "Toxic if inhaled", "H332": "Harmful if inhaled",
    "H334": "May cause allergy or asthma symptoms or breathing difficulties if inhaled",
    "H335": "May cause respiratory irritation", "H336": "May cause drowsiness or dizziness",
    "H340": "May cause genetic defects", "H341": "Suspected of causing genetic defects",
    "H350": "May cause cancer", "H350i": "May cause cancer by inhalation",
    "H351": "Suspected of causing cancer", "H360": "May damage fertility or the unborn child",
    "H360F": "May damage fertility", "H360D": "May damage the unborn child",
    "H360FD": "May damage fertility. May damage the unborn child",
    "H360Fd": "May damage fertility. Suspected of damaging the unborn child",
    "H360Df": "May damage the unborn child. Suspected of damaging fertility",
    "H361": "Suspected of damaging fertility or the unborn child", "H361f": "Suspected of damaging fertility",
    "H361d": "Suspected of damaging the unborn child", "H361fd": "Suspected of damaging fertility and the unborn child",
    "H362": "May cause harm to breast-fed children", "H370": "Causes damage to organs",
    "H371": "May cause damage to organs",
    "H372": "Causes damage to organs through prolonged or repeated exposure",
    "H373": "May cause damage to organs through prolonged or repeated exposure",
    "H400": "Very toxic to aquatic life", "H410": "Very toxic to aquatic life with long lasting effects",
    "H411": "Toxic to aquatic life with long lasting effects", "H412": "Harmful to aquatic life with long lasting effects",
    "H413": "May cause long lasting harmful effects to aquatic life",
    "H420": "Harms public health and the environment by destroying ozone in the upper atmosphere",
    "EUH001": "Explosive when dry", "EUH006": "Explosive with or without contact with air",
    "EUH014": "Reacts violently with water", "EUH018": "In use, may form flammable/explosive vapour-air mixture",
    "EUH019": "May form explosive peroxides", "EUH044": "Risk of explosion if heated under confinement",
    "EUH029": "Contact with water liberates toxic gas", "EUH031": "Contact with acids liberates toxic gas",
    "EUH032": "Contact with acids liberates very toxic gas",
    "EUH066": "Repeated exposure may cause skin dryness or cracking", "EUH070": "Toxic by eye contact",
    "EUH071": "Corrosive to the respiratory tract",
    "EUH201": "Contains lead. Should not be used on surfaces liable to be chewed or sucked by children.",
    "EUH201A": "Contains lead.",
    "EUH202": "Cyanoacrylate. Danger. Bonds skin and eyes in seconds. Keep out of the reach of children.",
    "EUH203": "Contains chromium(VI). May produce an allergic reaction.",
    "EUH204": "Contains isocyanates. May produce an allergic reaction.",
    "EUH205": "Contains epoxy constituents. May produce an allergic reaction.",
    "EUH206": "Warning! Do not use together with other products. May release dangerous gases (chlorine).",
    "EUH207": "Warning! Contains cadmium. Dangerous fumes are formed during use. See information supplied by the manufacturer. Comply with the safety instructions.",
    "EUH208": "Contains (name of sensitising substance). May produce an allergic reaction.",
    "EUH209": "Can become highly flammable in use.", "EUH209A": "Can become flammable in use.",
    "EUH210": "Safety data sheet available on request.",
    "EUH211": "Warning! Hazardous respirable droplets may be formed when sprayed. Do not breathe spray or mist.",
    "EUH212": "Warning! Hazardous respirable dust may be formed when used. Do not breathe dust.",
    "EUH401": "To avoid risks to human health and the environment, comply with the instructions for use."
}

def apply_h_phrase_tooltip(text_content):
    """
    Cerca frasi H/EUH nel testo e le racchiude in uno <span> con l'attributo title e la classe 'h-phrase-tooltip'.
    Questa funzione è usata per le pagine di dettaglio generate da Python.
    La gestione dei tooltip per la tabella principale (index.html) è affidata a JavaScript.
    """
    if not isinstance(text_content, str):
        return str(text_content) if pd.notna(text_content) else ''
    
    pattern = r'\b(EUH\d{3}[A-Za-z0-9]*|H\d{3}[A-Za-z0-9]*)\b'
    
    matches = sorted(list(set(re.findall(pattern, text_content))), key=len, reverse=True)
    
    modified_text = text_content
    for h_code in matches:
        definition = h_phrases_definitions.get(h_code, None)
        if definition:
            modified_text = re.sub(r'\b' + re.escape(h_code) + r'\b', f'<span class="h-phrase-tooltip" title="{definition}">{h_code}</span>', modified_text)
    
    return modified_text


def generate_html_pages(csv_file='dati_solvente.csv', output_folder='output_pages', 
                        logo_header_path='PNG/IconGreenSolventFlow3.png', 
                        favicon_path='PNG/flaviconGreenSolventFlow2.png', 
                        css_path='Style.css', structure_images_folder='PNG/structures',
                        js_folder='JS'):
    """
    Genera una pagina HTML principale e pagine di dettaglio per ogni riga di un file CSV.
    Include la logica per la tabella principale dinamica (ricerca, ordinamento, visibilità colonne).

    Args:
        csv_file (str): Il percorso del file CSV contenente i dati.
        output_folder (str): La cartella in cui salvare le pagine HTML di dettaglio.
        logo_header_path (str): Il percorso dell'immagine del logo nell'header (relativo alla root del sito).
        favicon_path (str): Il percorso del file favicon (relativo alla root del sito).
        css_path (str): Il percorso del file CSS (relativo alla root del sito).
        structure_images_folder (str): La cartella delle immagini di struttura (relativa alla root del sito).
        js_folder (str): La cartella dei file JavaScript (relativo alla root del sito).
    """
    df = pd.read_csv(csv_file, sep=';', decimal=',', na_values=[''], encoding='latin-1', low_memory=False)

    df.columns = df.columns.str.strip()

    if 'Unnamed: 20' in df.columns and df['Unnamed: 20'].isnull().all():
        df = df.drop(columns=['Unnamed: 20'])
        print("[INFO] Rimossa colonna 'Unnamed: 20' in quanto vuota.")
    
    initial_rows = len(df)
    df['Name'] = df['Name'].replace('', pd.NA)
    df.dropna(subset=['Name'], inplace=True)
    if len(df) < initial_rows:
        print(f"[ATTENZIONE] Rimosse {initial_rows - len(df)} righe con valore 'Name' mancante o vuoto.")
    else:
        print(f"[INFO] Trovate {len(df)} righe valide nel CSV.")

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"[INFO] Creata cartella di output: '{output_folder}'")
    else:
        print(f"[INFO] Cartella di output '{output_folder}' già esistente.")

    data_for_js = df.fillna('').to_dict(orient='records')

    for item in data_for_js:
        unique_id_raw = str(item['Name']).strip()
        unique_id_filename = re.sub(r'[\\/:*?"<>|\s]', '_', unique_id_raw)
        unique_id_filename = re.sub(r'__+', '_', unique_id_filename)
        unique_id_filename = unique_id_filename.strip('_')
        
        if not unique_id_filename:
            unique_id_filename = f"unnamed_entry_{data_for_js.index(item)}"
            print(f"[AVVISO] Il nome '{unique_id_raw}' ha generato un nome file vuoto. Usato fallback: '{unique_id_filename}'")
        
        item['detail_url'] = f"{output_folder}/{unique_id_filename}.html"
        
        if 'Structure' in item and item['Structure']:
            item['structure_img_src'] = f"{structure_images_folder}/{item['Structure']}"
        else:
            item['structure_img_src'] = ''

    excluded_columns_from_selection = ['detail_url', 'structure_img_src']
    all_potential_index_columns = [col for col in df.columns.tolist() if col not in excluded_columns_from_selection]
    
    default_visible_index_columns = ['Name', 'CAS', 'Structure']

    json_data_for_js = json.dumps(data_for_js, ensure_ascii=False)
    json_all_columns = json.dumps(all_potential_index_columns, ensure_ascii=False)
    json_default_columns = json.dumps(default_visible_index_columns, ensure_ascii=False)
    json_h_phrases_definitions = json.dumps(h_phrases_definitions, ensure_ascii=False)

    print("\n[INFO] Generazione delle pagine di dettaglio...")
    for index, row in df.iterrows():
        unique_id_raw = str(row['Name']).strip()
        unique_id_filename = re.sub(r'[\\/:*?"<>|\s]', '_', unique_id_raw)
        unique_id_filename = re.sub(r'__+', '_', unique_id_filename)
        unique_id_filename = unique_id_filename.strip('_')
        
        if not unique_id_filename:
            unique_id_filename = f"unnamed_entry_{index}"
            print(f"[AVVISO] Il nome '{unique_id_raw}' ha generato un nome file vuoto. Usato fallback: '{unique_id_filename}'")

        detail_filename = os.path.join(output_folder, f"{unique_id_filename}.html")

        detail_html_header = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SolventFlow - {unique_id_raw}</title>
    <link rel="shortcut icon" href="../{favicon_path}" type="image/x-icon">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Agdasima:wght@400;700&family=Megrim&display=swap" rel="stylesheet">
    <link rel="stylesheet" href="../{css_path}">
</head>
<body>
    <header>
        <a href="../index.html"> <img class="LogoHeader" src="../{logo_header_path}" alt="icon">
        </a>
    </header>
    <div class="detail-header-container">
        <a href="../index.html" class="back-to-list-button">
            <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" class="arrow-icon-left">
                <path fill-rule="evenodd" d="M12.53 16.28a.75.75 0 0 1-1.06 0l-7.5-7.5a.75.75 0 0 1 1.06-1.06L12 14.69l6.97-6.97a.75.75 0 1 1 1.06 1.06l-7.5 7.5Z" clip-rule="evenodd" />
            </svg>
            <span class="button-text">Back to List</span>
        </a>
        <h1>{unique_id_raw}</h1>
    </div>
    <table class="detailTable"> <tbody>
'''
        detail_table_rows = ""
        for col_name, value in row.items():
            if col_name in ['detail_url', 'structure_img_src']:
                continue

            display_value = ""
            if col_name == 'Structure':
                if pd.notna(value) and str(value).strip() != '':
                    detail_image_src = f"../{structure_images_folder}/{value}"
                    display_value = f'<img src="{detail_image_src}" alt="Structure of {unique_id_raw}" class="detail-image">'
                else:
                    display_value = "No structure available"
            else:
                display_value = "" if pd.isna(value) else str(value)
                display_value = apply_h_phrase_tooltip(display_value)
            
            detail_table_rows += f'''<tr>
                <th>{col_name}</th>
                <td>{display_value}</td>
            </tr>
'''
        
        detail_html_footer = f'''
        </tbody>
    </table>
    <button id="scrollToTopBtn" title="Back to Top">
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" class="arrow-icon">
        <path fill-rule="evenodd" d="M11.47 2.47a.75.75 0 0 1 1.06 0l7.5 7.5a.75.75 0 1 1-1.06 1.06L12 4.81V18.75a.75.75 0 0 1-1.5 0V4.81L4.03 11.03a.75.75 0 0 1-1.06-1.06l7.5-7.5Z" clip-rule="evenodd" />
    </svg>
</button>
    <script src="../{js_folder}/scroll.js"></script> </body>
</html>
'''
        try:
            with open(detail_filename, 'w', encoding='utf-8') as f:
                f.write(detail_html_header + detail_table_rows + detail_html_footer)
        except OSError as e:
            print(f"[ERRORE] Impossibile scrivere il file di dettaglio '{detail_filename}': {e}")
            continue
    print("[INFO] Pagine di dettaglio generate.")

    print("[INFO] Generazione del file 'index.html'...")

    full_main_html_content = f'''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>SolventFlow</title>
    <link rel="shortcut icon" href="{favicon_path}" type="image/x-icon">
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Agdasima:wght@400;700&family=Megrim&display=swap" rel="stylesheet">
    <link rel="stylesheet" href="{css_path}">
</head>
<body>
    <header>
        <a href="index.html"> <img class="LogoHeader" src="{logo_header_path}" alt="icon">
        </a>
        <li><a href="mailto:davide.frigatti@gmail.com" target="_blank">Send us an Email</a></li>
    </header>
    <div class="search-container">
        <input type="text" id="searchInput" placeholder="Search by name, CAS, etc.">
        <button id="clearSearchBtn" class="clear-search-btn" style="display: none;">&times;</button>
    </div>
        </div>

    <button id="floatingButton" class="floating-button">
        Filters
    </button>

    <div id="floatingPanel" class="floating-panel">
        <button class="close-button" id="closePanel">&times;</button>
        <h2>Filter your results</h2>
        <button class="panel-button" id="toggleColumnsBtn">Select Columns</button>
        <button class="panel-button">Opzione 2</button>
        <button class="panel-button">Opzione 3</button>
        <button class="panel-button">Contattaci</button>
    </div>

    <script src="JS/panel.js"></script>
        

    <div id="column-selection-panel" class="column-selection-panel">
        <h2>Select Columns to Display</h2>
        <div id="columnCheckboxes">
            {{}} </div>
        <button id="applyColumnsBtn">Apply</button>
    </div>

    <h1>List of Solvents</h1>
    <table class="table" id="mainTable">
        <thead>
            <tr id="tableHeaderRow">
                </tr>
        </thead>
        <tbody id="tableBody">
            </tbody>
    </table>
    <button id="scrollToTopBtn" title="Back to Top">
    <svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="currentColor" class="arrow-icon">
        <path fill-rule="evenodd" d="M11.47 2.47a.75.75 0 0 1 1.06 0l7.5 7.5a.75.75 0 1 1-1.06 1.06L12 4.81V18.75a.75.75 0 0 1-1.5 0V4.81L4.03 11.03a.75.75 0 0 1-1.06-1.06l7.5-7.5Z" clip-rule="evenodd" />
    </svg>
</button>
    <script>
        // Dati iniettati direttamente dal Python per il JavaScript
        const ALL_DATA = {json_data_for_js};
        const ALL_COLUMNS = {json_all_columns};
        const DEFAULT_VISIBLE_COLUMNS = {json_default_columns};
        const H_PHRASES_DEFINITIONS = {json_h_phrases_definitions};
        const STRUCTURE_IMAGES_FOLDER = '{structure_images_folder}';

        // Aggiungi questi console.log per la diagnostica
        console.log("DEFAULT_VISIBLE_COLUMNS (da Python nel HTML):", DEFAULT_VISIBLE_COLUMNS);
        console.log("Contenuto di H_PHRASES_DEFINITIONS:", H_PHRASES_DEFINITIONS);
        console.log("Script block custom diagnostics loaded.");
    </script>
    <script src="{js_folder}/table_interactivity.js"></script>
    <script src="{js_folder}/scroll.js"></script>
</body>
</html>
'''
    try:
        with open('index.html', 'w', encoding='utf-8') as f:
            f.write(full_main_html_content)
        print("[INFO] File 'index.html' generated successfully.")
    except OSError as e:
        print(f"[ERRORE] Impossibile scrivere il file 'index.html': {e}")

    print("[*] Website generated successfully!")
    
# Esegui la funzione per generare il sito
generate_html_pages()