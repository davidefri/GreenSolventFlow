// JS/table_interactivity.js

// Le variabili ALL_DATA, ALL_COLUMNS, DEFAULT_VISIBLE_COLUMNS, H_PHRASES_DEFINITIONS, STRUCTURE_IMAGES_FOLDER
// sono iniettate direttamente dall'HTML generato dal Python.

const searchInput = document.getElementById('searchInput');
const tableBody = document.getElementById('tableBody');
const tableHeaderRow = document.getElementById('tableHeaderRow');
const columnSelectionPanel = document.getElementById('column-selection-panel');
const columnCheckboxesDiv = document.getElementById('columnCheckboxes');
const applyColumnsBtn = document.getElementById('applyColumnsBtn');
const toggleColumnsBtn = document.getElementById('toggleColumnsBtn');
const mainTable = document.getElementById('mainTable');

let currentData = [...ALL_DATA]; // Copia i dati originali per permettere modifiche (es. ordinamento, filtro)
let visibleColumns = DEFAULT_VISIBLE_COLUMNS || [];
let sortColumn = 'Name'; // Colonna di default per l'ordinamento
let sortDirection = 'asc'; // Direzione di default

/**
 * Applica i tooltip alle frasi H/EUH nel testo.
 * @param {string} text Il testo contenente le frasi H/EUH.
 * @returns {string} Il testo con i tooltip applicati.
 */
function applyHPhraseTooltipToText(text) {
    if (typeof text !== 'string') {
        return text; // Ritorna il valore così com'è se non è una stringa
    }
    let modifiedText = text;
    // Regex per trovare codici H (es. H300, H350i) o EUH (es. EUH202, EUH208)
    // Ordina i match per lunghezza decrescente per evitare sostituzioni parziali (es. H360FD prima di H360F)
    const matches = Array.from(new Set(modifiedText.match(/\b(EUH\d{3}[A-Za-z0-9]*|H\d{3}[A-Za-z0-9]*)\b/g) || []))
                             .sort((a, b) => b.length - a.length);

    for (const hCode of matches) {
        const definition = H_PHRASES_DEFINITIONS[hCode];
        if (definition) {
            // Usa una funzione di sostituzione per gestire più occorrenze e prevenire problemi con caratteri speciali
            modifiedText = modifiedText.replace(new RegExp('\\b' + hCode.replace(/[.*+?^${}()|[\]\\]/g, '\\$&') + '\\b', 'g'), 
                                                `<span class="h-phrase-tooltip" title="${definition}">${hCode}</span>`);
        }
    }
    return modifiedText;
}


/**
 * Popola la tabella HTML con i dati filtrati e ordinati.
 */
function populateTable() {
    tableBody.innerHTML = ''; // Pulisci il corpo della tabella
    tableHeaderRow.innerHTML = ''; // Pulisci l'intestazione della tabella

    // Crea le intestazioni della tabella
    visibleColumns.forEach(column => {
        const th = document.createElement('th');
        th.textContent = column;
        th.setAttribute('data-column', column); // Aggiungi un attributo per identificare la colonna
        th.addEventListener('click', () => sortTable(column));
        
        // Aggiungi la classe per l'indicatore di ordinamento
        if (column === sortColumn) {
            th.classList.add(sortDirection === 'asc' ? 'sorted-asc' : 'sorted-desc');
        }
        tableHeaderRow.appendChild(th);
    });

    currentData.forEach(rowData => {
        const row = document.createElement('tr');
        visibleColumns.forEach(column => {
            const cell = document.createElement('td');
            let cellContent = rowData[column];

            if (column === 'Name') {
                const link = document.createElement('a');
                link.href = rowData.detail_url;
                link.textContent = cellContent;
                cell.appendChild(link);
            } else if (column === 'Structure' && rowData.structure_img_src) {
                const img = document.createElement('img');
                img.src = rowData.structure_img_src;
                img.alt = `Structure of ${rowData.Name}`;
                img.classList.add('index-image'); // Applica la classe per lo stile
                cell.appendChild(img);
            } else {
                // Applica tooltip alle frasi H/EUH in tutte le altre colonne testuali
                cell.innerHTML = applyHPhraseTooltipToText(cellContent);
            }
            row.appendChild(cell);
        });
        tableBody.appendChild(row);
    });
}

/**
 * Filtra i dati della tabella in base all'input di ricerca.
 */
function filterTable() {
    const searchTerm = searchInput.value.toLowerCase();
    currentData = ALL_DATA.filter(row => {
        // Cerca in tutte le colonne visibili per ora
        // Potresti ottimizzare per cercare solo in colonne specifiche se necessario
        return visibleColumns.some(column => {
            const value = String(row[column] || '').toLowerCase();
            return value.includes(searchTerm);
        });
    });

    // Ordina i dati filtrati e poi ripopola la tabella
    sortTable(sortColumn, sortDirection, false); // Ordina i dati filtrati
    populateTable(); // QUESTA È LA RIGA AGGIUNTA/VERIFICATA!
}

/**
 * Ordina la tabella in base alla colonna cliccata.
 * @param {string} column La colonna su cui ordinare.
 * @param {string} [direction] Direzione di ordinamento ('asc' o 'desc'). Se non specificato, inverte la direzione corrente.
 * @param {boolean} [repopulate=true] Se ripopolare la tabella dopo l'ordinamento.
 */
function sortTable(column, direction = null, repopulate = true) {
    if (column === sortColumn) {
        sortDirection = direction || (sortDirection === 'asc' ? 'desc' : 'asc');
    } else {
        sortColumn = column;
        sortDirection = direction || 'asc';
    }

    currentData.sort((a, b) => {
        let valA = a[sortColumn];
        let valB = b[sortColumn];

        // Gestione dei valori numerici e mancanti
        if (typeof valA === 'string' && !isNaN(parseFloat(valA)) && !isNaN(parseFloat(valB))) {
            valA = parseFloat(valA);
            valB = parseFloat(valB);
        } else {
            valA = String(valA || '').toLowerCase();
            valB = String(valB || '').toLowerCase();
        }

        if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
        if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
        return 0;
    });

    if (repopulate) {
        populateTable();
    }
}


/**
 * Inizializza il pannello di selezione delle colonne con le checkbox.
 */
function initializeColumnSelection() {
    columnCheckboxesDiv.innerHTML = ''; // Pulisci
    ALL_COLUMNS.forEach(column => {
        const div = document.createElement('div');
        div.classList.add('column-checkbox-item');

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.id = `col-${column.replace(/\s+/g, '-')}`; // ID univoco per il checkbox
        checkbox.value = column;
        if (visibleColumns.has(column)) {
            checkbox.checked = true;
        }

        const label = document.createElement('label');
        label.htmlFor = checkbox.id;
        label.textContent = column;

        div.appendChild(checkbox);
        div.appendChild(label);
        columnCheckboxesDiv.appendChild(div);
    });
}

/**
 * Gestisce l'applicazione delle colonne selezionate.
 */
function applySelectedColumns() {
    const newVisibleColumns = new Set();
    columnCheckboxesDiv.querySelectorAll('input[type="checkbox"]:checked').forEach(checkbox => {
        newVisibleColumns.add(checkbox.value);
    });
    visibleColumns = newVisibleColumns;
    populateTable(); // Ripopola la tabella con le nuove colonne visibili
    toggleColumnSelectionPanel(); // Chiudi il pannello
}

/**
 * Toggle la visibilità del pannello di selezione colonne e dell'overlay.
 */
function toggleColumnSelectionPanel() {
    const overlay = document.getElementById('overlay');
    if (columnSelectionPanel.style.display === 'block') {
        columnSelectionPanel.style.display = 'none';
        overlay.style.display = 'none';
    } else {
        initializeColumnSelection(); // Re-inizializza per riflettere lo stato corrente
        columnSelectionPanel.style.display = 'block';
        overlay.style.display = 'block';
    }
}

// --- Event Listeners ---
searchInput.addEventListener('input', filterTable);
applyColumnsBtn.addEventListener('click', applySelectedColumns);
toggleColumnsBtn.addEventListener('click', (e) => {
    e.preventDefault(); // Previene il comportamento di default del link
    toggleColumnSelectionPanel();
});

// Aggiungi un overlay per chiudere il pannello cliccando fuori
const overlay = document.createElement('div');
overlay.id = 'overlay';
document.body.appendChild(overlay);
overlay.addEventListener('click', toggleColumnSelectionPanel);


// Inizializzazione della pagina
document.addEventListener('DOMContentLoaded', () => {
    populateTable(); // Popola la tabella iniziale
    // La selezione delle colonne viene inizializzata solo quando il pannello viene aperto
});

// Imposta lo stato iniziale dell'ordinamento per la colonna di default
sortTable(sortColumn, sortDirection, false); // Ordina i dati senza ripopolare subito