// JS/table_interactivity.js

// Le variabili ALL_DATA, ALL_COLUMNS, DEFAULT_VISIBLE_COLUMNS, H_PHRASES_DEFINITIONS, STRUCTURE_IMAGES_FOLDER
// sono iniettate direttamente nell'HTML generato dal Python.
// Non devono essere dichiarate di nuovo qui (es. con 'var', 'let', 'const'),
// altrimenti si verificherebbe un errore "has already been declared".
// Saranno disponibili globalmente grazie alla dichiarazione nel tag <script> inline dell'HTML.

// Recupera i riferimenti agli elementi DOM una volta sola all'inizio
const searchInput = document.getElementById('searchInput');
const tableBody = document.getElementById('tableBody');
const tableHeaderRow = document.getElementById('tableHeaderRow');
const columnSelectionPanel = document.getElementById('column-selection-panel');
const columnCheckboxesDiv = document.getElementById('columnCheckboxes');
const applyColumnsBtn = document.getElementById('applyColumnsBtn');
const toggleColumnsBtn = document.getElementById('toggleColumnsBtn');
const mainTable = document.getElementById('mainTable');

// --- INIZIALIZZAZIONE DELLE VARIABILI DI STATO ---
// Queste variabili gestiscono lo stato della tabella e devono essere dichiarate in questo file.
// Assicurati che DEFAULT_VISIBLE_COLUMNS e ALL_DATA siano stati correttamente iniettati come array/oggetti JS validi nell'HTML.
let visibleColumns = DEFAULT_VISIBLE_COLUMNS || []; // Deve essere un array per usare .some() e .includes()
let currentData = [...ALL_DATA]; // Copia i dati originali per permettere modifiche (es. ordinamento, filtro)
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
    // Utilizza un Set temporaneo per gestire i match unici, poi converti in array per ordinare
    const matches = Array.from(new Set(modifiedText.match(/\b(EUH\d{3}[A-Za-z0-9]*|H\d{3}[A-Za-z0-9]*)\b/g) || []))
                                 .sort((a, b) => b.length - a.length);

    for (const hCode of matches) {
        const definition = H_PHRASES_DEFINITIONS[hCode];
        if (definition) {
            // Utilizza una RegExp globale per sostituire tutte le occorrenze
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
        th.setAttribute('data-column', column);
        th.addEventListener('click', () => sortTable(column));

        if (column === sortColumn) {
            th.classList.add(sortDirection === 'asc' ? 'sorted-asc' : 'sorted-desc');
        }
        tableHeaderRow.appendChild(th);
    });

    // Popola le righe della tabella
    currentData.forEach(rowData => {
        const row = document.createElement('tr');
        visibleColumns.forEach(column => {
            const cell = document.createElement('td');
            let cellContent = rowData[column];

            if (column === 'Name') {
                const link = document.createElement('a');
                link.href = rowData.detail_url; // Assumi che detail_url esista in rowData
                link.textContent = cellContent;
                cell.appendChild(link);
            } else if (column === 'Structure' && rowData.structure_img_src) {
                const img = document.createElement('img');
                img.src = rowData.structure_img_src; // Assumi che structure_img_src esista in rowData
                img.alt = `Structure of ${rowData.Name}`;
                img.classList.add('index-image'); // Classe per lo stile CSS delle immagini
                cell.appendChild(img);
            } else {
                // Applica i tooltip solo se il contenuto è una stringa
                cell.innerHTML = applyHPhraseTooltipToText(cellContent);
            }
            row.appendChild(cell);
        });
        tableBody.appendChild(row);
    });
}

/**
 * Filtra i dati della tabella in base all'input di ricerca corrente.
 * Questa funzione è fondamentale per mantenere la coerenza dello stato.
 */
function filterTable() {
    const searchTerm = searchInput.value.toLowerCase();
    // Filtra sempre partendo dai dati originali ALL_DATA
    currentData = ALL_DATA.filter(row => {
        // visibleColumns è un array, quindi il metodo .some() funziona correttamente
        return visibleColumns.some(column => {
            const value = String(row[column] || '').toLowerCase(); // Assicurati che il valore sia una stringa
            return value.includes(searchTerm);
        });
    });

    // Dopo il filtro, ordina i dati e poi ripopola la tabella
    sortTable(sortColumn, sortDirection, false); // Ordina i dati filtrati (senza ripopolare subito)
    populateTable(); // Popola la tabella con i dati filtrati e ordinati
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

        // Gestione dell'ordinamento numerico per stringhe che possono essere numeri
        if (typeof valA === 'string' && typeof valB === 'string' &&
            !isNaN(parseFloat(valA)) && !isNaN(parseFloat(valB)) &&
            valA.trim() !== '' && valB.trim() !== '') { // Controlla anche che non siano stringhe vuote
            valA = parseFloat(valA);
            valB = parseFloat(valB);
        } else {
            // Conversione a stringa e minuscolo per ordinamento alfabetico standard
            valA = String(valA || '').toLowerCase();
            valB = String(valB || '').toLowerCase();
        }

        // Logica di confronto
        if (valA < valB) return sortDirection === 'asc' ? -1 : 1;
        if (valA > valB) return sortDirection === 'asc' ? 1 : -1;
        return 0; // Se i valori sono uguali
    });

    if (repopulate) {
        populateTable();
    }
}


/**
 * Inizializza il pannello di selezione delle colonne con le checkbox.
 */
function initializeColumnSelection() {
    columnCheckboxesDiv.innerHTML = ''; // Pulisci i checkbox esistenti

    // Genera una checkbox per ogni colonna disponibile in ALL_COLUMNS
    ALL_COLUMNS.forEach(column => {
        const div = document.createElement('div');
        div.classList.add('column-checkbox-item');

        const checkbox = document.createElement('input');
        checkbox.type = 'checkbox';
        checkbox.id = `col-${column.replace(/\s+/g, '-')}`; // ID unico per il checkbox
        checkbox.value = column;
        // Se la colonna è attualmente visibile, spunta la checkbox
        if (visibleColumns.includes(column)) {
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
 * Gestisce l'applicazione delle colonne selezionate dall'utente.
 */
function applySelectedColumns() {
    const newVisibleColumns = new Set(); // Usa un Set temporaneo per raccogliere i valori unici dai checkbox
    columnCheckboxesDiv.querySelectorAll('input[type="checkbox"]:checked').forEach(checkbox => {
        newVisibleColumns.add(checkbox.value);
    });
    // Converti il Set in un Array e assegnarlo a visibleColumns per coerenza
    visibleColumns = Array.from(newVisibleColumns);

    // --- PUNTO CHIAVE PER LA COERENZA DELLO STATO ---
    // Dopo aver modificato le colonne visibili, è necessario ricalcolare l'intera tabella.
    // Chiamando filterTable(), si riapplica il filtro di ricerca corrente,
    // l'ordinamento, e infine si ridisegna la tabella con populateTable().
    filterTable();

    toggleColumnSelectionPanel(); // Chiudi il pannello dopo l'applicazione
}

/**
 * Toggle la visibilità del pannello di selezione colonne e dell'overlay.
 */
function toggleColumnSelectionPanel() {
    const overlay = document.getElementById('overlay');
    if (columnSelectionPanel.style.display === 'block') {
        // Se il pannello è visibile, nascondilo
        columnSelectionPanel.style.display = 'none';
        overlay.style.display = 'none';
    } else {
        // Se il pannello è nascosto, mostralo e inizializza i checkbox
        initializeColumnSelection(); // Re-inizializza i checkbox per riflettere lo stato corrente
        columnSelectionPanel.style.display = 'block';
        overlay.style.display = 'block';
    }
}

// --- Event Listeners ---
// Collega gli eventi ai rispettivi elementi DOM.
// Ho aggiunto controlli 'if (elemento)' per maggiore robustezza nel caso in cui un elemento non venga trovato.

if (searchInput) {
    searchInput.addEventListener('input', filterTable); // 'input' è più reattivo di 'keyup'
} else {
    console.error("Elemento 'searchInput' non trovato. La ricerca non funzionerà.");
}

if (applyColumnsBtn) {
    applyColumnsBtn.addEventListener('click', applySelectedColumns);
} else {
    console.error("Elemento 'applyColumnsBtn' non trovato.");
}

if (toggleColumnsBtn) {
    toggleColumnsBtn.addEventListener('click', (e) => {
        e.preventDefault(); // Previene il comportamento di default del link (es. scroll a #)
        toggleColumnSelectionPanel();
    });
} else {
    console.error("Elemento 'toggleColumnsBtn' non trovato.");
}

// Aggiungi un overlay per chiudere il pannello cliccando fuori.
// Crea l'overlay se non esiste già nel DOM.
let overlay = document.getElementById('overlay');
if (!overlay) {
    overlay = document.createElement('div');
    overlay.id = 'overlay';
    document.body.appendChild(overlay);
}
overlay.addEventListener('click', toggleColumnSelectionPanel);


// --- INIZIALIZZAZIONE COMPLETA DELLA PAGINA ---
// Esegui la logica iniziale solo quando il DOM è completamente caricato.
document.addEventListener('DOMContentLoaded', () => {
    // A questo punto, ALL_DATA e DEFAULT_VISIBLE_COLUMNS (iniettati nell'HTML) sono disponibili.
    // Inizializza le variabili di stato con i valori globali.
    // Se queste variabili non fossero state inizializzate sopra, andrebbero inizializzate qui.
    // Dato che sono già inizializzate nel global scope del JS, queste righe sono più per chiarezza/reset.
    currentData = [...ALL_DATA];
    visibleColumns = DEFAULT_VISIBLE_COLUMNS || [];

    // Chiamata iniziale per popolare la tabella con lo stato predefinito (nessun filtro, ordinamento predefinito).
    // filterTable() gestirà l'applicazione del filtro di ricerca iniziale (se presente),
    // l'ordinamento predefinito e la popolare della tabella.
    filterTable();
});